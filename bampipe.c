/*
 * bampipe.c v11.5 — 高性能 BAM fixmate+sort+markdup 单管道 (并行SAM解析)
 *
 * === 架构 ===
 *   Phase 1a: Reader  — stdin SAM → hts_getline → 4线程并行sam_parse1 → QNAME配对 → SPSC ring buf
 *   Phase 1b: Fixmate — fixmate(MQ/MC/ms/ISIZE) → 填充sort缓冲 → HANDOFF
 *   Phase 1c: Sort    — LSD基数排序(idx-permute, 8-pass 8-bit) → temp BAM + .idx
 *   Phase 2:  Merge   — per-chr k-way heap merge+markdup (idx预读共享, 句柄复用)
 *   Phase 3:  Concat  — BGZF raw block 拼接 → 最终输出
 *
 * === v11.5 fix (2026-05-18) ==============================================
 *   Bug: Ru-2.sam 最长 CIGAR=52→MC标签≤56B, 加 MQ(7)+ms(7)=70B,
 *   BAM_HEADROOM=64 不足以容纳, bam_aux_append 调用 realloc(b->data)
 *   而 data 指向 slab 内部, 导致 "realloc(): invalid pointer" 崩溃.
 *   修复: BAM_HEADROOM 从 64 扩大到 128 (安全裕量), 零 malloc 开销.
 * ========================================================================
 * === v11.3 fix ===========================================================
 *   Bug: 主线程用 ready==3 判断批次完成, 但无代际区分. 新 dispatch 后 worker
 *   若快速完成, ready=3 被主线程误判为旧批次已完成, 导致跳过整个批次 (8192条).
 *   修复: 每个 worker 使用 task_id / task_ack 握手机制:
 *     dispatch: 主线程等 task_ack==0 后, task_id++, 设 task_ack=1, 发信号
 *     worker:   等 task_ack==1 后, 清 task_ack=0, 取 task_id, 解析
 *     complete: worker 设 task_done = my_id, 主线程等 task_done >= pending_id
 *   握手保证主线程不会覆盖未确认的分配, 消除代际跳过的可能.
 *   流水线: dispatch-then-wait (先派发后等待), 工人解析新批次与主线程处理旧结果重叠.
 * ========================================================================
 *   Slab arena: 64MB slab块分配bam1_t, 64B headroom防realloc, 消除80M次malloc/free
 *   LSD基数排序: 64-bit key→8-pass 8-bit LSD对idx排序→permute, O(n)
 *   unmapped预计算: srec_t缓存tid<0标志, srec_cmp/hent_lt避免访问bam1_t
 *   -C BGZF MT压缩: 临时BAM用wb+tpool, merge读端+tpool解压, per-chr输出+tpool
 *   文件句柄复用 + idx缓存 + .idx主线程预读 + FH round-robin + buf pool
 *
 * === 性能 (80M reads, 12chr, Xeon 80核, SSD, htslib-1.23.1) ===
 *   samtools 3-pipe @8 : 4m05s
 *   bampipe v9 @24 -C  : 3m12s (Phase1=146s, Merge=42s)
 *   比 samtools 快 53s (22%)
 *   注意: v11.3 的 gen 计数器保持 dispatch-then-wait, VLA 拷贝开销 < 流水线收益.
 *
 * === 适用规模 ===
 *   人基因组(1-22,X,Y):24chr ✓ | 800M reads: MAX_TMP_FILES=256 ✓
 *   1000+ chr组装: .idx架构与chr数无关 ✓
 *
 * === 编译 ===
 *   gcc -O3 -march=native -flto -funroll-loops -finline-functions \
 *       -fno-semantic-interposition \
 *       -I./samtools-1.23.1/include -L./samtools-1.23.1/lib \
 *       -o bampipe bampipe.c -lhts -ldeflate -lpthread -lz -lm \
 *       -Wl,-rpath,./samtools-1.23.1/lib
 */

/* ---- 头文件 ---- */
#include "htslib/sam.h"          // BAM/SAM 读写
#include "htslib/hts.h"          // htsFile, hts_set_opt
#include "htslib/bgzf.h"         // bgzf_tell/seek/raw_read/write
#include "htslib/kstring.h"      // kstring (构建 MC 标签)
#include "htslib/hts_os.h"       // hts_pos_t
#include "htslib/thread_pool.h"  // hts_tpool (BGZF 多线程压缩)
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <time.h>
#include <inttypes.h>
#include <pthread.h>

/* ---- 常量 ---- */
#define DEFAULT_MEM_MB   4096    // 默认每个线程内存 (4G)
#define MAX_TMP_FILES    256     // 最大临时 BAM 文件数
#define MD_MIN_QUALITY   15      // markdup 最低质量分数阈值
#define FIXMATE_RB_SIZE  2048    // reader→fixmate ring buffer 容量
#define RAWCPY_BUF_SZ    (16*1024*1024) // Phase 3 BGZF 拼接缓冲区
#define SLAB_SIZE       (64*1024*1024)  // 64MB slabs
#define BAM_HEADROOM     128             // MC(≤56)+MQ(7)+ms(7)=≤70B, 128B安全裕量防realloc

/* slab 分配器: 连续大块分配 bam1_t+data, 无逐条malloc */
typedef struct slab { struct slab *next; uint8_t *base, *cur, *end; } slab_t;
static slab_t *arena_head;
static void *arena_alloc(size_t sz) {
    sz = (sz + 7) & ~7;
    slab_t *s = arena_head;
    if (!s || s->cur + sz > s->end) {
        s = (slab_t*)malloc(sizeof(slab_t) + SLAB_SIZE);
        s->base = s->cur = (uint8_t*)(s + 1); s->end = s->base + SLAB_SIZE; s->next = arena_head; arena_head = s;
    }
    void *p = s->cur; s->cur += sz; return p;
}
static void arena_free_all(void) { slab_t *s = arena_head; while (s) { slab_t *n = s->next; free(s); s = n; } arena_head = NULL; }

/* 从 arena 分配 bam1_t + data 连续块, 不经过 malloc, 且数据区有 BAM_HEADROOM 裕量 */
static bam1_t *arena_dup_bam(bam1_t *src) {
    size_t sz = sizeof(bam1_t) + src->l_data + BAM_HEADROOM;
    bam1_t *b = (bam1_t*)arena_alloc(sz);
    *b = *src;
    b->data = (uint8_t*)b + sizeof(bam1_t);
    b->m_data = src->l_data + BAM_HEADROOM;
    memcpy(b->data, src->data, src->l_data);
    return b;
}

/* ================================================================
 * 1. 排序记录 & 比较函数
 * ================================================================ */
typedef struct { bam1_t *b; hts_pos_t pos; int32_t tid; uint8_t rev, unmapped, in_arena; } srec_t;

/* 坐标排序比较：先 tid，再 pos，再 rev；unmapped 排最后 */
static inline int srec_cmp(const srec_t *a, const srec_t *b) {
    if (a->unmapped != b->unmapped) return a->unmapped - b->unmapped;
    if (a->unmapped) return strcmp(bam_get_qname(a->b), bam_get_qname(b->b));
    if (a->tid != b->tid) return (a->tid < b->tid) ? -1 : 1;
    if (a->pos != b->pos) return (a->pos < b->pos) ? -1 : 1;
    if (a->rev != b->rev) return (a->rev < b->rev) ? -1 : 1;
    return 0;
}
static int srec_qcmp(const void *a, const void *b) { return srec_cmp((const srec_t*)a, (const srec_t*)b); }
static void srec_sort(srec_t *r, size_t n) { if (n > 1) qsort(r, n, sizeof(srec_t), srec_qcmp); }

/* LSD基数排序: 对idx数组排序(不触碰srec_t字段), 然后permute */
static void radix_sort_srec(srec_t *r, size_t n, int n_chr) {
    if (n <= 1) return;
    srec_t *orig=r; srec_t *tmp_r=malloc(n*sizeof(srec_t));
    /* 分离 unmapped */
    size_t lo=0, hi=n-1;
    while(lo<=hi){while(lo<n&&!orig[lo].unmapped)lo++;while(hi>0&&orig[hi].unmapped)hi--;if(lo<hi){srec_t t=orig[lo];orig[lo]=orig[hi];orig[hi]=t;lo++;hi--;}}
    size_t nm=lo;
    if (nm > 1) {
        /* 建 64-bit key + idx array */
        uint64_t *keys=malloc(nm*sizeof(uint64_t));
        uint32_t *idx=malloc(nm*sizeof(uint32_t));
        for(size_t i=0;i<nm;i++){keys[i]=((uint64_t)(orig[i].tid+1)<<40)|((uint64_t)(orig[i].pos+1)<<1)|orig[i].rev; idx[i]=i;}
        /* LSD sort on idx (8 passes, 8-bit) */
        uint32_t *tmp_idx=malloc(nm*sizeof(uint32_t));
        int cnt[256];
        for(int pass=0;pass<8;pass++){
            memset(cnt,0,sizeof(cnt));
            for(size_t i=0;i<nm;i++) cnt[(keys[idx[i]]>>(pass*8))&0xFF]++;
            size_t total=0; for(int j=0;j<256;j++){size_t c=cnt[j];cnt[j]=total;total+=c;}
            for(size_t i=0;i<nm;i++) tmp_idx[cnt[(keys[idx[i]]>>(pass*8))&0xFF]++]=idx[i];
            uint32_t *sw=idx;idx=tmp_idx;tmp_idx=sw;
        }
        /* swap back if needed, then permute */
        if(idx!=tmp_idx){free(tmp_idx);tmp_idx=idx;}else{free(idx);idx=tmp_idx;}
        for(size_t i=0;i<nm;i++) tmp_r[i]=orig[idx[i]];
        memcpy(orig,tmp_r,nm*sizeof(srec_t));
        free(idx); free(keys);
    }
    if(n-nm>1) qsort(orig+nm,n-nm,sizeof(srec_t),srec_qcmp);
    free(tmp_r); (void)n_chr;
}

/* ================================================================
 * 2. 最小堆 (k-way merge 用)
 * ================================================================ */
typedef struct { srec_t rec; int fi; } hent_t;  // fi = 来源文件索引

static inline int hent_lt(const hent_t *a, const hent_t *b) {
    if (a->rec.unmapped != b->rec.unmapped) return a->rec.unmapped < b->rec.unmapped;
    if (a->rec.unmapped) return strcmp(bam_get_qname(a->rec.b), bam_get_qname(b->rec.b)) < 0;
    if (a->rec.tid != b->rec.tid) return a->rec.tid < b->rec.tid;
    if (a->rec.pos != b->rec.pos) return a->rec.pos < b->rec.pos;
    if (a->rec.rev != b->rec.rev) return a->rec.rev < b->rec.rev;
    return 0;
}

typedef struct { hent_t *d; int sz, cap; } mheap_t;
static void mh_init(mheap_t *h, int c) { h->d = malloc(c * sizeof(hent_t)); h->sz = 0; h->cap = c; }
static void mh_free(mheap_t *h) { free(h->d); h->sz = h->cap = 0; }

static void mh_push(mheap_t *h, hent_t e) {
    if (h->sz >= h->cap) { h->cap *= 2; h->d = realloc(h->d, h->cap * sizeof(hent_t)); }
    int i = h->sz++; h->d[i] = e;
    while (i > 0) { int p = (i - 1) / 2; if (hent_lt(&h->d[i], &h->d[p])) { hent_t t = h->d[i]; h->d[i] = h->d[p]; h->d[p] = t; i = p; } else break; }
}
static hent_t mh_pop(mheap_t *h) {
    hent_t r = h->d[0]; h->d[0] = h->d[--h->sz]; int i = 0;
    while (1) { int l = 2*i+1, r = 2*i+2, s = i; if (l < h->sz && hent_lt(&h->d[l], &h->d[s])) s = l; if (r < h->sz && hent_lt(&h->d[r], &h->d[s])) s = r;
        if (s != i) { hent_t t = h->d[i]; h->d[i] = h->d[s]; h->d[s] = t; i = s; } else break; }
    return r;
}

/* ================================================================
 * 3. 工具函数: 质量分数、unclip 坐标、CIGAR 字符串化
 * ================================================================ */
static int64_t calc_score(bam1_t *b) { int64_t s = 0; uint8_t *q = bam_get_qual(b); for (int i = 0; i < b->core.l_qseq; i++) if (q[i] >= MD_MIN_QUALITY) s += q[i]; return s; }
static int64_t unclip_start(bam1_t *b) { if (b->core.n_cigar > 0) { uint32_t *c = bam_get_cigar(b); int op = bam_cigar_op(c[0]); if (op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP) return (int64_t)b->core.pos - bam_cigar_oplen(c[0]); } return b->core.pos; }
static int64_t unclip_end(bam1_t *b) { if (b->core.n_cigar > 0) { uint32_t *c = bam_get_cigar(b); int op = bam_cigar_op(c[b->core.n_cigar - 1]); hts_pos_t e = bam_endpos(b); if (op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP) e += bam_cigar_oplen(c[b->core.n_cigar - 1]); return e; } return bam_endpos(b); }
static int cigar2str(const bam1_t *b, kstring_t *s) { if (b->core.n_cigar == 0) return (kputc('*', s) == EOF) ? -1 : 0; const uint32_t *c = bam_get_cigar(b); for (int i = 0; i < b->core.n_cigar; i++) { if (kputw(bam_cigar_oplen(c[i]), s) == EOF) return -1; if (kputc(bam_cigar_opchr(c[i]), s) == EOF) return -1; } return 0; }
static size_t parse_size(const char *s) { char *e; double v = strtod(s, &e); if (*e == 'G' || *e == 'g') return (size_t)(v * 1073741824ULL); if (*e == 'M' || *e == 'm') return (size_t)(v * 1048576ULL); if (*e == 'K' || *e == 'k') return (size_t)(v * 1024ULL); return (size_t)(v * 1048576ULL); }
static int64_t mc_us(int mpos, const char *mc) { if (!mc || mc[0] == '*' || mc[0] == 0) return mpos; const char *p = mc; while (*p >= '0' && *p <= '9') p++; if (*p == 'H' || *p == 'S') { int h = (int)strtol(mc, NULL, 10); return mpos - h; } return mpos; }
static int64_t mc_ue(int mpos, const char *mc) { if (!mc || mc[0] == '*' || mc[0] == 0) return mpos; int64_t e = mpos; const char *p = mc, *lo = NULL; int ll = 0; while (*p) { int ln = (int)strtol(p, (char**)&p, 10); if (!*p) break; char op = *p++; lo = (p - 1); ll = ln; if (op == 'M' || op == 'D' || op == 'N' || op == '=' || op == 'X') e += ln; } if (lo && (*lo == 'H' || *lo == 'S')) e += ll; return e; }

/* ================================================================
 * 4. markdup: compact 16-byte key (32b tid/mtid/pos/mpos + ori + lm)
 * ================================================================ */
typedef struct { int32_t tid, mtid, pos, mpos; int8_t ori; } mdk_t;
static inline int mdk_eq(mdk_t a, mdk_t b) { return a.tid==b.tid && a.mtid==b.mtid && a.pos==b.pos && a.mpos==b.mpos && a.ori==b.ori; }
static inline uint64_t mdk_hash(mdk_t k) { return (uint64_t)k.tid ^ ((uint64_t)k.pos<<11) ^ ((uint64_t)k.mtid<<22) ^ ((uint64_t)k.mpos<<33) ^ ((uint64_t)k.ori<<44); }

static int64_t unclip_start_bam(bam1_t *b) {
    uint32_t *c = bam_get_cigar(b);
    int64_t cl = 0;
    for (int i = 0; i < b->core.n_cigar; i++) {
        int op = bam_cigar_op(c[i]);
        if (op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP) cl += bam_cigar_oplen(c[i]); else break;
    }
    return b->core.pos - cl;  /* bampipe compatible */
}
static int64_t unclip_end_bam(bam1_t *b) {
    int64_t e = bam_endpos(b);
    uint32_t *c = bam_get_cigar(b);
    for (int i = b->core.n_cigar - 1; i >= 0; i--) {
        int op = bam_cigar_op(c[i]);
        if (op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP) e += bam_cigar_oplen(c[i]); else break;
    }
    return e;
}
static int64_t five_prime(bam1_t *b) { return (b->core.flag & BAM_FREVERSE) ? unclip_end_bam(b) : unclip_start_bam(b); }
static int64_t mate_five_prime(bam1_t *b) { 
    if (b->core.flag & BAM_FMUNMAP) return b->core.mpos;
    uint8_t *d = bam_aux_get(b, "MC"); 
    if (!d) return b->core.mpos;
    const char *mc = bam_aux2Z(d);
    if (!mc || mc[0]=='*') return b->core.mpos;
    return (b->core.flag & BAM_FMREVERSE) ? mc_ue(b->core.mpos, mc) : mc_us(b->core.mpos, mc);
}

static mdk_t make_pair_key(bam1_t *b) {
    mdk_t k; int32_t tt=b->core.tid, mt=b->core.mtid;
    int64_t p1=five_prime(b), p2=mate_five_prime(b);
    int rev=(b->core.flag&BAM_FREVERSE)?1:0, mrev=(b->core.flag&BAM_FMREVERSE)?1:0, l;
    if (tt!=mt) l=(tt<mt); else if (p1!=p2) l=(p1<p2); else l=(rev<mrev);
    int lr=l?rev:mrev, rr=l?mrev:rev;
    k.tid=tt+1; k.mtid=mt+1; k.pos=(int32_t)(l?p1:p2); k.mpos=(int32_t)(l?p2:p1); k.ori=(lr<<1)|rr;
    return k;
}

#define PBCAP 64  // 位置缓冲初始容量
typedef struct { bam1_t **b; int64_t *sc; int n, cap, best; mdk_t key; } pbuf_t;
static void pb_init(pbuf_t *p) { p->cap = PBCAP; p->b = malloc(p->cap * sizeof(bam1_t*)); p->sc = malloc(p->cap * sizeof(int64_t)); p->n = 0; p->best = -1; p->key.tid = -2; }
static void pb_clear(pbuf_t *p) { for (int i = 0; i < p->n; i++) bam_destroy1(p->b[i]); p->n = 0; p->best = -1; p->key.tid = -2; }
static void pb_free(pbuf_t *p) { pb_clear(p); free(p->b); free(p->sc); }
static void pb_add(pbuf_t *p, bam1_t *b) { if (p->n >= p->cap) { p->cap *= 2; p->b = realloc(p->b, p->cap * sizeof(bam1_t*)); p->sc = realloc(p->sc, p->cap * sizeof(int64_t)); } p->b[p->n] = b; p->sc[p->n] = calc_score(b); if (p->best < 0 || p->sc[p->n] > p->sc[p->best]) p->best = p->n; p->n++; }
static void pb_flush(pbuf_t *p, samFile *o, sam_hdr_t *h, int rd, int64_t *w) { for (int i = 0; i < p->n; i++) { if (i != p->best && !rd) p->b[i]->core.flag |= BAM_FDUP; if (!(rd && (p->b[i]->core.flag & BAM_FDUP))) { sam_write1(o, h, p->b[i]); (*w)++; } } pb_clear(p); }

/* ================================================================
 * 5. Fixmate 操作: 同位配对、标签计算、ISIZE
 * ================================================================ */
static void f_unm(bam1_t *s,  bam1_t *d) { if ((d->core.flag & BAM_FUNMAP) && !(s->core.flag & BAM_FUNMAP)) { d->core.tid = s->core.tid; d->core.pos = s->core.pos; } }
static void f_mate(bam1_t *s, bam1_t *d) {
    d->core.mtid = s->core.tid; d->core.mpos = s->core.pos;
    if (s->core.flag & BAM_FREVERSE) d->core.flag |= BAM_FMREVERSE; else d->core.flag &= ~BAM_FMREVERSE;
    if (s->core.flag & BAM_FUNMAP) d->core.flag |= BAM_FMUNMAP;
}
static void f_aux_replace(bam1_t *d, const char tag[2], char type, int len, uint8_t *data) {
    uint8_t *t = bam_aux_get(d, tag);
    if (t) bam_aux_del(d, t);
    bam_aux_append(d, tag, type, len, data);
}
static int  f_mc(bam1_t *s, bam1_t *d)   {
    char buf[128]; kstring_t x = {0, 128, buf};                      // 栈预分配, 避免 realloc
    if (cigar2str(s, &x) < 0) return -1;
    { uint8_t *t = bam_aux_get(d, "MC"); if (t) bam_aux_del(d, t); }
    int r = bam_aux_append(d, "MC", 'Z', ks_len(&x) + 1, (uint8_t*)ks_str(&x));
    if (x.s != buf) free(x.s);
    return (r < 0) ? -1 : 0;
}
static void f_mq(bam1_t *s, bam1_t *d)   { uint32_t m = s->core.qual; f_aux_replace(d, "MQ", 'i', sizeof(m), (uint8_t*)&m); }
static void f_ms(bam1_t *s, bam1_t *d)   { int64_t sc = calc_score(s); uint32_t s32 = (uint32_t)sc; f_aux_replace(d, "ms", 'i', sizeof(s32), (uint8_t*)&s32); }
static void f_isize(bam1_t *a, bam1_t *b) {
    if (a->core.tid != b->core.tid || (a->core.flag & (BAM_FUNMAP | BAM_FMUNMAP)) || (b->core.flag & (BAM_FUNMAP | BAM_FMUNMAP))) { a->core.isize = b->core.isize = 0; return; }
    hts_pos_t a5 = (a->core.flag & BAM_FREVERSE) ? bam_endpos(a) : a->core.pos;
    hts_pos_t b5 = (b->core.flag & BAM_FREVERSE) ? bam_endpos(b) : b->core.pos;
    a->core.isize = (int32_t)(b5 - a5); b->core.isize = (int32_t)(a5 - b5);
}
static int  fixpair(bam1_t *a, bam1_t *b)   { f_unm(a,b); f_unm(b,a); f_mate(a,b); f_mate(b,a); f_isize(a,b); f_mq(a,b); f_mq(b,a); if (f_mc(a,b) < 0 || f_mc(b,a) < 0) return -1; f_ms(a,b); f_ms(b,a); return 0; }
static int  fixpair_light(bam1_t *a, bam1_t *b) { f_unm(a,b); f_unm(b,a); f_mate(a,b); f_mate(b,a); f_isize(a,b); if (f_mc(a,b) < 0 || f_mc(b,a) < 0) return -1; return 0; }
static void fixsingle(bam1_t *b) { if (!(b->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY))) { b->core.mtid = -1; b->core.mpos = -1; b->core.isize = 0; b->core.flag &= ~(BAM_FMREVERSE | BAM_FMUNMAP); } }

/* ================================================================
 * 6. Temp BAM 写出 + per-chr 索引 (.idx)
 *    wb0 = 未压缩 (默认), -C 开启 wb + MT 压缩/解压
 *    8MB BGZF 块减少系统调用
 * ================================================================ */
typedef struct { int32_t chr; int64_t offset; int64_t count; } idx_ent_t;
static int flush_buf_idx(srec_t *r, size_t n, const char *pf, int kidx, sam_hdr_t *h, char **out_bam, htsThreadPool *pool) {
    if (n == 0) return 0;
    int n_chr = sam_hdr_nref(h);
    radix_sort_srec(r, n, n_chr);                                    // 1. LSD基数排序 (idx permute)
    int64_t *cnt = calloc(n_chr + 1, sizeof(int64_t));
    for (size_t i = 0; i < n; i++) { int c = r[i].tid; if (r[i].unmapped || c < 0 || c >= n_chr) cnt[n_chr]++; else cnt[c]++; }
    int n_present = 0; for (int c = 0; c <= n_chr; c++) if (cnt[c] > 0) n_present++;
    char *bam_fn = malloc(strlen(pf) + 64);
    char *idx_fn = malloc(strlen(pf) + 68);
    { int bl = sprintf(bam_fn, "%s.%06d.bam", pf, kidx); sprintf(idx_fn, "%s.%06d.bam.idx", pf, kidx); (void)bl; }
    idx_ent_t *idx = malloc(n_present * sizeof(idx_ent_t)); int ei = 0; size_t ri = 0;
    int do_compress = (pool && pool->pool);                         // -C 标志启用时压缩
    htsFile *fp = hts_open(bam_fn, do_compress ? "wb" : "wb0");    // 3. 写入
    if (!fp) { free(bam_fn); free(idx_fn); free(cnt); free(idx); return -1; }
    hts_set_opt(fp, HTS_OPT_BLOCK_SIZE, 8 * 1024 * 1024);
    if (do_compress) hts_set_opt(fp, HTS_OPT_THREAD_POOL, pool);   // BGZF 多线程压缩
    if (sam_hdr_write(fp, h) < 0) { hts_close(fp); unlink(bam_fn); free(bam_fn); free(idx_fn); free(cnt); free(idx); return -1; }
    for (int c = 0; c <= n_chr; c++) {
        if (cnt[c] == 0) continue;
        if (do_compress) bgzf_flush(fp->fp.bgzf);                  // 确保 MT 压缩完成, bgzf_tell 准确
        idx[ei].chr = c; idx[ei].count = cnt[c]; idx[ei].offset = bgzf_tell(fp->fp.bgzf);
        for (int64_t j = 0; j < cnt[c]; j++) { if (sam_write1(fp, h, r[ri].b) < 0) { hts_close(fp); unlink(bam_fn); free(bam_fn); free(idx_fn); free(cnt); free(idx); return -1; } ri++; }
        ei++;
    }
    hts_close(fp); free(cnt);
    FILE *ifp = fopen(idx_fn, "wb");
    if (ifp) { int32_t ne = n_present; fwrite(&ne, sizeof(int32_t), 1, ifp); fwrite(idx, sizeof(idx_ent_t), n_present, ifp); fclose(ifp); }
    free(idx); free(idx_fn); *out_bam = bam_fn; return 0;
}

/* ================================================================
 * 7. SPSC 无锁 ring buffer (reader → fixmate)
 * ================================================================ */
typedef struct { bam1_t *a, *b; int is_supsec, in_arena; } fpair_t;
typedef struct { fpair_t *ents; volatile int head, tail; int size, mask, done; } fm_ring_t;
static void fmrb_init(fm_ring_t *r, int sz) { int s = 1; while (s < sz) s *= 2; r->ents = malloc(s * sizeof(fpair_t)); r->head = r->tail = 0; r->size = s; r->mask = s - 1; r->done = 0; }
static void fmrb_free(fm_ring_t *r) { free(r->ents); }
static inline int fmrb_full(fm_ring_t *r) { return (r->tail - r->head) >= r->size; }
static inline int fmrb_empty(fm_ring_t *r) { return r->head == r->tail; }
static void fmrb_put(fm_ring_t *r, fpair_t e) { while (fmrb_full(r)) { __sync_synchronize(); } int t = r->tail & r->mask; r->ents[t] = e; __sync_synchronize(); r->tail++; }
static int  fmrb_get(fm_ring_t *r, fpair_t *e) { while (fmrb_empty(r) && !r->done) { __sync_synchronize(); } if (fmrb_empty(r) && r->done) return 0; int h = r->head & r->mask; *e = r->ents[h]; __sync_synchronize(); r->head++; return 1; }
static void fmrb_finish(fm_ring_t *r) { r->done = 1; }

/* ================================================================
 * 8. Merge 线程: 按染色体并行 k-way heap 合并 + markdup
 *    - 动态染色体分配 (atomic fetch_add) 自动负载均衡
 *    - v4: 临时 BAM 每线程只 open 一次，跨染色体复用句柄
 *          每个染色体仅 bgzf_seek 到 .idx 记录的 byte 偏移
 *    - BGZF 线程池 (hts_tpool) 加速 per-chr 输出压缩
 * ================================================================ */
typedef struct { char **bam_fns; int n_files, n_chr; sam_hdr_t *hdr; char *pfx; int rd; volatile int *p_next_chr; char **chr_out; int64_t *chr_wr; pthread_mutex_t *out_mtx; htsThreadPool *pool; volatile int *p_merge_go; pthread_mutex_t *p_merge_mtx; pthread_cond_t *p_merge_cond; idx_ent_t **idx; int *idx_n; int **chr_off; } mthr_ctx_t;

static void *merge_thr_fn(void *arg) {
    mthr_ctx_t *ctx = (mthr_ctx_t*)arg; int n_chr = ctx->n_chr;

    // 等待 Phase1 完成信号
    pthread_mutex_lock(ctx->p_merge_mtx);
    while (!*(ctx->p_merge_go)) pthread_cond_wait(ctx->p_merge_cond, ctx->p_merge_mtx);
    pthread_mutex_unlock(ctx->p_merge_mtx);

    idx_ent_t **idx = ctx->idx; int *idx_n = ctx->idx_n; int **chr_off = ctx->chr_off;

    htsFile **chr_fp = calloc(ctx->n_files, sizeof(htsFile*));
    bam1_t **reuse_b = calloc(ctx->n_files, sizeof(bam1_t*));  /* bam1_t pool for merge loop */
    for (int i = 0; i < ctx->n_files; i++) {
        chr_fp[i] = hts_open(ctx->bam_fns[i], "rb"); if (!chr_fp[i]) continue;
        if (ctx->pool && ctx->pool->pool) hts_set_opt(chr_fp[i], HTS_OPT_THREAD_POOL, ctx->pool); // MT 解压
        sam_hdr_t *th = sam_hdr_read(chr_fp[i]); if (th) sam_hdr_destroy(th);
    }
    
    /* per-thread hash table for markdup (pair keys + seen_count) */
#define HT_SZ 8388608
#define HT_MASK (HT_SZ-1)
    mdk_t *htk = calloc(HT_SZ, sizeof(mdk_t));
    int8_t *htcnt = calloc(HT_SZ, 1);  /* seen_f1 bits */
    int *htu = malloc(HT_SZ * sizeof(int)); int htu_n = 0;
    for (int i=0; i<HT_SZ; i++) htk[i].tid = -2;

    while (1) {
        int c = __sync_fetch_and_add(ctx->p_next_chr, 1);
        if (c > n_chr) break;
        int is_unmapped = (c == n_chr);
        int nfc = 0; for (int i = 0; i < ctx->n_files; i++) if (chr_off[i][c] >= 0) nfc++;
        if (nfc == 0) continue;
        char *chr_fn = malloc(strlen(ctx->pfx) + 64); sprintf(chr_fn, "%s._ch%04d.bam", ctx->pfx, c);
        htsFile *chr_out = hts_open(chr_fn, "wb"); if (!chr_out) { free(chr_fn); continue; }
        if (sam_hdr_write(chr_out, ctx->hdr) < 0) { hts_close(chr_out); unlink(chr_fn); free(chr_fn); continue; }
        hts_set_opt(chr_out, HTS_OPT_BLOCK_SIZE, 8 * 1024 * 1024);
        if (ctx->pool) hts_set_opt(chr_out, HTS_OPT_THREAD_POOL, ctx->pool);

        mheap_t mh; mh_init(&mh, nfc > 0 ? nfc : 64);               // 预分配避免 realloc
        for (int i = 0; i < ctx->n_files; i++) {
            int j = chr_off[i][c]; if (j < 0 || !chr_fp[i]) continue;
            bgzf_seek(chr_fp[i]->fp.bgzf, idx[i][j].offset, SEEK_SET);
            reuse_b[i] = bam_init1();
            if (sam_read1(chr_fp[i], ctx->hdr, reuse_b[i]) >= 0) { hent_t e; e.rec.tid = reuse_b[i]->core.tid; e.rec.pos = reuse_b[i]->core.pos; e.rec.rev = (reuse_b[i]->core.flag & BAM_FREVERSE) ? 1 : 0; e.rec.unmapped = (reuse_b[i]->core.tid < 0) ? 1 : 0; e.rec.b = reuse_b[i]; e.fi = i; mh_push(&mh, e); }
        }
        if (is_unmapped) {                                          // Unmapped: 直接写出（不做 markdup）
            int64_t uwr = 0;
            while (mh.sz > 0) { hent_t t = mh_pop(&mh); sam_write1(chr_out, ctx->hdr, t.rec.b); uwr++; int fi = t.fi;
                bam1_t *nb = t.rec.b; int more = 0; if (chr_fp[fi] && sam_read1(chr_fp[fi], ctx->hdr, nb) >= 0 && ((nb->core.tid < 0 && c == n_chr) || nb->core.tid == c)) { hent_t ne; ne.rec.tid = is_unmapped ? (nb->core.tid < 0 ? n_chr : nb->core.tid) : c; ne.rec.pos = nb->core.pos; ne.rec.rev = (nb->core.flag & BAM_FREVERSE) ? 1 : 0; ne.rec.unmapped = (nb->core.tid < 0) ? 1 : 0; ne.rec.b = nb; ne.fi = fi; mh_push(&mh, ne); } }
            hts_close(chr_out); pthread_mutex_lock(ctx->out_mtx); ctx->chr_wr[c] = uwr; ctx->chr_out[c] = chr_fn; pthread_mutex_unlock(ctx->out_mtx);
        } else {                                                    // Mapped: k-way merge + markdup (samtools pair-key, no f1 in key)
            int64_t wr = 0;
            for (int i=0; i<htu_n; i++) { int idx=htu[i]; htk[idx].tid=-2; htcnt[idx]=0; }
            htu_n = 0;

            while (mh.sz > 0) { hent_t t = mh_pop(&mh); bam1_t *b = t.rec.b; int fi = t.fi;
                if (t.rec.unmapped) { sam_write1(chr_out, ctx->hdr, b); wr++; }
                else {
                    mdk_t k; uint8_t *mc = bam_aux_get(b, "MC");
                    if (mc && b->core.mtid >= 0 && !(b->core.flag & BAM_FMUNMAP)) k = make_pair_key(b);
                    else { k.tid=b->core.tid+1; k.mtid=0;
                        k.pos=(int32_t)((b->core.flag&BAM_FREVERSE)?unclip_end_bam(b):unclip_start_bam(b));
                        k.mpos=(int32_t)b->core.mpos; k.ori=((b->core.flag&BAM_FREVERSE)?2:0)|((b->core.flag&BAM_FMREVERSE)?1:0); }
                    uint32_t idx = (uint32_t)((uint64_t)k.tid ^ ((uint64_t)k.pos<<4) ^ ((uint64_t)k.mtid<<12) ^ ((uint64_t)k.mpos<<20) ^ ((uint64_t)k.ori<<30)) & HT_MASK;
                    int coll = 0;
                    while (htk[idx].tid != -2 && !mdk_eq(htk[idx], k) && ++coll < 128) idx = (idx + 1) & HT_MASK;
                    if (coll < 128 && htu_n < HT_MASK) {
                        if (htk[idx].tid == -2) {
                            htk[idx] = k; htu[htu_n++] = idx;
                            htcnt[idx] = (b->core.flag & BAM_FREAD1) ? 2 : 1;
                        } else {
                            int f1_bit = (b->core.flag & BAM_FREAD1) ? 2 : 1;
                            if (htcnt[idx] & f1_bit) {
                                if (!ctx->rd) b->core.flag |= BAM_FDUP;
                            } else { htcnt[idx] |= f1_bit; }
                        }
                    }
                    sam_write1(chr_out, ctx->hdr, b); wr++;
                }
                if (chr_fp[fi] && sam_read1(chr_fp[fi], ctx->hdr, b) >= 0 && ((b->core.tid < 0 && c == n_chr) || b->core.tid == c)) { hent_t ne; ne.rec.tid = b->core.tid; ne.rec.pos = b->core.pos; ne.rec.rev = (b->core.flag & BAM_FREVERSE) ? 1 : 0; ne.rec.unmapped = (b->core.tid < 0) ? 1 : 0; ne.rec.b = b; ne.fi = fi; mh_push(&mh, ne); } }
            hts_close(chr_out);
            pthread_mutex_lock(ctx->out_mtx); ctx->chr_wr[c] = wr; ctx->chr_out[c] = chr_fn; pthread_mutex_unlock(ctx->out_mtx);
        }
        mh_free(&mh);
    }
    free(htk); free(htcnt); free(htu);
    for (int i = 0; i < ctx->n_files; i++) { if (reuse_b[i]) bam_destroy1(reuse_b[i]); if (chr_fp[i]) hts_close(chr_fp[i]); }
    free(reuse_b); free(chr_fp);
    return NULL;  // idx/idx_n/chr_off shared, freed by main
}

/* ================================================================
 * 9. Sort worker 线程: 从 fixmate 获取排序后缓冲区，排序并写临时 BAM
 * ================================================================ */
typedef struct { srec_t *buf; size_t n; int idx; char *pfx; sam_hdr_t *hdr; htsThreadPool *pool; char *out_fn; srec_t *reuse_buf; int ready, quit; pthread_mutex_t mtx; pthread_cond_t cond; } swork_t;
static void *swork_fn(void *arg) { swork_t *w = (swork_t*)arg;
    while (1) { pthread_mutex_lock(&w->mtx); while (w->ready != 1 && !w->quit) pthread_cond_wait(&w->cond, &w->mtx); if (w->quit && w->ready != 1) { pthread_mutex_unlock(&w->mtx); break; }
        srec_t *buf = w->buf; size_t n = w->n; int idx = w->idx; char *pfx = w->pfx; sam_hdr_t *hdr = w->hdr; htsThreadPool *pool = w->pool; w->ready = 2; pthread_mutex_unlock(&w->mtx);
        char *fn = NULL; flush_buf_idx(buf, n, pfx, idx, hdr, &fn, pool);
        for (size_t j = 0; j < n; j++) if (!buf[j].in_arena) bam_destroy1(buf[j].b);
        pthread_mutex_lock(&w->mtx); w->out_fn = fn; w->reuse_buf = buf; w->ready = 3; pthread_cond_signal(&w->cond); pthread_mutex_unlock(&w->mtx); }
    return NULL;
}

/* ================================================================
 * 10. Fixmate 线程: 从 ring buffer 取配对，执行 fixmate
 *     v8: buf pool 复用避免频繁 malloc/free, cigar2str 栈预分配
 * ================================================================ */
typedef struct { fm_ring_t *rb; swork_t *sw; int nw; char **tfn; int *p_ntmp, *p_next_idx; sam_hdr_t *hdr; char *pfx; size_t mrec; int64_t *fixed; int no_aux; } fixmate_ctx_t;
static void *fixmate_thr_fn(void *arg) {
    fixmate_ctx_t *fctx = (fixmate_ctx_t*)arg; fm_ring_t *rb = fctx->rb; swork_t *sw = fctx->sw; int nw = fctx->nw; char **tfn = fctx->tfn; int *ntmp = fctx->p_ntmp;
    srec_t *buf_pool[8] = {NULL}; int pool_n = 0; int last_w = 0;
    srec_t *buf = malloc(fctx->mrec * sizeof(srec_t)); size_t nbuf = 0;
#define FH { int w; for(w=last_w; w<nw; w++) { pthread_mutex_lock(&sw[w].mtx); if(sw[w].ready!=1&&sw[w].ready!=2) { if(sw[w].ready==3) { char *fn=sw[w].out_fn; srec_t *rbuf=sw[w].reuse_buf; sw[w].ready=0; pthread_mutex_unlock(&sw[w].mtx); if(fn) tfn[(*ntmp)++]=fn; if(rbuf && pool_n<8) buf_pool[pool_n++]=rbuf; } else pthread_mutex_unlock(&sw[w].mtx); break; } pthread_mutex_unlock(&sw[w].mtx); } if(w==nw) { for(w=0; w<last_w; w++) { pthread_mutex_lock(&sw[w].mtx); if(sw[w].ready!=1&&sw[w].ready!=2) { if(sw[w].ready==3) { char *fn=sw[w].out_fn; srec_t *rbuf=sw[w].reuse_buf; sw[w].ready=0; pthread_mutex_unlock(&sw[w].mtx); if(fn) tfn[(*ntmp)++]=fn; if(rbuf && pool_n<8) buf_pool[pool_n++]=rbuf; } else pthread_mutex_unlock(&sw[w].mtx); break; } pthread_mutex_unlock(&sw[w].mtx); } } if(w==last_w||w==nw) { w=0; pthread_mutex_lock(&sw[0].mtx); while(sw[0].ready==1||sw[0].ready==2) pthread_cond_wait(&sw[0].cond, &sw[0].mtx); if(sw[0].ready==3) { char *fn=sw[0].out_fn; srec_t *rbuf=sw[0].reuse_buf; sw[0].ready=0; pthread_mutex_unlock(&sw[0].mtx); if(fn) tfn[(*ntmp)++]=fn; if(rbuf && pool_n<8) buf_pool[pool_n++]=rbuf; } else pthread_mutex_unlock(&sw[0].mtx); } last_w=(w+1)%nw; srec_t *old_buf=buf; size_t old_n=nbuf; if(pool_n>0) { buf=buf_pool[--pool_n]; } else { buf=malloc(fctx->mrec*sizeof(srec_t)); } nbuf=0; pthread_mutex_lock(&sw[w].mtx); sw[w].buf=old_buf; sw[w].n=old_n; sw[w].idx=(*fctx->p_next_idx)++; sw[w].out_fn=NULL; sw[w].reuse_buf=NULL; sw[w].ready=1; pthread_cond_signal(&sw[w].cond); pthread_mutex_unlock(&sw[w].mtx); }
    fpair_t pr; while (fmrb_get(rb, &pr))
        if (pr.is_supsec)    { fixsingle(pr.a); if (nbuf >= fctx->mrec) { FH } buf[nbuf].tid = pr.a->core.tid; buf[nbuf].pos = pr.a->core.pos; buf[nbuf].rev = (pr.a->core.flag & BAM_FREVERSE) ? 1 : 0; buf[nbuf].unmapped = (pr.a->core.tid < 0) ? 1 : 0; buf[nbuf].in_arena = pr.in_arena; buf[nbuf].b = pr.a; nbuf++; }
        else if (!pr.b)      { fixsingle(pr.a); if (nbuf >= fctx->mrec) { FH } buf[nbuf].tid = pr.a->core.tid; buf[nbuf].pos = pr.a->core.pos; buf[nbuf].rev = (pr.a->core.flag & BAM_FREVERSE) ? 1 : 0; buf[nbuf].unmapped = (pr.a->core.tid < 0) ? 1 : 0; buf[nbuf].in_arena = pr.in_arena; buf[nbuf].b = pr.a; nbuf++; }
        else                 { int r = fctx->no_aux ? fixpair_light(pr.a, pr.b) : fixpair(pr.a, pr.b); (*fctx->fixed) += 2; if (r < 0) { continue; } if (nbuf + 2 > fctx->mrec) { FH } buf[nbuf].tid = pr.a->core.tid; buf[nbuf].pos = pr.a->core.pos; buf[nbuf].rev = (pr.a->core.flag & BAM_FREVERSE) ? 1 : 0; buf[nbuf].unmapped = (pr.a->core.tid < 0) ? 1 : 0; buf[nbuf].in_arena = pr.in_arena; buf[nbuf].b = pr.a; nbuf++; buf[nbuf].tid = pr.b->core.tid; buf[nbuf].pos = pr.b->core.pos; buf[nbuf].rev = (pr.b->core.flag & BAM_FREVERSE) ? 1 : 0; buf[nbuf].unmapped = (pr.b->core.tid < 0) ? 1 : 0; buf[nbuf].in_arena = pr.in_arena; buf[nbuf].b = pr.b; nbuf++; }
    if (nbuf > 0) { FH }
#undef FH
    return NULL;
}

/* ================================================================
 * 11. Parallel SAM parser workers (task_id / task_ack handshake)
 *     sam_parse1 is thread-safe (reads header, writes bam1_t)
 *     dispatch: main thread waits for task_ack==0, then increments
 *               task_id, sets task_ack=1, signals worker.
 *     worker:   wakes when task_ack==1, clears task_ack=0, reads task_id,
 *               parses assigned chunk, sets task_done = my_id.
 *     wait:     main thread waits until task_done >= pending_id[wi],
 *               guaranteeing this specific task is finished.
 *     handshake prevents overwriting unacknowledged assignments,
 *     eliminating any possibility of worker skipping a generation.
 * ================================================================ */
#define PARSE_BATCH 8192
#define MAX_LINE_SZ 512

typedef struct {
    char **lines;          // Share of SAM lines
    bam1_t **results;      // Share of result bam1_t pointers
    int n_recs;            // Records for this worker
    sam_hdr_t *hdr;        // Shared header
    int ready;             // 0=idle, 1=assigned, 2=working, 3=done
    int task_id;           // Task number for current dispatch
    int task_ack;          // 0=not acknowledged, 1=pending ack
    int task_done;         // Highest completed task_id
    int quit;
    pthread_mutex_t mtx;
    pthread_cond_t cond;
} pworker_t;

static void *parser_thr_fn(void *arg) {
    pworker_t *pw = (pworker_t*)arg;
    while (1) {
        pthread_mutex_lock(&pw->mtx);
        while (!pw->task_ack && !pw->quit)
            pthread_cond_wait(&pw->cond, &pw->mtx);
        if (pw->quit) { pthread_mutex_unlock(&pw->mtx); break; }
        pw->task_ack = 0;
        int my_id = pw->task_id;
        int n = pw->n_recs;
        char **lines = pw->lines;
        bam1_t **res = pw->results;
        sam_hdr_t *hdr = pw->hdr;
        pw->ready = 2;
        pthread_mutex_unlock(&pw->mtx);

        for (int i = 0; i < n; i++) {
            kstring_t ks = {0, 0, NULL};
            ks.s = lines[i];
            ks.l = strlen(lines[i]);
            ks.m = ks.l + 1;
            res[i] = bam_init1();
            if (sam_parse1(&ks, hdr, res[i]) < 0) { bam_destroy1(res[i]); res[i] = NULL; }
        }
        pthread_mutex_lock(&pw->mtx);
        pw->task_done = my_id;
        pw->ready = 3;
        pthread_cond_signal(&pw->cond);
        pthread_mutex_unlock(&pw->mtx);
    }
    return NULL;
}

/* ================================================================
 * 12. 主函数
 *     流程: Reader → Fixmate → Sort Pool → Merge Pool → Concat
 *     输入: 必须为 name-sorted SAM (bwa mem 默认输出即 name-sorted)
 *     输出: coordinate-sorted, fixmated, mark-duplicated BAM
 *
 *     参数说明:
 *       -@ INT     线程数 (N个sort+N个merge+fixmate+reader, 默认8)
 *                  建议设为 CPU 核心数
 *       -m SIZE    每个 sort 线程最大内存 (默认4G, 支持K/M/G后缀)
 *                  总内存约 = -m × (-@ 的一半)
 *       -T PREFIX  临时文件前缀 (默认 /tmp/bampipe_tmp)
 *                  建议指向 SSD 或 tmpfs
 *       -o FILE    输出 BAM 文件 (默认 stdout)
 *       -r         物理移除重复 (默认仅标记 BAM_FDUP flag)
 *       --no-PG    不添加 @PG header 行
 *       -h         显示帮助
 *
 *     使用示例:
 *       bwa mem ref.fa R1.fq R2.fq | bampipe -@ 24 -m 4G -o out.bam
 *       cat name_sorted.sam | bampipe -@ 12 -m 2G -T /ssd/tmp -o out.bam
 * ================================================================ */
int main(int argc, char **argv) {
    int n_thr = 8, compress_temp = 1, no_aux = 0; size_t max_m = (size_t)DEFAULT_MEM_MB * 1048576ULL;
    char *pf = NULL, *ofn = NULL, *ifn = NULL; int rd = 0, npg = 0, explicit_T = 0;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-@") && i + 1 < argc) n_thr = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-m") && i + 1 < argc) max_m = parse_size(argv[++i]);
        else if (!strcmp(argv[i], "-T") && i + 1 < argc) { free(pf); pf = strdup(argv[++i]); explicit_T = 1; }
        else if (!strcmp(argv[i], "-o") && i + 1 < argc) ofn = argv[++i];
        else if (!strcmp(argv[i], "-r")) rd = 1;
        else if (!strcmp(argv[i], "-i") && i + 1 < argc) ifn = argv[++i];
        else if (!strcmp(argv[i], "-C")) compress_temp = 1;
        else if (!strcmp(argv[i], "-A")) no_aux = 1;
        else if (!strcmp(argv[i], "--no-PG")) npg = 1;
        else if (!strcmp(argv[i], "-h")) { fprintf(stderr, "bampipe v11.5 — fixmate+sort+markdup + parallel SAM parsing + BGZF default\nUsage: bampipe [options] < in.sam > out.bam\n -i FILE  Input BAM/SAM file (default: stdin)\n -@ INT   Threads [8]\n -m SIZE  Memory per thread [4G]\n -T PREFIX  Temp prefix (default: same dir as -o, or ./bampipe_tmp)\n -o FILE  Output BAM\n -r       Remove duplicates\n -A       Skip MQ/ms tags (faster)\n --no-PG  No @PG header\n -h       This help\nNote: BGZF compression is enabled by default.\n"); free(pf); return 0; }
    }
    if (!pf) pf = strdup("bampipe_tmp");
    if (!explicit_T && ofn && strcmp(ofn, "-") && strcmp(ofn, "/dev/stdout") && strcmp(ofn, "/dev/null")) {
        char *slash = strrchr(ofn, '/');
        if (slash) {
            int dlen = (int)(slash - ofn);
            char *dir = malloc(dlen + 16);
            memcpy(dir, ofn, dlen); dir[dlen] = '\0';
            free(pf); pf = malloc(dlen + 32);
            sprintf(pf, "%s/bampipe_tmp", dir);
            free(dir);
        }
    }
    if (n_thr < 1) n_thr = 1; if (max_m < 1048576) max_m = 1048576;
    samFile *in = ifn ? hts_open(ifn, "r") : hts_open("-", "r");
    if (!in) { free(pf); return 1; }
    sam_hdr_t *hdr = sam_hdr_read(in); if (!hdr) { hts_close(in); free(pf); return 1; }
    sam_hdr_remove_line_id(hdr, "HD", "SO", NULL); sam_hdr_add_line(hdr, "HD", "VN", "1.6", "SO", "coordinate", NULL);
    if (!npg) sam_hdr_add_pg(hdr, "bampipe", "ID", "bampipe", "PN", "bampipe", "VN", "10.0", NULL);

    /* ---- Phase 2 预备 + 共享 BGZF 线程池 ---- */
    int nw = n_thr; if (nw < 2) nw = 2; if (nw > 48) nw = 48;
    int n_chr = sam_hdr_nref(hdr); int nm = nw; if (nm < 2) nm = 2; if (nm > n_chr + 1) nm = n_chr + 1; if (nm > 32) nm = 32;

    hts_tpool *tpool = NULL; htsThreadPool tp = {NULL, 0};
    if (nm > 1 || compress_temp) {
        int psz = compress_temp ? nm : (nm / 2);   // 压缩模式用满, 非压缩用半
        if (psz < 2) psz = 2; if (psz > 24) psz = 24;
        tpool = hts_tpool_init(psz);
        if (tpool) { tp.pool = tpool; tp.qsize = psz * 2; }
    }
    htsThreadPool *shared_tp = tpool ? &tp : NULL;
    htsThreadPool *compress_pool = compress_temp ? shared_tp : NULL;  // -C: Phase1 压缩用
    htsThreadPool *merge_pool = shared_tp;                            // Phase2/3 始终用

    volatile int merge_go = 0;
    pthread_mutex_t merge_mtx = PTHREAD_MUTEX_INITIALIZER;
    pthread_cond_t  merge_cond = PTHREAD_COND_INITIALIZER;

    mthr_ctx_t *mctx = calloc(nm, sizeof(mthr_ctx_t)); pthread_t *mt = malloc(nm * sizeof(pthread_t));
    for (int t = 0; t < nm; t++) { mctx[t].n_chr = n_chr; mctx[t].hdr = hdr; mctx[t].pfx = pf; mctx[t].rd = rd;
        mctx[t].p_merge_go = &merge_go; mctx[t].p_merge_mtx = &merge_mtx; mctx[t].p_merge_cond = &merge_cond;
        pthread_create(&mt[t], NULL, merge_thr_fn, &mctx[t]); }

    /* ---- Phase 1: Fixmate + Sort ---- */
    size_t mrec = max_m / nw / sizeof(srec_t); if (mrec > 30000000) mrec = 30000000; if (mrec < 200000) mrec = 200000;
    int ntmp = 0, next_idx = 0; char **tfn = malloc(MAX_TMP_FILES * sizeof(char*)); int64_t tin = 0, fixed = 0; time_t t0 = time(NULL);

    swork_t *sw = calloc(nw, sizeof(swork_t)); pthread_t *swt = malloc(nw * sizeof(pthread_t));
    for (int w = 0; w < nw; w++) { sw[w].pfx = pf; sw[w].hdr = hdr; sw[w].pool = compress_pool; pthread_mutex_init(&sw[w].mtx, NULL); pthread_cond_init(&sw[w].cond, NULL); pthread_create(&swt[w], NULL, swork_fn, &sw[w]); }

    fm_ring_t fmrb; fmrb_init(&fmrb, FIXMATE_RB_SIZE);
    fixmate_ctx_t fctx; fctx.rb = &fmrb; fctx.sw = sw; fctx.nw = nw; fctx.tfn = tfn; fctx.p_ntmp = &ntmp; fctx.p_next_idx = &next_idx; fctx.hdr = hdr; fctx.pfx = pf; fctx.mrec = mrec; fctx.fixed = &fixed; fctx.no_aux = no_aux;
    pthread_t fm_thr; pthread_create(&fm_thr, NULL, fixmate_thr_fn, &fctx);

    /* ---- Parallel SAM reader with batch-based parsing ---- */
    int np = nw / 4; if (np < 1) np = 1; if (np > 4) np = 4;
    pworker_t *pwork = calloc(np, sizeof(pworker_t));
    pthread_t *pthrs = malloc(np * sizeof(pthread_t));
    for (int i = 0; i < np; i++) { pwork[i].hdr = hdr;
        pthread_mutex_init(&pwork[i].mtx, NULL); pthread_cond_init(&pwork[i].cond, NULL);
        pthread_create(&pthrs[i], NULL, parser_thr_fn, &pwork[i]); }

    /* double-buffered batch reader: fill one batch while parsers work on the other */
#define NBUF 2
    typedef struct { char buf[PARSE_BATCH * MAX_LINE_SZ]; char *lines[PARSE_BATCH];
             bam1_t *bams[PARSE_BATCH]; int n, used; } batch_t;
    batch_t *bbuf = calloc(NBUF, sizeof(batch_t));
    int cur = 0, pending = -1, n_lines = 0, buf_used = 0; bam1_t *prev = NULL;
    int pending_id[np]; for (int i = 0; i < np; i++) pending_id[i] = 0;
    kstring_t ks = {0, 0, NULL};
    while (hts_getline(in, '\n', &ks) >= 0) {
        if (ks.l == 0 || ks.s[0] == '@' || ks.s[0] == '\n') continue;
        int len = (int)ks.l;
        if (len > 0 && ks.s[len-1] == '\n') len--;
        if (len > 0 && ks.s[len-1] == '\r') len--;
        /* process previous batch result before starting to fill */
        if (n_lines >= PARSE_BATCH || buf_used + len + 1 > (int)sizeof(bbuf[0].buf)) {
            int disp = cur; cur = 1 - cur;
            bbuf[disp].n = n_lines; bbuf[disp].used = buf_used;
            /* dispatch new batch FIRST — workers start parsing while main reads/processes */
            int per_w = n_lines / np;
            int wait_id[np];
            for (int wi = 0; wi < np; wi++) {
                int st = wi * per_w, en = (wi == np-1) ? n_lines : (wi+1) * per_w;
                pthread_mutex_lock(&pwork[wi].mtx);
                while (pwork[wi].task_ack)
                    pthread_cond_wait(&pwork[wi].cond, &pwork[wi].mtx);
                pwork[wi].task_id++;
                wait_id[wi] = pwork[wi].task_id;
                pwork[wi].task_ack = 1;
                pwork[wi].lines = bbuf[disp].lines + st;
                pwork[wi].results = bbuf[disp].bams + st;
                pwork[wi].n_recs = en - st;
                pwork[wi].ready = 1;
                pthread_cond_signal(&pwork[wi].cond);
                pthread_mutex_unlock(&pwork[wi].mtx);
            }
            /* then wait for previously dispatched batch (parallelism: workers already busy on disp) */
            if (pending >= 0) {
                for (int wi = 0; wi < np; wi++) {
                    pthread_mutex_lock(&pwork[wi].mtx);
                    while (pwork[wi].task_done < pending_id[wi])
                        pthread_cond_wait(&pwork[wi].cond, &pwork[wi].mtx);
                    pthread_mutex_unlock(&pwork[wi].mtx);
                }
                for (int i = 0; i < bbuf[pending].n; i++, tin++) {
                bam1_t *r = bbuf[pending].bams[i]; if (!r) { continue; }
                int sup = (r->core.flag & BAM_FSUPPLEMENTARY), sec = (r->core.flag & BAM_FSECONDARY);
                    if (sup || sec) { fpair_t pr; pr.a = arena_dup_bam(r); pr.b = NULL; pr.is_supsec = 1; pr.in_arena = 1; fmrb_put(&fmrb, pr); bam_destroy1(r); }
                    else if (!prev) { prev = arena_dup_bam(r); bam_destroy1(r); }
                    else if (strcmp(bam_get_qname(prev), bam_get_qname(r)) == 0) { fpair_t pr; pr.a = prev; pr.b = arena_dup_bam(r); pr.is_supsec = 0; pr.in_arena = 1; fmrb_put(&fmrb, pr); bam_destroy1(r); prev = NULL; }
                    else { fpair_t pr; pr.a = prev; pr.b = NULL; pr.is_supsec = 0; pr.in_arena = 1; fmrb_put(&fmrb, pr); prev = arena_dup_bam(r); bam_destroy1(r); }
                }
            }
            for (int wi = 0; wi < np; wi++) pending_id[wi] = wait_id[wi];
            pending = disp;  /* disp is now being parsed, results pending */
            n_lines = 0; buf_used = 0;
        }
        memcpy(bbuf[cur].buf + buf_used, ks.s, (size_t)len);
        bbuf[cur].buf[buf_used + len] = '\0';
        bbuf[cur].lines[n_lines] = bbuf[cur].buf + buf_used;
        buf_used += len + 1; n_lines++;
    }
    /* Flush final batch */
    if (n_lines > 0) {
        /* dispatch remaining lines first — workers start while main processes pending */
        bbuf[cur].n = n_lines; bbuf[cur].used = buf_used;
        int per_w = n_lines / np;
        int flush_id[np];
        for (int wi = 0; wi < np; wi++) {
            int st = wi * per_w, en = (wi == np-1) ? n_lines : (wi+1) * per_w;
            pthread_mutex_lock(&pwork[wi].mtx);
            while (pwork[wi].task_ack)
                pthread_cond_wait(&pwork[wi].cond, &pwork[wi].mtx);
            pwork[wi].task_id++;
            flush_id[wi] = pwork[wi].task_id;
            pwork[wi].task_ack = 1;
            pwork[wi].lines = bbuf[cur].lines + st;
            pwork[wi].results = bbuf[cur].bams + st;
            pwork[wi].n_recs = en - st;
            pwork[wi].ready = 1;
            pthread_cond_signal(&pwork[wi].cond);
            pthread_mutex_unlock(&pwork[wi].mtx);
        }
        /* then wait for previous pending batch */
        if (pending >= 0) {
            for (int wi = 0; wi < np; wi++) {
                pthread_mutex_lock(&pwork[wi].mtx);
                while (pwork[wi].task_done < pending_id[wi])
                    pthread_cond_wait(&pwork[wi].cond, &pwork[wi].mtx);
                pthread_mutex_unlock(&pwork[wi].mtx);
            }
            for (int i = 0; i < bbuf[pending].n; i++, tin++) {
                bam1_t *r = bbuf[pending].bams[i]; if (!r) { continue; }
                if (r->core.flag & (BAM_FSUPPLEMENTARY|BAM_FSECONDARY)) { fpair_t pr; pr.a = arena_dup_bam(r); pr.b = NULL; pr.is_supsec = 1; pr.in_arena = 1; fmrb_put(&fmrb, pr); bam_destroy1(r); }
                else if (!prev) { prev = arena_dup_bam(r); bam_destroy1(r); }
                else if (strcmp(bam_get_qname(prev), bam_get_qname(r)) == 0) { fpair_t pr; pr.a = prev; pr.b = arena_dup_bam(r); pr.is_supsec = 0; pr.in_arena = 1; fmrb_put(&fmrb, pr); bam_destroy1(r); prev = NULL; }
                else { fpair_t pr; pr.a = prev; pr.b = NULL; pr.is_supsec = 0; pr.in_arena = 1; fmrb_put(&fmrb, pr); prev = arena_dup_bam(r); bam_destroy1(r); }
            }
        }
        for (int wi = 0; wi < np; wi++) pending_id[wi] = flush_id[wi];
        pending = cur;
    }
    /* process last pending batch */
    if (pending >= 0) {
        for (int wi = 0; wi < np; wi++) {
            pthread_mutex_lock(&pwork[wi].mtx);
            while (pwork[wi].task_done < pending_id[wi])
                pthread_cond_wait(&pwork[wi].cond, &pwork[wi].mtx);
            pthread_mutex_unlock(&pwork[wi].mtx);
        }
            for (int i = 0; i < bbuf[pending].n; i++, tin++) {
                bam1_t *r = bbuf[pending].bams[i]; if (!r) { continue; }
                int sup = (r->core.flag & BAM_FSUPPLEMENTARY), sec = (r->core.flag & BAM_FSECONDARY);
            if (sup || sec) { fpair_t pr; pr.a = arena_dup_bam(r); pr.b = NULL; pr.is_supsec = 1; pr.in_arena = 1; fmrb_put(&fmrb, pr); bam_destroy1(r); }
            else if (!prev) { prev = arena_dup_bam(r); bam_destroy1(r); }
            else if (strcmp(bam_get_qname(prev), bam_get_qname(r)) == 0) { fpair_t pr; pr.a = prev; pr.b = arena_dup_bam(r); pr.is_supsec = 0; pr.in_arena = 1; fmrb_put(&fmrb, pr); bam_destroy1(r); prev = NULL; }
            else { fpair_t pr; pr.a = prev; pr.b = NULL; pr.is_supsec = 0; pr.in_arena = 1; fmrb_put(&fmrb, pr); prev = arena_dup_bam(r); bam_destroy1(r); }
        }
    }
    if (prev) { fpair_t pr; pr.a = prev; pr.b = NULL; pr.is_supsec = 0; pr.in_arena = 1; fmrb_put(&fmrb, pr); }
    /* Stop parser threads */
    for (int i = 0; i < np; i++) { pthread_mutex_lock(&pwork[i].mtx); pwork[i].quit = 1; pthread_cond_signal(&pwork[i].cond); pthread_mutex_unlock(&pwork[i].mtx); }
    for (int i = 0; i < np; i++) { pthread_join(pthrs[i], NULL); pthread_mutex_destroy(&pwork[i].mtx); pthread_cond_destroy(&pwork[i].cond); }
    free(pwork); free(pthrs);
    fmrb_finish(&fmrb); pthread_join(fm_thr, NULL);

    for (int w = 0; w < nw; w++) { pthread_mutex_lock(&sw[w].mtx); while (sw[w].ready == 1 || sw[w].ready == 2) pthread_cond_wait(&sw[w].cond, &sw[w].mtx); if (sw[w].ready == 3) { char *fn = sw[w].out_fn; sw[w].ready = 0; if (fn) tfn[ntmp++] = fn; } sw[w].quit = 1; pthread_cond_signal(&sw[w].cond); pthread_mutex_unlock(&sw[w].mtx); }
    for (int w = 0; w < nw; w++) { pthread_join(swt[w], NULL); pthread_mutex_destroy(&sw[w].mtx); pthread_cond_destroy(&sw[w].cond); }
    free(sw); free(swt); fmrb_free(&fmrb);
    time_t t1 = time(NULL); fprintf(stderr, "[bampipe] Phase1: %ld in, %ld fixmated, %d tmps (N=%d), %.0fs\n", (long)tin, (long)fixed, ntmp, nw, difftime(t1, t0));

    /* ---- Phase 2: merge 线程已就绪, 主线程预读 .idx 共享 ---- */
    idx_ent_t **shared_idx = calloc(ntmp, sizeof(idx_ent_t*));
    int *shared_idx_n = calloc(ntmp, sizeof(int));
    int **shared_chr_off = calloc(ntmp, sizeof(int*));
    for (int i = 0; i < ntmp; i++) {
        shared_chr_off[i] = malloc((n_chr + 1) * sizeof(int));
        for (int c = 0; c <= n_chr; c++) shared_chr_off[i][c] = -1;
        char ifn[4096]; sprintf(ifn, "%s.idx", tfn[i]);
        FILE *ifp = fopen(ifn, "rb"); if (!ifp) continue;
        int32_t ne; if (fread(&ne, sizeof(int32_t), 1, ifp) != 1) { fclose(ifp); continue; }
        shared_idx[i] = malloc(ne * sizeof(idx_ent_t));
        if (fread(shared_idx[i], sizeof(idx_ent_t), ne, ifp) != (size_t)ne) { fclose(ifp); continue; }
        shared_idx_n[i] = ne; fclose(ifp);
        for (int j = 0; j < ne; j++) if (shared_idx[i][j].chr <= n_chr) shared_chr_off[i][shared_idx[i][j].chr] = j;
    }
    char **chr_out = calloc(n_chr + 1, sizeof(char*)); int64_t *chr_wr = calloc(n_chr + 1, sizeof(int64_t));
    pthread_mutex_t out_mtx = PTHREAD_MUTEX_INITIALIZER; volatile int next_chr = 0;
    for (int t = 0; t < nm; t++) { mctx[t].bam_fns = tfn; mctx[t].n_files = ntmp;
        mctx[t].idx = shared_idx; mctx[t].idx_n = shared_idx_n; mctx[t].chr_off = shared_chr_off;
        mctx[t].pool = merge_pool; mctx[t].p_next_chr = &next_chr; mctx[t].chr_out = chr_out; mctx[t].chr_wr = chr_wr; mctx[t].out_mtx = &out_mtx; }
    merge_go = 1; pthread_cond_broadcast(&merge_cond);
    for (int t = 0; t < nm; t++) pthread_join(mt[t], NULL);
    free(mctx); free(mt); if (tpool) hts_tpool_destroy(tpool);
    for (int i = 0; i < ntmp; i++) { free(shared_idx[i]); free(shared_chr_off[i]); }
    free(shared_idx); free(shared_idx_n); free(shared_chr_off); time_t t2 = time(NULL);

    /* ---- Phase 3: BGZF block-level concat ---- */
    samFile *out = ofn ? hts_open(ofn, "wb") : hts_open("-", "wb");
    if (out) { hts_set_opt(out, HTS_OPT_BLOCK_SIZE, 4 * 1024 * 1024); if (n_thr > 1) hts_set_threads(out, n_thr); sam_hdr_write(out, hdr); bgzf_flush(out->fp.bgzf);
        int64_t total_wr = 0; for (int c = 0; c <= n_chr; c++) total_wr += chr_wr[c]; uint8_t *cpybuf = malloc(RAWCPY_BUF_SZ);
        for (int c = 0; c <= n_chr; c++) { if (!chr_out[c]) continue; samFile *cf = hts_open(chr_out[c], "rb"); if (!cf) { free(chr_out[c]); continue; } sam_hdr_t *th = sam_hdr_read(cf); if (th) sam_hdr_destroy(th); bgzf_flush(cf->fp.bgzf); ssize_t nr; while ((nr = bgzf_raw_read(cf->fp.bgzf, cpybuf, RAWCPY_BUF_SZ)) > 0) bgzf_raw_write(out->fp.bgzf, cpybuf, (size_t)nr); hts_close(cf); unlink(chr_out[c]); free(chr_out[c]); }
        free(cpybuf); hts_close(out); time_t t3 = time(NULL); fprintf(stderr, "[bampipe] merge: %ld recs, %.0fs, total: %.0fs\n", (long)total_wr, difftime(t3, t1), difftime(t3, t0)); }
    for (int i = 0; i < ntmp; i++) { char ifn[4096]; sprintf(ifn, "%s.idx", tfn[i]); unlink(ifn); unlink(tfn[i]); free(tfn[i]); }
    free(tfn); free(chr_out); free(chr_wr); sam_hdr_destroy(hdr); hts_close(in); free(pf);
    free(bbuf); arena_free_all(); return 0;
}
