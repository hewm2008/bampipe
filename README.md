# bampipe — 高性能 fixmate + sort + markdup 单管道 (v11)

[![Language](https://img.shields.io/badge/language-C-blue.svg)]()

将 `samtools fixmate`、`samtools sort`、`samtools markdup` 三步融合为**单次数据遍历**，
消除 3× parse + 3× serialize 的冗余。默认启用 BGZF 多线程压缩（`-C`）。

**bampipe @24 仅 2m38s**，比 samtools 三管道 @24（4m49s）快 **45%（-131s）**。

## 设计思路

```
samtools:  SAM → [fixmate] → BAM → [sort] → BAM → [markdup] → BAM
           3x parse + 3x serialize

bampipe:    SAM → [fixmate → sort → markdup] → BAM
           1x parse + 1x serialize
```

### 架构 (v11)

```
Phase 1a: Reader  线程  — stdin/file 读 SAM, 批量缓存, 4线程并行 sam_parse1 → QNAME配对 → ring buf
Phase 1b: Fixmate 线程  — 取配对 → fixmate → 填充排序缓冲 → HANDOFF
Phase 1c: Sort    N线程  — LSD基数排序(idx-permute) → 写 BGZF 压缩 temp BAM + .idx (-C 默认)
Phase 2:  Merge   M线程  — 按染色体并行 k-way heap merge + markdup → 写 per-chr BAM
Phase 3:  Concat        — BGZF 块级拼接 per-chr BAM → 最终输出
```

### 核心优化

| 优化点 | 说明 | 收益 |
|--------|------|------|
| 单次解析 | 1× parse + 1× serialize vs samtools 3× | CPU -66% |
| **并行 SAM 解析** | 4线程并行 `sam_parse1` | Phase1 大幅加速 |
| **-C 默认启用** | BGZF 多线程压缩 temp BAM，CPU 利用率 3.6→8.1 核 | Phase1 +45% |
| 双缓冲 reader | reader 填充下批时 parsers 跑上批 | 流水线重叠 |
| Fixmate 线程化 | reader 和 fixmate 分离 | 流水线重叠 |
| 染色体并行合并 | 各染色体分线程合并 | Merge 充分并行 |
| Slab arena | 64MB 大块 bump alloc, 消除 80M malloc/free | Phase1 -66s |
| LSD 基数排序 | 8-pass 8-bit idx-permute, O(n) 排序 | 与 qsort 持平 |
| BGZF 线程池压缩 | 共享 hts_tpool 加速输出 | ~25s |

## 性能

### 测试环境

- CPU: Xeon Gold 6133 80 核 @ 2.50 GHz
- 内存: 502 GB
- 磁盘: SSD (本地) + 网络存储
- 输入: BWA-mem paired-end SAM, 80,198,610 reads, 12 条染色体, 32 GB

### 对比

| 方案 | 墙钟 | Phase1 | Merge | CPU% | 有效核数 | 内存 |
|------|:----:|:------:|:-----:|:----:|:--------:|:----:|
| **bampipe @24 (-C 默认)** | **2m38s** | 119s | 39s | 807% | 8.1 | 36.4GB |
| samtools 3-pipe @24 | 4m49s | — | — | 610% | 6.1 | 33.2GB |
| **bampipe @16 (-C 默认)** | **2m51s** | 129s | 38s | 763% | 7.6 | 36.3GB |
| **bampipe @8 (-C 默认)** | **3m57s** | 174s | 60s | 570% | 5.7 | 36.7GB |
| sambamba markdup @24 (BAM) | 2m27s | — | — | 1182% | 11.8 | 2.0GB |

> bampipe 默认启用 BGZF 压缩（-C）。sambamba 仅做 markdup（不做 fixmate），输入为 6 GB BAM。

### 正确性

| 检查项 | 结果 |
|--------|:----:|
| 记录数 | 80,198,610 ✓ (与 samtools 完全一致) |
| Flag 分布 | 与 samtools 一致 |
| ISIZE/MC/MQ/ms 标签 | 与 samtools 完全一致 |
| BAM 索引 | `samtools index` 正常生成 ✓ |

## 编译 & 安装

### 依赖

- htslib 1.23.1（含 libdeflate）
- zlib
- pthread
- GCC ≥ 11（LTO 优化）

### 编译

```sh
bash make.sh
```

或手动：

```sh
gcc -O3 -march=native -flto -funroll-loops -finline-functions \
    -fno-semantic-interposition \
    -I./samtools-1.23.1/include -L./samtools-1.23.1/lib \
    -o bampipe bampipe.c -lhts -ldeflate -lpthread -lz -lm \
    -Wl,-rpath,./samtools-1.23.1/lib
```

## 使用

```sh
# pipe SAM from BWA-MEM (推荐)
bwa mem ref.fa R1.fq R2.fq | ./bampipe -@ 24 -m 4G -o out.bam

# 读取 SAM 文件
./bampipe -i input.sam -@ 24 -m 4G -o out.bam

# 读取 BAM 文件 (最快，跳过 SAM 文本解析)
./bampipe -i input.bam -@ 24 -m 4G -o out.bam

# 跳过 MQ/ms 标签 (更快)
cat input.sam | ./bampipe -@ 24 -m 4G -A -o output.bam
```

### 参数

| 参数 | 说明 | 默认 |
|------|------|:----:|
| `-@ INT` | 线程数 (sort + merge + reader) | 8 |
| `-m SIZE` | 每个 sort 线程最大内存 (总内存 ≈ `-m × -@`) | 4G |
| `-i FILE` | 输入文件 (SAM/BAM/CRAM, 自动检测格式) | stdin |
| `-o FILE` | 输出 BAM | stdout |
| `-T PREFIX` | 临时文件前缀 | /tmp/bampipe_tmp |
| `-r` | 物理移除重复 (默认仅标记 BAM_FDUP) | — |
| `-A` | 跳过 MQ/ms 标签 (更快, 不生成 aux 标签) | — |
| `--no-PG` | 不添加 @PG header | — |
| `-h` | 帮助 | — |

> `-C` 已默认启用（BGZF 多线程压缩）。无需手动传入。

## 版本演进

```
v3:    fixmate 线程分离 + BGZF tpool                → 4m16s
v6:    unmapped 预计算 + -C BGZF 压缩               → 3m51s
v8:    FH 轮询 + idx 共享 + buf pool + 编译优化      → 4m11s
v9:    Slab arena + LSD 基数排序                     → 3m12s
v9+:   -i FILE 文件输入 + BAM 输入                   → 2m19s
v10:   并行 SAM 解析 (sam_parse1 4线程) + in_arena   → 2m03s
v11:   -C 默认启用 + 双缓冲 reader                   → 2m38s①
```

① v11 性能数据源自网络存储上的 big_test.sam（32 GB），DEVNOTES v10 数据源自缓存热状态。

## 优缺点

### 优势

1. **单次解析** — 1× parse + 1× serialize vs samtools 3×
2. **-C 默认启用** — BGZF 多线程压缩，CPU 利用率 8.1 核
3. **快 45%** @24 — 比 samtools 三管道 @24 快 **45%（158s vs 289s）**
4. **代码紧凑** — 约 650 行 C，易于理解和修改
5. **单进程** — 无进程间通信开销
6. **灵活输入** — 支持 pipe SAM / 文件 SAM / 文件 BAM

### 限制

1. **输入必须是 name-sorted** — 不支持 coordinate-sorted 输入
2. **不支持 optical duplicate detection** — 仅位置驱动的去重
3. **需要 htslib-1.23.1** — 依赖特定版本的 BGZF API
4. **临时文件** — Phase 1 产生 BGZF 压缩临时 BAM（`-C` 默认，减少磁盘占用）

## 注意事项

1. `-m` 是**每个 sort 线程**的内存上限，总可用内存 ≈ `-m × -@`
2. 临时文件默认使用 BGZF 压缩（`-C` 默认）
3. 染色体数通过 BAM header 动态获取，支持任意参考基因组
4. 输出 BAM 的 DUP flag 与 samtools 一致
5. `-i -` 从 stdin 读取，格式由 htslib 自动检测（SAM/BAM/CRAM）

## 许可

MIT — 继承自 samtools/htslib.

## 讨论与联系

- :email: hewm2008@gmail.com / hewm2008@qq.com
- **QQ 群: 125293663**

