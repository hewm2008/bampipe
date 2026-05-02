# bampipe — 高性能 fixmate + sort + markdup 单管道程序

[![Language](https://img.shields.io/badge/language-C-blue.svg)]()

将 `samtools fixmate`、`samtools sort`、`samtools markdup` 三个步骤融合为**单次数据遍历**，
消除 3× parse + 3× serialize 的冗余，利用多线程并行显著缩短 wall-time。

在 8000 万条 reads 测试中：**bampipe @24 = 4m16s**，比 samtools 三管道 @8（4m38s）快 **22 秒**。

## 设计思路

```
samtools:  SAM → [fixmate] → BAM → [sort] → BAM → [markdup] → BAM
           3x parse + 3x serialize

bampipe:     SAM → [fixmate → sort → markdup] → BAM
           1x parse + 1x serialize
```

### 架构

```
Phase 1a: Reader  线程  — stdin 读 SAM，按 QNAME 配对 → ring buffer
Phase 1b: Fixmate 线程  — 取配对 → fixmate → 填充排序缓冲 → HANDOFF
Phase 1c: Sort    N线程  — qsort 缓冲区 → 写 wb0 临时 BAM + .idx 索引
Phase 2:  Merge   M线程  — 按染色体并行 k-way heap merge + markdup → 写 per-chr BAM
Phase 3:  Concat        — BGZF 块级拼接 per-chr BAM → 最终输出
```

### 核心优化

| 优化点 | 说明 | 收益 |
|--------|------|------|
| 单次解析 | 1× parse + 1× serialize vs samtools 的 3× | **~45% CPU** |
| Fixmate 线程化 | reader 和 fixmate 分离到两个线程，CPU 重叠 | **~80s** |
| 染色体并行合并 | 12 个染色体分别由不同线程合并 | **~30s** |
| 8MB BGZF 写入块 | 减少系统调用（samtools 用 2MB） | **~22s** |
| BGZF 线程池压缩 | 共享 hts_tpool 加速 per-chr 输出 | **~25s** |
| 消除双次 bam_dup1 | fixmate 中 read B 不再重复复制 | **~10s** |

## 性能

### 测试环境

- CPU: Xeon 64 核
- 内存: 512GB
- 磁盘: 本地 SSD
- 输入: BWA-mem paired-end SAM, 80,198,610 reads, 12 条染色体

### 对比

| 方案 | 线程数 | 时间 | vs samtools |
|------|:------:|:----:|:-----------:|
| samtools 3-pipe | @8 (24 有效) | 4m38s | — |
| **bampipe** | **@24** | **4m16s** | **-22s** |
| samtools 3-pipe | @4 (12 有效) | 7m02s | — |
| bampipe | @12 | 5m56s | -66s |

### 正确性

- 记录数: 80,198,610 ✓ (与 samtools 完全一致)
- Flag 差异: 55 / 199,972 = 0.027% (来自 DUP 标记的 tiebreaker 差异)
- ISIZE/MC/MQ/ms 标签: 与 samtools 完全一致
- BAM 索引: `samtools index` 正常生成

## 编译 & 安装

### 依赖

- htslib 1.23.1（含 libdeflate）
- zlib
- pthread
- GCC ≥ 11（LTO 优化）

```sh
# 编译
gcc -O3 -march=native -flto -funroll-loops -finline-functions \
    -I./samtools-1.23.1/include \
    -L./samtools-1.23.1/lib \
    -o bampipe bampipe.c \
    -lhts -ldeflate -lpthread -lz -lm \
    -Wl,-rpath,./samtools-1.23.1/lib
```

### 使用

```sh
# 基本用法（从 stdin 读取 name-sorted SAM，输出到 stdout）
cat input.sam | ./bampipe -@ 24 -m 4G -o output.bam

# 自定义临时目录（建议用 SSD 所在路径）
cat input.sam | ./bampipe -@ 24 -m 4G -T /ssd/tmp -o output.bam

# 移除重复（而非标记 DUP flag）
cat input.sam | ./bampipe -@ 24 -m 4G -r -o dedup.bam

# 直接接 BWA
bwa mem ref.fa R1.fq R2.fq | ./bampipe -@ 24 -m 8G -o aligned.sort.dup.bam
```

### 参数

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `-@ INT` | 线程数 | 8 |
| `-m SIZE` | 每个线程最大排序内存（如 4G, 500M） | 4G |
| `-T PREFIX` | 临时文件前缀 | /tmp/bampipe_tmp |
| `-o FILE` | 输出 BAM 文件 | stdout |
| `-r` | 移除重复（不标记 DUP flag） | 仅标记 |
| `--no-PG` | 不添加 @PG 头行 | — |
| `-h` | 帮助 | — |

> **输入必须是 name-sorted SAM**（`bwa mem` 的默认输出）。

### 验证

```sh
# 生成 samtools 参考输出
cat small_test.sam | samtools fixmate -m -@ 4 - - | \
    samtools sort -@ 4 -m 500M - | \
    samtools markdup -@ 4 -s - ref.bam

# 生成 bampipe 输出
cat small_test.sam | ./bampipe -@ 4 -m 2G -o test.bam

# 比较 flag（期望 ~55 条差异）
diff <(samtools view ref.bam | awk '{print $1"\t"$2}' | sort) \
     <(samtools view test.bam | awk '{print $1"\t"$2}' | sort) | wc -l
```

## 优缺点

### 优势

1. **单次解析** — 节省 3× 的解析开销
2. **4m16s @24** — 比 samtools 三管道更快
3. **代码紧凑** — 约 200 行 C，易于理解和修改
4. **单进程** — 无进程间通信开销
5. **内存灵活** — 总可用内存 = `-m × -@`，充分利用大内存机器

### 限制

1. **输入必须是 name-sorted** — 不支持 coordinate-sorted 输入
2. **不支持 optical duplicate detection** — 仅位置驱动的去重
3. **需要 htslib-1.23.1** — 依赖特定版本的 BGZF API
4. **临时文件** — Phase 1 产生未压缩临时 BAM（约 30GB / 8000 万 reads）
5. **无增量模式** — 必须完整处理所有数据

## 注意事项

1. `-m` 是**每个线程**的内存上限（与 samtools -m 语义一致），总可用内存 ≈ `-m × -@`
2. 临时文件使用 **wb0 模式**（未压缩），需要足够的磁盘空间
3. 染色体数通过 BAM header 动态获取，支持任意参考基因组
4. 输出 BAM 的 DUP flag 与 samtools 有 0.027% 的差异（tiebreaker），不影响下游分析

## 许可

MIT — 继承自 [samtools](https://github.com/samtools/samtools) 和 [htslib](https://github.com/samtools/htslib)。

## 讨论与联系

- :email: hewm2008@gmail.com / hewm2008@qq.com
- 加入 **QQ 群: 125293663**

