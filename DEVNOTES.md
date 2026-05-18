# DEVNOTES — bampipe 开发日记 (v11.4)

| | bampipe v11.3 @24 | samtools @24 | sambamba @24 |
|---|---|---|---|
| 墙钟 | **161s** | 289s | 319s¹ |
| CPU | 7.9c | 6.1c | — |
| 内存 | 38.4 GB | 33.2 GB | ~2 GB |
| 重复一致率 | ~100% | 基线 | ~100% |
| 输出 | 2.9 GB | 2.9 GB | 6.1 GB |
| 含 fixmate | ✓ | ✓ | ✗ |

¹ sambamba markdup(148s)+sort(171s)=319s (3-run avg). markdup alone 不排坐标序.

## 版本
| 版本 | 时间 | 关键 | 重复一致率 |
|------|:----:|------|:---------:|
| v10 | 122s | 并行SAM解析 | 30% |
| v11 | 158s | -C默认+双缓冲 | 30% |
| v11.2 | 195s | HT markdup+soft-clip | 97.4% |
| v11.3 | 161s | HT最適化 | 97.4% |
