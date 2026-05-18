#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
HTS_DIR="$SCRIPT_DIR/samtools-1.23.1"

# ---- 1. 检查 GCC 版本 (需要 ≥ 11 以支持 LTO) ----
GCC_VER=$(gcc -dumpversion | cut -d. -f1)
if [ "$GCC_VER" -lt 11 ]; then
    echo "错误: GCC 版本过低 (当前: $(gcc --version | head -1))"
    echo "  需要 GCC ≥ 11 (bampipe 使用 -flto 链接时优化)"
    exit 1
fi
echo "[check] GCC $(gcc --version | head -1) ✓"

# ---- 2. 检查 htslib 依赖 ----
if [ ! -d "$HTS_DIR/include" ] || [ ! -d "$HTS_DIR/lib" ]; then
    echo "错误: 依赖目录 samtools-1.23.1 不存在，请先编译 htslib-1.23.1"
    echo "  cd samtools-1.23.1 && make"
    exit 1
fi

# ---- 3. 编译 (忽略 htslib unused-result warning) ----
echo "[build] 编译 bampipe v11.4 ..."
gcc -O3 -march=native -flto -funroll-loops -finline-functions \
    -fno-semantic-interposition \
    -w \
    -I"$HTS_DIR/include" \
    -L"$HTS_DIR/lib" \
    -o "$SCRIPT_DIR/bampipe" "$SCRIPT_DIR/bampipe.c" \
    -lhts -ldeflate -lpthread -lz -lm \
    -Wl,-rpath,"$HTS_DIR/lib"

# ---- 4. 验证编译产物 ----
if [ -x "$SCRIPT_DIR/bampipe" ]; then
    echo "[build] 完成 → $SCRIPT_DIR/bampipe"
    echo "[build] $(file "$SCRIPT_DIR/bampipe")"
    echo ""
    $SCRIPT_DIR/bampipe -h
    echo ""
else
    echo "[build] 失败: 未生成可执行文件"
    exit 1
fi
