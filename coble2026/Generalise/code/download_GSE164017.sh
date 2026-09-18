#!/usr/bin/env bash
set -euo pipefail

DATA_DIR="data"
BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE164nnn/GSE164017/suppl"

mkdir -p "$DATA_DIR"

wget -c -O "$DATA_DIR/GSE164017_RAW.tar" "$BASE/GSE164017_RAW.tar"
wget -c -O "$DATA_DIR/GSE164017_features.tsv.gz" "$BASE/GSE164017_features.tsv.gz"

echo "done: $(ls -lh $DATA_DIR)"