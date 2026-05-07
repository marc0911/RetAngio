#!/bin/bash
set -euo pipefail

REPO="/work/PRTNR/CHUV/HOJG/mschwab2/retangio/repo"
NAS_BASE="/nas/PRTNR/CHUV/HOJG/mschwab2/retangio/LTS/repo_snapshots"
STAMP="$(date +%F_%H%M%S)"
DEST="${NAS_BASE}/${STAMP}"

mkdir -p "${DEST}"

rsync -avh \
  --exclude 'data/extracted' \
  --exclude 'data/raw_archives' \
  --exclude '.snakemake' \
  --exclude '*.h5' \
  "${REPO}/" "${DEST}/"

echo
echo "Snapshot saved to:"
echo "${DEST}"
echo "${DEST}" > "${NAS_BASE}/LATEST_PATH.txt"
echo
echo "Latest snapshot written to:"
echo "${NAS_BASE}/LATEST_PATH.txt"
