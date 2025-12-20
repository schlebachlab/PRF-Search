#!/usr/bin/env bash
set -euo pipefail

# Recombine split LFS parts into the original .tsv.gz
# Usage:
#   ./recombine_peptide2parent.sh
# or
#   ./recombine_peptide2parent.sh output.tsv.gz

OUT="${1:-not_trimmed_with_all_zf_v101o_fs_pep_all_msfrag_nme_always_peptide2parent.tsv.gz}"
PARTS=( "${OUT}.part-00" "${OUT}.part-01" "${OUT}.part-02" )

for p in "${PARTS[@]}"; do
  if [[ ! -f "$p" ]]; then
    echo "Missing part: $p" >&2
    exit 1
  fi
done

cat "${PARTS[@]}" > "$OUT"
echo "Wrote: $OUT"

# Optional integrity checks (if sha256sum is available)
if command -v sha256sum >/dev/null 2>&1; then
  echo "sha256: $(sha256sum "$OUT" | awk '{print $1}')  $OUT"
fi
