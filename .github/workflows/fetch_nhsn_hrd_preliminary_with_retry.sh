#!/usr/bin/env bash
set -euo pipefail

DATA_DIR="data/interim/cases/NHSN-HRD_archive/preliminary"
MAX_RETRIES="${MAX_RETRIES:-6}" # Use env vars if set, otherwise defaults
WAIT_SECONDS="${WAIT_SECONDS:-3600}"

echo "Starting NHSN HRD preliminary data fetch"
echo "Max retries: $MAX_RETRIES"
echo "Wait between retries: ${WAIT_SECONDS}s"

mkdir -p "$DATA_DIR"

for ((i=1; i<=MAX_RETRIES; i++)); do
  echo "=================================================="
  echo "Attempt $i of $MAX_RETRIES"

  python data/conversion/cases/fetch-format_NHSN-HRD-data.py --preliminary True

  NEWEST_FILE=$(ls -t "$DATA_DIR"/NHSN-HRD_reference-date-*_gathered-*.parquet.gzip | head -n 1)

  if [[ -z "$NEWEST_FILE" ]]; then
    echo "::error::No NHSN HRD files found after data fetch"
    exit 1
  fi

  echo "Newest file:"
  echo "  $NEWEST_FILE"

  REF_DATE=$(basename "$NEWEST_FILE" | sed -E 's/.*reference-date-([0-9]{4}-[0-9]{2}-[0-9]{2})_.*/\1/')
  echo "Extracted reference date: $REF_DATE"

  MATCHING_FILES=$(ls "$DATA_DIR"/NHSN-HRD_reference-date-"$REF_DATE"_gathered-*.parquet.gzip 2>/dev/null | wc -l)

  if [[ "$MATCHING_FILES" -eq 1 ]]; then
    echo "New reference date detected — continuing workflow."
    exit 0
  fi

  echo "Reference date already exists ($MATCHING_FILES files)."

  if [[ "$i" -lt "$MAX_RETRIES" ]]; then
    echo "Waiting $WAIT_SECONDS seconds before retry..."
    sleep "$WAIT_SECONDS"
  else
    echo "::error::No new NHSN HRD reference date after $MAX_RETRIES attempts"
    exit 1
  fi
done
