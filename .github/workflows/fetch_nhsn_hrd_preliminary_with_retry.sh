#!/usr/bin/env bash
set -euo pipefail

DATA_DIR="data/interim/cases/NHSN-HRD_archive/preliminary"
MAX_RETRIES="${MAX_RETRIES:-6}" # Use env vars if set, otherwise defaults
WAIT_SECONDS="${WAIT_SECONDS:-3600}"
MAX_PY_FAILURES="${MAX_PY_FAILURES:-3}" # Give up after this many python script errors
PY_FAILURES=0

echo "Starting NHSN HRD preliminary data fetch"
echo "Max retries: $MAX_RETRIES"
echo "Max python script failures before giving up: $MAX_PY_FAILURES"
echo "Wait between retries: ${WAIT_SECONDS}s"

mkdir -p "$DATA_DIR"

for ((i=1; i<=MAX_RETRIES; i++)); do
  echo "=================================================="
  echo "Attempt $i of $MAX_RETRIES"
  echo "Downloading data"

  # run python script to fetch data --> if it fails, counts as a retry (intermittent errors have happened)
  if ! python data/conversion/cases/fetch-format_NHSN-HRD-data.py --preliminary True; then
    PY_FAILURES=$((PY_FAILURES + 1))
    echo "::warning::Python data-fetch script failed on attempt $i of $MAX_RETRIES (python failure $PY_FAILURES of $MAX_PY_FAILURES)."

    if [[ "$PY_FAILURES" -ge "$MAX_PY_FAILURES" ]]; then
      echo "::error::Python data-fetch script has failed $PY_FAILURES times — giving up instead of retrying further."
      exit 1
    fi

    if [[ "$i" -lt "$MAX_RETRIES" ]]; then
      echo "Waiting $WAIT_SECONDS seconds before retrying..."
      sleep "$WAIT_SECONDS"
      continue
    else
      echo "::error::Python data-fetch script failed after $MAX_RETRIES attempts."
      exit 1
    fi
  fi

  NEWEST_FILE=$(ls -t "$DATA_DIR"/NHSN-HRD_reference-date-*_gathered-*.parquet.gzip | head -n 1)

  if [[ -z "$NEWEST_FILE" ]]; then
    echo "::error::No NHSN HRD files found after data fetch"
    exit 1
  fi

  echo "Filename:"
  echo "  $NEWEST_FILE"

  REF_DATE=$(basename "$NEWEST_FILE" | sed -E 's/.*reference-date-([0-9]{4}-[0-9]{2}-[0-9]{2})_.*/\1/')
  echo "Reference date:"
  echo "  $REF_DATE"

  MATCHING_FILES=$(ls "$DATA_DIR"/NHSN-HRD_reference-date-"$REF_DATE"_gathered-*.parquet.gzip 2>/dev/null | wc -l)
  OTHER_FILES=$((MATCHING_FILES - 1))

  if [[ "$MATCHING_FILES" -eq 1 ]]; then
    echo "New reference date detected — continuing workflow."
    exit 0
  else
    echo "Data for reference date already exists ($OTHER_FILES other files were found)."
    echo "Removing:"
    echo "  $NEWEST_FILE"
    rm -f "$NEWEST_FILE"
  fi

  if [[ "$i" -lt "$MAX_RETRIES" ]]; then
    echo "Waiting $WAIT_SECONDS seconds before retrying..."
    sleep "$WAIT_SECONDS"
  else
    echo "::error::No new NHSN HRD reference date after $MAX_RETRIES attempts"
    exit 1
  fi
done
