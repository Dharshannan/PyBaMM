#!/usr/bin/env bash
# run.sh TAG COND [VAR=VALUE ...]
#   COND: 25 | 45 | 25p1595 | 25p2080  (rpt_soc_plots.py configs; always --aligned)
#   Extra VAR=VALUE pairs are passed as environment (e.g. C64_POS_LAM_MULT=8).
#   Outputs land in this folder; "TAG exit=N" is appended to done.txt.
HERE="$(cd "$(dirname "$0")" && pwd)"
TAG=$1; COND=$2; shift 2
PY=/c/Users/ds3420/AppData/Local/anaconda3/envs/PyBaMM_Pressure/python.exe
( cd "$HERE" && env RPT_OUT_DIR="$HERE" RPT_TAG="$TAG" "$@" \
    "$PY" ../rpt_soc_plots.py "$COND" --aligned > "log_$TAG.txt" 2>&1
  echo "$TAG exit=$?" >> "$HERE/done.txt" )
