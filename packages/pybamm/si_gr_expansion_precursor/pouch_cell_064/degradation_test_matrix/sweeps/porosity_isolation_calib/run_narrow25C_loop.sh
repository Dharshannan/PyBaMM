#!/bin/bash
export C64_SOLVENT_MULT=5
export C64_SI_CRIT_STRESS=6e6
export C64_SI_LAM_EAC=40000
export C64_WIDTH=0.0003
export C64_F0=0.85
export C64_SI_LAM_ISO_EXPONENT=15.0
export C64_SI_ESEI=38000
export C64_GR_ESEI=38000
export C64_AMBIENT_T=298.15
export C64_EFC_LIMIT=420
export C64_MAX_CYCLES=1500
export C64_CHUNK_CYCLES=300
export C64_RESUME_PKL=narrow25C_chunked_state.pkl
PY="/c/Users/ds3420/AppData/Local/anaconda3/envs/PyBaMM_Pressure/python.exe"
for i in $(seq 1 8); do
  echo "=== CHUNK $i ===" >> log_narrow25C_chunked.txt
  "$PY" run_chunked_narrow25C.py >> log_narrow25C_chunked.txt 2>&1
  if grep -q "OVERALL STOP" log_narrow25C_chunked.txt; then
    echo "=== LOOP DONE (overall stop reached) ===" >> log_narrow25C_chunked.txt
    break
  fi
done
