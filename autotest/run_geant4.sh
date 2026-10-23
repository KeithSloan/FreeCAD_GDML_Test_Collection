#!/usr/bin/env bash
# Optional Geant4 stage: load every exported GDML of every version in a run
# with Munther's load_gdml_color (MH_G4_Example/MH_G01) using -printCM, record
# whether Geant4 accepted it, then compare Geant4's volume / centre of mass /
# inertia with FreeCAD's (autotest/geant4_cog.py).
#
#   autotest/run_geant4.sh <run_dir>
#
# Writes <run_dir>/<version>/geant4.jsonl, geant4_out/<file>.out,
# geant4_cog.jsonl and <run_dir>/geant4_cog.md.
# Status per file: OK, FAILED (non-zero exit), ABORT ("Aborting" in output,
# e.g. G4Exception), TIMEOUT.
#
# Environment:
#   AT_G4_BIN      load_gdml_color to use (default: on PATH, else
#                  MH_G4_Example/MH_G01/build-conda/, else .../build/)
#   AT_G4_CONDA    conda env prefix with Geant4 (default ~/miniconda3/envs/g4);
#                  activated to set Geant4's data paths
#   AT_G4_ENV      geant4.sh to source instead (non-conda Geant4 installs)
#   (neither is used if Geant4's environment is already set)
#   AT_G4_MAC      batch macro (default MH_G4_Example/Geant4_COG/batch.mac)
#   AT_G4_TIMEOUT  secs per file (default 600)
#   AT_G4_JOBS     files run at once (default 4)
#   AT_G4_LOGNSIM  log10 of Monte Carlo points per volume (default: Geant4's 6)
set -u
. "$(dirname "$0")/fclib.sh"
RUN_DIR="${1:?usage: run_geant4.sh <run_dir>}"
RUN_DIR="$(cd "$RUN_DIR" && pwd)"
G4_TIMEOUT="${AT_G4_TIMEOUT:-600}"
G4_JOBS="${AT_G4_JOBS:-4}"
G4_MAC="${AT_G4_MAC:-$AT_TC/MH_G4_Example/Geant4_COG/batch.mac}"

G4_BIN="${AT_G4_BIN:-}"
if [ -z "$G4_BIN" ]; then
    if command -v load_gdml_color >/dev/null 2>&1; then
        G4_BIN="$(command -v load_gdml_color)"
    else
        for b in build-conda build; do
            G4_BIN="$AT_TC/MH_G4_Example/MH_G01/$b/load_gdml_color"
            [ -x "$G4_BIN" ] && break
        done
    fi
fi
if [ ! -x "$G4_BIN" ]; then
    echo "run_geant4: no load_gdml_color ($G4_BIN) - Geant4 stage skipped"
    echo "  build it: see MH_G4_Example/README, or set AT_G4_BIN"
    exit 0
fi
# Geant4 needs its data-set variables (G4ENSDFSTATEDATA etc.) at run time
if [ -z "${G4ENSDFSTATEDATA:-}" ]; then
    G4_CONDA="${AT_G4_CONDA:-$HOME/miniconda3/envs/g4}"
    set +u   # conda's activate.d scripts and geant4.sh read unset variables
    if [ -n "${AT_G4_ENV:-}" ]; then
        . "$AT_G4_ENV" >/dev/null 2>&1
    elif [ -d "$G4_CONDA/conda-meta" ]; then
        # conda's Geant4 sets its data paths on activation (its geant4.sh refuses)
        . "$(cd "$G4_CONDA/../.." && pwd)/bin/activate" "$G4_CONDA" >/dev/null 2>&1
    fi
    set -u
    if [ -z "${G4ENSDFSTATEDATA:-}" ]; then
        echo "run_geant4: Geant4 environment not set - set AT_G4_CONDA (conda env prefix) or AT_G4_ENV (geant4.sh)"
        exit 1
    fi
fi
[ -f "$G4_MAC" ] || { echo "run_geant4: batch macro $G4_MAC not found"; exit 1; }

echo "run_geant4: $G4_BIN"
echo "run_geant4: macro $G4_MAC, $G4_JOBS at once, timeout ${G4_TIMEOUT}s"
NSIM_ARG=""
[ -n "${AT_G4_LOGNSIM:-}" ] && NSIM_ARG="-logNsim=$AT_G4_LOGNSIM"
export G4_BIN G4_TIMEOUT NSIM_ARG

# one file: run in exported_gdml (GDML may refer to files beside it)
run_one() {
    g="$1"; out="$2"
    f="$(basename "$g")"
    t0=$(date +%s)
    ( cd "$(dirname "$g")" && \
      perl -e 'alarm shift; exec @ARGV or die' "$G4_TIMEOUT" \
          "$G4_BIN" "$f" -printCM $NSIM_ARG -m=batch.mac ) >"$out/$f.out" 2>&1
    rc=$?
    if [ "$rc" = 142 ]; then st=TIMEOUT
    elif grep -q "Aborting" "$out/$f.out"; then st=ABORT
    elif [ "$rc" != 0 ]; then st=FAILED
    else st=OK; fi
    printf '{"gdml": "%s", "status": "%s", "exit_code": %s, "secs": %s}\n' \
        "$f" "$st" "$rc" "$(( $(date +%s) - t0 ))" > "$out/$f.status"
    echo "  $st  $f"
}
export -f run_one

for vdir in "$RUN_DIR"/*/; do
    [ -d "$vdir/exported_gdml" ] || continue
    out="${vdir%/}/geant4_out"
    mkdir -p "$out"
    rm -f "$out"/*.status
    cp "$G4_MAC" "$vdir/exported_gdml/batch.mac"
    echo "run_geant4: $(basename "$vdir")"
    ls "$vdir"/exported_gdml/*.gdml 2>/dev/null | \
        xargs -P "$G4_JOBS" -I{} bash -c 'run_one "$1" "$2"' _ {} "$out"
    cat "$out"/*.status > "$vdir/geant4.jsonl" 2>/dev/null || : > "$vdir/geant4.jsonl"
done

python3 "$AT_DIR/geant4_cog.py" "$RUN_DIR" > /dev/null && \
    echo "run_geant4: COG comparison in $RUN_DIR/geant4_cog.md"
