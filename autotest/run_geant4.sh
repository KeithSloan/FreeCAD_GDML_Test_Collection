#!/usr/bin/env bash
# Optional Geant4 stage: load every exported GDML of every version in a run
# with load_gdml_color (built from MH_G4_Example/MH_G01) and record whether
# Geant4 accepted it.
#
#   autotest/run_geant4.sh <run_dir>
#
# Writes <run_dir>/<version>/geant4.jsonl and geant4_out/<file>.out
# Status per file: OK, FAILED (non-zero exit), ABORT ("Aborting" in output,
# e.g. G4Exception), TIMEOUT.  AT_G4_TIMEOUT secs per file (default 600).
# Needs the Geant4 environment sourced (geant4.sh) so load_gdml_color is on PATH.
set -u
. "$(dirname "$0")/fclib.sh"
RUN_DIR="${1:?usage: run_geant4.sh <run_dir>}"
G4_TIMEOUT="${AT_G4_TIMEOUT:-600}"

if ! command -v load_gdml_color >/dev/null 2>&1; then
    echo "run_geant4: load_gdml_color not on PATH - Geant4 stage skipped"
    exit 0
fi

for vdir in "$RUN_DIR"/*/; do
    [ -d "$vdir/exported_gdml" ] || continue
    out="$vdir/geant4_out"
    mkdir -p "$out"
    : > "$vdir/geant4.jsonl"
    cp "$AT_TC/batch.mac" "$vdir/exported_gdml/batch.mac"
    echo "run_geant4: $(basename "$vdir")"
    for g in "$vdir"/exported_gdml/*.gdml; do
        [ -f "$g" ] || continue
        f="$(basename "$g")"
        ( cd "$vdir/exported_gdml" && \
          perl -e 'alarm shift; exec @ARGV or die' "$G4_TIMEOUT" \
              load_gdml_color "$f" -printCM -m=batch.mac ) >"$out/$f.out" 2>&1
        rc=$?
        if [ "$rc" = 142 ]; then st=TIMEOUT
        elif grep -q "Aborting" "$out/$f.out"; then st=ABORT
        elif [ "$rc" != 0 ]; then st=FAILED
        else st=OK; fi
        printf '{"gdml": "%s", "status": "%s", "exit_code": %s}\n' "$f" "$st" "$rc" >> "$vdir/geant4.jsonl"
        echo "  $st  $f"
    done
done
