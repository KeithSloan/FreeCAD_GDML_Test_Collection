#!/usr/bin/env bash
# Run the GDML test collection against every FreeCAD in versions.conf, then
# write a comparison report. Safe to leave running unattended: a hung stage is
# killed, a crash is recorded against the file that caused it, and FreeCAD is
# relaunched to carry on with the next file.
#
#   autotest/run_versions.sh                 all versions in versions.conf
#   autotest/run_versions.sh fc-26.3rc2      just the named label(s)
#
# Environment (all optional):
#   AT_RESULTS          results root          (default <collection>/autotest_results)
#   AT_RUN_ID           run directory name    (default timestamp; reuse one to resume)
#   AT_STAGE_TIMEOUT    secs one stage may take before FreeCAD is killed (default 900)
#   AT_STARTUP_TIMEOUT  secs for FreeCAD to start the runner              (default 300)
#   AT_MAX_RESTARTS     relaunches per version after crashes/hangs        (default 25)
#   AT_STAGES           open,export,cog,import                            (default all)
#   AT_ONLY             glob(s) restricting files, comma-separated, e.g. 'Loft*,Sweep*'
#   AT_GEANT4=1         also run autotest/run_geant4.sh on the exports
#   AT_HEADLESS_X=1     Linux only: run FreeCAD under xvfb-run
#   AT_SHARED_SESSION=1 test all files in one FreeCAD (faster; by default each
#                       file gets a fresh FreeCAD so results can't depend on
#                       which files were opened before it)
set -u
. "$(dirname "$0")/fclib.sh"

STAGE_TIMEOUT="${AT_STAGE_TIMEOUT:-900}"
STARTUP_TIMEOUT="${AT_STARTUP_TIMEOUT:-300}"
MAX_RESTARTS="${AT_MAX_RESTARTS:-25}"
RUN_ID="${AT_RUN_ID:-$(date +%Y%m%d-%H%M%S)}"
RUN_DIR="${AT_RESULTS:-$AT_TC/autotest_results}/$RUN_ID"
PY="$(at_python)"
mkdir -p "$RUN_DIR"

if [ "$(uname)" = Darwin ]; then
    mtime() { stat -f %m "$1"; }
else
    mtime() { stat -c %Y "$1"; }
fi
say() { echo "[$(date +%H:%M:%S)] $*" | tee -a "$RUN_DIR/run.log"; }

if pgrep -x FreeCAD >/dev/null 2>&1 || pgrep -x freecad >/dev/null 2>&1; then
    say "WARNING: a FreeCAD is already running. If its 'single instance' option is on,"
    say "         the test script may be handed to it instead. Quit it if results look empty."
fi

VLIST="$RUN_DIR/versions.txt"
at_versions > "$VLIST"

run_one() {
    local label="$1" path="$2"
    local out="$RUN_DIR/$label" gui log pid rc reason launched now hb age recorded restarts=0
    mkdir -p "$out"
    log="$out/freecad.log"
    gui="$(at_gui_bin "$path")"
    if [ ! -x "$gui" ]; then say "$label: no executable at $gui - skipped"; return; fi
    if [ -f "$out/DONE" ]; then say "$label: already complete - skipped"; return; fi

    # A heartbeat left by an earlier, interrupted run marks the file that
    # process died on: record it before launching, or it looks like a hang.
    if [ -f "$out/current.json" ]; then
        "$PY" "$AT_DIR/record_abort.py" "$out" "CRASH" "previous-run" >>"$RUN_DIR/run.log" 2>&1
        tail -1 "$RUN_DIR/run.log"
    fi

    while :; do
        say "$label: launching FreeCAD (restart $restarts)"
        rm -f "$out/runner_error.txt" "$out/env.json"
        if [ "${AT_HEADLESS_X:-0}" = 1 ]; then
            AT_OUT="$out" AT_TC="$AT_TC" AT_CONSOLE_LOG="$out/console.log" \
                xvfb-run -a -s "-screen 0 1600x1200x24" \
                "$gui" --log-file "$out/console.log" "$AT_DIR/gdml_autotest_runner.py" </dev/null >>"$log" 2>&1 &
        else
            AT_OUT="$out" AT_TC="$AT_TC" AT_CONSOLE_LOG="$out/console.log" \
                "$gui" --log-file "$out/console.log" "$AT_DIR/gdml_autotest_runner.py" </dev/null >>"$log" 2>&1 &
        fi
        pid=$!
        launched=$(date +%s)
        reason=""
        while kill -0 "$pid" 2>/dev/null; do
            sleep 5
            now=$(date +%s)
            hb="$out/current.json"
            if [ -f "$hb" ]; then
                age=$(( now - $(mtime "$hb") ))
                if [ "$age" -gt "$STAGE_TIMEOUT" ]; then
                    reason="TIMEOUT"
                    say "$label: stage stuck for ${age}s - killing FreeCAD"
                fi
            elif [ ! -f "$out/env.json" ] && [ $(( now - launched )) -gt "$STARTUP_TIMEOUT" ]; then
                reason="STARTUP_TIMEOUT"
                say "$label: runner never started - killing FreeCAD"
            fi
            if [ -n "$reason" ]; then
                pkill -9 -P "$pid" 2>/dev/null
                kill -9 "$pid" 2>/dev/null
                break
            fi
        done
        wait "$pid" 2>/dev/null
        rc=$?
        # Each launch rewrites console.log and report_view.txt: keep them all.
        [ -f "$out/console.log" ] && cat "$out/console.log" >> "$out/console_all.log"
        [ -f "$out/report_view.txt" ] && cat "$out/report_view.txt" >> "$out/report_view_all.txt"
        rm -f "$out/report_view.txt"

        if [ -f "$out/DONE" ]; then say "$label: finished"; break; fi
        if [ -f "$out/NEXT" ] && [ -z "$reason" ]; then
            # Runner tested one file and asked for a fresh FreeCAD
            # (default: one file per process). Not a crash.
            rm -f "$out/NEXT"
            continue
        fi
        if [ -f "$out/runner_error.txt" ]; then
            say "$label: runner aborted - see $out/runner_error.txt"; break
        fi
        [ -z "$reason" ] && reason="CRASH"
        recorded=1
        if [ "$reason" != "STARTUP_TIMEOUT" ]; then
            "$PY" "$AT_DIR/record_abort.py" "$out" "$reason" "$rc" >>"$RUN_DIR/run.log" 2>&1 && recorded=0
            tail -1 "$RUN_DIR/run.log"
        fi
        if [ "$recorded" != 0 ]; then
            say "$label: FreeCAD exited (rc=$rc) before reaching a test file."
            say "        Open this FreeCAD once by hand to clear any first-run dialogs, then rerun with AT_RUN_ID=$RUN_ID"
            break
        fi
        restarts=$((restarts + 1))
        if [ "$restarts" -gt "$MAX_RESTARTS" ]; then
            say "$label: too many restarts - giving up"; break
        fi
    done
}

exec 3< "$VLIST"
while IFS='|' read -r label path <&3; do
    if [ $# -gt 0 ]; then
        case " $* " in *" $label "*) ;; *) continue ;; esac
    fi
    run_one "$label" "$path"
done
exec 3<&-

if [ "${AT_GEANT4:-0}" = 1 ]; then
    "$AT_DIR/run_geant4.sh" "$RUN_DIR" 2>&1 | tee -a "$RUN_DIR/run.log"
fi

"$PY" "$AT_DIR/compare.py" "$RUN_DIR" | tee -a "$RUN_DIR/run.log"
say "report: $RUN_DIR/report.md"
