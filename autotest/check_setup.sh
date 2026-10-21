#!/usr/bin/env bash
# Report, for each FreeCAD in versions.conf: version, Python, OCC, user Mod dir,
# which GDML workbench it will load (and git commit), and whether Fasteners is
# installed (some test files need it).
#
#   autotest/check_setup.sh            report only
#   autotest/check_setup.sh --link     also symlink <Mod>/GDML -> $GDML_REPO
#                                      when that Mod dir has no GDML at all
#
# GDML_REPO defaults to ~/Workbenches/GDML
set -u
. "$(dirname "$0")/fclib.sh"

LINK=0
[ "${1:-}" = "--link" ] && LINK=1
GDML_REPO="${GDML_REPO:-$HOME/Workbenches/GDML}"

PROBE_DIR="$(mktemp -d)"
PROBE="$PROBE_DIR/gdml_probe.py"
cat > "$PROBE" <<'EOF'
import os, sys, FreeCAD
print("AT_VERSION=" + ".".join(str(v) for v in FreeCAD.Version()[:3]) + " " + str(FreeCAD.Version()[3]))
print("AT_PYTHON=" + sys.version.split()[0])
try:
    import Part
    print("AT_OCC=" + str(Part.OCC_VERSION))
except Exception as e:
    print("AT_OCC=? " + str(e))
mod = os.path.join(FreeCAD.getUserAppDataDir(), "Mod")
print("AT_MODDIR=" + mod)
try:
    import freecad.gdml
    print("AT_GDML=" + os.path.realpath(os.path.dirname(freecad.gdml.__file__)))
except Exception as e:
    print("AT_GDML=NOT IMPORTABLE: " + str(e))
try:
    import lxml
    print("AT_LXML=yes " + lxml.__version__)
except Exception:
    print("AT_LXML=NO - GDML export falls back to xml.etree")
names = sorted(os.listdir(mod)) if os.path.isdir(mod) else []
print("AT_FASTENERS=" + ("yes" if any(n.lower().startswith("fasteners") for n in names) else "no"))
print("AT_MODS=" + ",".join(names))
EOF

at_versions | while IFS='|' read -r label path; do
    echo "=== $label  ($path)"
    if [ ! -e "$path" ]; then echo "    MISSING install"; continue; fi
    cmd="$(at_cmd_bin "$path")"
    if [ -z "$cmd" ]; then echo "    no freecadcmd found in bundle"; continue; fi
    out="$("$cmd" "$PROBE" </dev/null 2>&1)"
    get() { echo "$out" | sed -n "s/^$1=//p" | head -1; }
    echo "    FreeCAD : $(get AT_VERSION)"
    echo "    Python  : $(get AT_PYTHON)    OCC: $(get AT_OCC)"
    moddir="$(get AT_MODDIR)"
    echo "    Mod dir : $moddir"
    gdml="$(get AT_GDML)"
    echo "    GDML    : $gdml"
    case "$gdml" in
        /*) want="$(cd "$GDML_REPO/freecad/gdml" 2>/dev/null && pwd -P)"
            if [ -n "$want" ] && [ "$gdml" != "$want" ]; then
                echo "    !! NOT your checkout ($want) - probably a copy made by FreeCAD's"
                echo "       'migrate configuration' step. Replace $moddir/GDML with a symlink to $GDML_REPO"
            fi
            echo "    GDML git: $(git -C "$gdml" rev-parse --abbrev-ref HEAD 2>/dev/null)@$(git -C "$gdml" rev-parse --short HEAD 2>/dev/null)" ;;
    esac
    echo "    lxml    : $(get AT_LXML)"
    echo "    Fasteners: $(get AT_FASTENERS)"
    echo "    Mods    : $(get AT_MODS)"
    case "$gdml" in
        NOT*|"")
            if [ "$LINK" = 1 ] && [ -n "$moddir" ]; then
                if [ -e "$moddir/GDML" ] || [ -L "$moddir/GDML" ]; then
                    echo "    -> $moddir/GDML exists but did not import; not touching it"
                else
                    mkdir -p "$moddir" && ln -s "$GDML_REPO" "$moddir/GDML" \
                        && echo "    -> linked $moddir/GDML -> $GDML_REPO"
                fi
            else
                echo "    -> GDML not loadable here. Re-run with --link to symlink $GDML_REPO"
            fi ;;
    esac
    if [ -z "$(get AT_VERSION)" ]; then
        echo "    probe output:"; echo "$out" | tail -15 | sed 's/^/      /'
    fi
done
rm -rf "$PROBE_DIR"
