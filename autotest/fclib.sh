# Shared helpers for the autotest shell scripts (bash 3.2 compatible).

AT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AT_TC="$(cd "$AT_DIR/.." && pwd)"

# Print "label|path" for every active line of versions.conf
at_versions() {
    local conf="${AT_VERSIONS_CONF:-$AT_DIR/versions.conf}"
    sed -e 's/#.*//' "$conf" | while IFS='|' read -r label path; do
        label="$(echo "$label" | xargs)"
        path="$(echo "$path" | xargs)"
        [ -n "$label" ] && [ -n "$path" ] && echo "$label|$path"
    done
}

# GUI executable for an .app bundle or a direct executable
at_gui_bin() {
    local p="$1"
    if [ -d "$p" ] && [ -f "$p/Contents/Info.plist" ]; then
        local exe
        exe="$(/usr/bin/plutil -extract CFBundleExecutable raw "$p/Contents/Info.plist" 2>/dev/null)"
        [ -z "$exe" ] && exe="FreeCAD"
        echo "$p/Contents/MacOS/$exe"
    else
        echo "$p"
    fi
}

# Console (no GUI) executable, used only for the setup check
at_cmd_bin() {
    local p="$1" c
    if [ -d "$p" ]; then
        for c in "$p/Contents/Resources/bin/freecadcmd" "$p/Contents/Resources/bin/FreeCADCmd" \
                 "$p/Contents/MacOS/freecadcmd" "$p/Contents/MacOS/FreeCADCmd"; do
            [ -x "$c" ] && { echo "$c"; return; }
        done
    fi
    for c in "$(dirname "$p")/freecadcmd" "$(dirname "$p")/FreeCADCmd"; do
        [ -x "$c" ] && { echo "$c"; return; }
    done
    echo ""
}

at_python() {
    if command -v python3 >/dev/null 2>&1; then echo python3; else echo python; fi
}
