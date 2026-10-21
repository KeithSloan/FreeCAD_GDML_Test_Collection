# ***************************************************************************
# *   GDML workbench automated test runner                                  *
# *                                                                         *
# *   Runs INSIDE the FreeCAD GUI (launched by run_versions.sh):            *
# *       FreeCAD gdml_autotest_runner.py                                   *
# *   The GUI is required because exportGDML takes the world volume from    *
# *   FreeCADGui.Selection and calcCenterOfMass walks the tree widget.      *
# *                                                                         *
# *   Configuration comes from environment variables:                       *
# *     AT_OUT       output directory for this FreeCAD version (required)   *
# *     AT_TC        root of FreeCAD_GDML_Test_Collection (required)        *
# *     AT_FC_SRC    FC files dir     (default $AT_TC/FC_test_files)        *
# *     AT_GDML_SRC  GDML files dir   (default $AT_TC/test_gdml_files)      *
# *     AT_STAGES    comma list: open,export,cog,import (default all)       *
# *     AT_ONLY      optional glob to restrict files, e.g. "Loft*"          *
# *                                                                         *
# *   Results: $AT_OUT/results.jsonl (one JSON line per file), env.json,    *
# *   exported_gdml/, cogs/, roundtrip_gdml/, dialogs.log, DONE marker.     *
# *   Already-recorded files are skipped, so a relaunch after a crash       *
# *   resumes where it stopped.                                             *
# ***************************************************************************

import fnmatch
import glob
import importlib.util
import json
import os
import platform
import re
import subprocess
import sys
import time
import traceback

import FreeCAD
import FreeCADGui

from PySide import QtCore, QtWidgets

OUT = os.environ.get("AT_OUT", "")
TC = os.environ.get("AT_TC", "")
FC_SRC = os.environ.get("AT_FC_SRC") or os.path.join(TC, "FC_test_files")
GDML_SRC = os.environ.get("AT_GDML_SRC") or os.path.join(TC, "test_gdml_files")
STAGES = set(s.strip() for s in os.environ.get("AT_STAGES", "open,export,cog,import").split(",") if s.strip())
ONLY = os.environ.get("AT_ONLY", "")
CONSOLE_LOG = os.environ.get("AT_CONSOLE_LOG", "")

RESULTS = os.path.join(OUT, "results.jsonl")
CURRENT = os.path.join(OUT, "current.json")
DIALOGS = os.path.join(OUT, "dialogs.log")

_state = {"file": None, "stage": None, "dialogs": []}
_timer = None


# -------------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------------

def log(msg):
    line = f"[autotest {time.strftime('%H:%M:%S')}] {msg}"
    print(line)
    try:
        FreeCAD.Console.PrintMessage(line + "\n")
    except Exception:
        pass


def heartbeat(file=None, stage=None):
    """Written at every stage so the shell watchdog can spot hangs and crashes."""
    if file is not None:
        _state["file"] = file
    _state["stage"] = stage
    with open(CURRENT, "w", encoding="utf-8") as f:
        json.dump({"file": _state["file"], "stage": stage, "t": time.time()}, f)


def appendResult(rec):
    with open(RESULTS, "a", encoding="utf-8") as f:
        f.write(json.dumps(rec) + "\n")


def doneKeys():
    keys = set()
    if os.path.isfile(RESULTS):
        with open(RESULTS, encoding="utf-8") as f:
            for line in f:
                try:
                    keys.add(json.loads(line)["key"])
                except Exception:
                    pass
    return keys


def processEvents():
    try:
        QtWidgets.QApplication.processEvents()
    except Exception:
        pass


def dismissModals():
    """Close any modal dialog so an unattended run never blocks.
    Every dismissed dialog is logged against the file being processed."""
    w = QtWidgets.QApplication.activeModalWidget()
    if w is None:
        return
    title = w.windowTitle()
    text = ""
    if isinstance(w, QtWidgets.QMessageBox):
        text = w.text()
    entry = {"file": _state["file"], "stage": _state["stage"], "title": title, "text": text[:500]}
    _state["dialogs"].append(entry)
    with open(DIALOGS, "a", encoding="utf-8") as f:
        f.write(json.dumps(entry) + "\n")
    try:
        w.reject()
    except Exception:
        w.close()


def gitCommit(path):
    try:
        d = path if os.path.isdir(path) else os.path.dirname(path)
        r = subprocess.run(["git", "-C", d, "rev-parse", "--short", "HEAD"],
                           capture_output=True, text=True, timeout=10)
        b = subprocess.run(["git", "-C", d, "rev-parse", "--abbrev-ref", "HEAD"],
                           capture_output=True, text=True, timeout=10)
        return (b.stdout.strip() + "@" + r.stdout.strip()) if r.returncode == 0 else None
    except Exception:
        return None


def writeEnv():
    import freecad.gdml as gdmlpkg
    try:
        import Part
        occ = Part.OCC_VERSION
    except Exception:
        occ = None
    gdmlDir = os.path.dirname(os.path.realpath(gdmlpkg.__file__))
    try:
        import freecad.gdml.exportGDML as exportGDML
        xmlBackend = getattr(exportGDML, "XML_IO_VERSION", "?")
    except Exception as e:
        xmlBackend = f"exportGDML import failed: {e}"
    env = {
        "freecad_version": ".".join(str(v) for v in FreeCAD.Version()[:3]),
        "freecad_version_full": [str(v) for v in FreeCAD.Version()],
        "python": sys.version.split()[0],
        "occ": occ,
        "platform": platform.platform(),
        "machine": platform.machine(),
        "gdml_path": gdmlDir,
        "gdml_commit": gitCommit(gdmlDir),
        "gdml_xml_backend": xmlBackend,
        "user_app_data": FreeCAD.getUserAppDataDir(),
        "test_collection_commit": gitCommit(TC),
        "stages": sorted(STAGES),
        "started": time.strftime("%Y-%m-%d %H:%M:%S"),
    }
    with open(os.path.join(OUT, "env.json"), "w", encoding="utf-8") as f:
        json.dump(env, f, indent=2)
    log(f"FreeCAD {env['freecad_version']} Python {env['python']} OCC {occ}")
    log(f"GDML from {gdmlDir} ({env['gdml_commit']}) XML backend {xmlBackend}")
    return env


def loadCOGModule():
    path = os.path.join(TC, "calcCenterOfMass.py")
    spec = importlib.util.spec_from_file_location("calcCenterOfMass", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["calcCenterOfMass"] = mod
    spec.loader.exec_module(mod)
    return mod


def invalidObjects(doc):
    bad = []
    for o in doc.Objects:
        try:
            if "Invalid" in o.State or not o.isValid():
                bad.append(o.Label)
        except Exception:
            pass
    return bad


def rootPart(doc):
    for obj in doc.RootObjects:
        if obj.TypeId == "App::Part":
            return obj
    return None


def select(obj):
    FreeCADGui.Selection.clearSelection()
    FreeCADGui.Selection.addSelection(obj)


def closeAll():
    for name in list(FreeCAD.listDocuments().keys()):
        try:
            FreeCAD.closeDocument(name)
        except Exception:
            pass
    processEvents()


def nonEmpty(path):
    return bool(path) and os.path.isfile(path) and os.path.getsize(path) > 0


# -------------------------------------------------------------------------
# Console capture: a file that needs an addon this FreeCAD can't import
# (Fasteners, Gears, ...) still opens, using the shapes stored in it, so it
# would look like a pass. FreeCAD reports the failed import on the console,
# so scan what was printed while the file was opened.
# GDML's own groups have no proxy, so a proxy check can't be used instead.
# -------------------------------------------------------------------------

_MISSING_RE = re.compile(r"no module named|blocked import", re.IGNORECASE)
_MODNAME_RE = re.compile(r"""(?:no module named|blocked import(?: of)?(?: module)?)\W*([A-Za-z_][\w.]*)""",
                         re.IGNORECASE)


def _reportView():
    try:
        mw = FreeCADGui.getMainWindow()
        dock = mw.findChild(QtWidgets.QDockWidget, "Report view")
        w = dock.widget() if dock else None
        return w if w is not None and hasattr(w, "toPlainText") else None
    except Exception:
        return None


def consoleMark():
    logSize = os.path.getsize(CONSOLE_LOG) if CONSOLE_LOG and os.path.isfile(CONSOLE_LOG) else 0
    rv = _reportView()
    return (logSize, len(rv.toPlainText()) if rv is not None else 0)


def consoleSince(mark):
    text = ""
    if CONSOLE_LOG and os.path.isfile(CONSOLE_LOG):
        try:
            with open(CONSOLE_LOG, "rb") as f:
                f.seek(mark[0])
                text += f.read().decode("utf-8", "replace")
        except Exception:
            pass
    rv = _reportView()
    if rv is not None:
        try:
            text += "\n" + rv.toPlainText()[mark[1]:]
        except Exception:
            pass
    return text


def missingAddons(text):
    """(module names, matching console lines) for failed imports in text."""
    lines, mods = [], []
    for line in text.splitlines():
        if _MISSING_RE.search(line):
            line = line.strip()
            if line not in lines:
                lines.append(line[:300])
            m = _MODNAME_RE.search(line)
            if m:
                name = m.group(1).strip(".")
                if name and name not in mods:
                    mods.append(name)
    return mods, lines


def stopOwnMCPServer():
    """Keep test instances off the FreeCAD MCP port: if the FreeCADMCP addon
    auto-started its server in this process, stop it."""
    try:
        mod = sys.modules.get("rpc_server.rpc_server")
        if mod is not None and getattr(mod, "rpc_server_instance", None) is not None:
            mod.stop_rpc_server()
            log("stopped FreeCAD MCP server in this test instance")
    except Exception as e:
        log(f"could not stop MCP server: {e}")


def runStage(rec, stage, fn):
    """Run one stage, record ok/error/time; never raises."""
    heartbeat(stage=stage)
    t0 = time.time()
    try:
        rec[stage] = fn() or {}
        rec[stage].setdefault("ok", True)
    except BaseException:
        rec[stage] = {"ok": False, "error": traceback.format_exc()[-4000:]}
    rec[stage]["secs"] = round(time.time() - t0, 3)
    processEvents()
    return rec[stage]["ok"]


# -------------------------------------------------------------------------
# FC file tests: open/recompute -> export GDML -> COG
# -------------------------------------------------------------------------

def testFCFile(fcPath, cogMod):
    import freecad.gdml.exportGDML as exportGDML

    name = os.path.basename(fcPath)
    stem = name[:-6]
    rec = {"key": "fc:" + name, "kind": "fc", "file": name}
    _state["dialogs"] = []
    heartbeat(file=name, stage="start")
    ctx = {}

    def stOpen():
        mark = consoleMark()
        doc = FreeCAD.openDocument(fcPath)
        processEvents()
        ctx["doc"] = doc
        before = invalidObjects(doc)
        nRecomputed = doc.recompute()
        processEvents()
        after = invalidObjects(doc)
        part = rootPart(doc)
        ctx["part"] = part
        mods, lines = missingAddons(consoleSince(mark))
        rec["missing_addons"] = mods or (["?"] if lines else [])
        rec["missing_addon_messages"] = lines
        return {
            "objects": len(doc.Objects),
            "invalid_before_recompute": before,
            "invalid_after_recompute": after,
            "recomputed": nRecomputed,
            "root_part": part.Label if part else None,
            "ok": part is not None,
            "error": None if part else "no root App::Part found",
        }

    def stExport():
        part = ctx["part"]
        path = os.path.join(OUT, "exported_gdml", f"{stem}-{part.Label}.gdml")
        if os.path.exists(path):
            os.remove(path)
        select(part)
        exportGDML.export([part], path)
        ok = nonEmpty(path)
        return {"path": os.path.relpath(path, OUT), "bytes": os.path.getsize(path) if ok else 0,
                "ok": ok, "error": None if ok else "export produced no file"}

    def stCOG():
        part = ctx["part"]
        path = os.path.join(OUT, "cogs", f"{stem}.txt")
        select(part)
        cogMod.calcCenterOfMass(path)   # file is closed when the function returns
        ok = nonEmpty(path)
        text = open(path, encoding="utf-8").read() if ok else ""
        return {"path": os.path.relpath(path, OUT), "text": text,
                "ok": ok, "error": None if ok else "no COG output"}

    try:
        okOpen = runStage(rec, "open", stOpen)
        if okOpen and "export" in STAGES:
            runStage(rec, "export", stExport)
        if okOpen and "cog" in STAGES:
            runStage(rec, "cog", stCOG)
    finally:
        heartbeat(stage="close")
        closeAll()

    rec["dialogs"] = list(_state["dialogs"])
    rec["status"] = summarise(rec, ["open", "export", "cog"])
    appendResult(rec)
    log(f"{name}: {rec['status']}")


# -------------------------------------------------------------------------
# GDML file tests: import -> recompute -> re-export (round trip)
# -------------------------------------------------------------------------

def testGDMLFile(gdmlPath):
    import freecad.gdml.exportGDML as exportGDML
    from freecad.gdml.importGDML import processGDML

    name = os.path.basename(gdmlPath)
    stem = os.path.splitext(name)[0]
    rec = {"key": "gdml:" + name, "kind": "gdml", "file": name}
    _state["dialogs"] = []
    heartbeat(file=name, stage="start")
    ctx = {}

    def stImport():
        doc = FreeCAD.newDocument(stem.replace("-", "_").replace(".", "_"))
        processGDML(doc, False, gdmlPath, False, 1, 0)
        processEvents()
        doc.recompute()
        processEvents()
        part = rootPart(doc)
        ctx["part"] = part
        return {"objects": len(doc.Objects),
                "invalid_after_recompute": invalidObjects(doc),
                "root_part": part.Label if part else None,
                "ok": part is not None,
                "error": None if part else "no root App::Part after import"}

    def stRoundtrip():
        part = ctx["part"]
        path = os.path.join(OUT, "roundtrip_gdml", f"{stem}-{part.Label}.gdml")
        if os.path.exists(path):
            os.remove(path)
        select(part)
        exportGDML.export([part], path)
        ok = nonEmpty(path)
        return {"path": os.path.relpath(path, OUT), "bytes": os.path.getsize(path) if ok else 0,
                "ok": ok, "error": None if ok else "re-export produced no file"}

    try:
        if runStage(rec, "import", stImport):
            runStage(rec, "roundtrip", stRoundtrip)
    finally:
        heartbeat(stage="close")
        closeAll()

    rec["dialogs"] = list(_state["dialogs"])
    rec["status"] = summarise(rec, ["import", "roundtrip"])
    appendResult(rec)
    log(f"{name}: {rec['status']}")


def summarise(rec, stages):
    # Not a real test: FreeCAD fell back to the shapes stored in the file.
    if rec.get("missing_addons"):
        return "MISSING_ADDON"
    for st in stages:
        if st in rec and not rec[st].get("ok"):
            return f"FAIL_{st.upper()}"
    if rec.get("dialogs"):
        return "OK_WITH_DIALOGS"
    return "OK"


# -------------------------------------------------------------------------
# Main
# -------------------------------------------------------------------------

def selected(paths):
    paths = sorted(paths)
    if ONLY:
        paths = [p for p in paths if fnmatch.fnmatch(os.path.basename(p), ONLY)]
    return paths


def main():
    global _timer
    exitCode = 0
    try:
        if not OUT or not TC:
            raise RuntimeError("AT_OUT and AT_TC must be set")
        for d in ("exported_gdml", "cogs", "roundtrip_gdml"):
            os.makedirs(os.path.join(OUT, d), exist_ok=True)

        # No auto-recovery files: after a crash the relaunch must not stop
        # on FreeCAD's Document Recovery dialog.
        FreeCAD.ParamGet("User parameter:BaseApp/Preferences/Document").SetBool("SaveAutoRecovery", False)

        # Dialog watchdog: fires inside any modal exec() loop.
        _timer = QtCore.QTimer()
        _timer.timeout.connect(dismissModals)
        _timer.start(2000)

        # Workbench must be importable for document proxies to restore.
        import freecad.gdml  # noqa: F401
        stopOwnMCPServer()
        writeEnv()
        done = doneKeys()
        cogMod = loadCOGModule() if "cog" in STAGES else None

        if "open" in STAGES:
            for p in selected(glob.glob(os.path.join(FC_SRC, "*.FCStd"))):
                if "fc:" + os.path.basename(p) in done:
                    continue
                testFCFile(p, cogMod)

        if "import" in STAGES:
            for p in selected(glob.glob(os.path.join(GDML_SRC, "*.gdml"))):
                if "gdml:" + os.path.basename(p) in done:
                    continue
                testGDMLFile(p)

        with open(os.path.join(OUT, "DONE"), "w") as f:
            f.write(time.strftime("%Y-%m-%d %H:%M:%S\n"))
        if os.path.exists(CURRENT):
            os.remove(CURRENT)
        log("all tests finished")
    except BaseException:
        exitCode = 3
        msg = traceback.format_exc()
        log("runner aborted:\n" + msg)
        try:
            with open(os.path.join(OUT, "runner_error.txt"), "w") as f:
                f.write(msg)
        except Exception:
            pass
    finally:
        sys.stdout.flush()
        sys.stderr.flush()
        os._exit(exitCode)   # skip "save changes?" prompts on quit


# Defer until the GUI event loop is running so workbenches, the tree view
# and the main window are fully initialised.
QtCore.QTimer.singleShot(3000, main)
