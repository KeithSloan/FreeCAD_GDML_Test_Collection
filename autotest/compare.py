"""Compare GDML autotest results across FreeCAD versions and against the
committed baselines (FC_COGs/*.txt, FC_exported_gdml/*.gdml).

    python3 autotest/compare.py <run_dir>

Writes <run_dir>/report.md (for people) and <run_dir>/report.json (for agents).
Pure standard library; runs with any Python 3.8+.
"""
import glob
import json
import os
import re
import sys
import xml.etree.ElementTree as ET

HERE = os.path.dirname(os.path.abspath(__file__))
TC = os.path.dirname(HERE)

VOL_REL_TOL = 1e-3      # volume, relative
CM_REL_TOL = 1e-3       # centre of mass, relative to |cm| ...
CM_ABS_TOL = 0.02       # ... or this many mm (output is printed to 2 dp)
INERTIA_REL_TOL = 1e-2  # printed to 4 significant figures, may be Monte Carlo

COG_RE = re.compile(
    r'^(?P<name>.+?)\s+Volume:\s*(?P<vol>[-\d.eE+]+)\s*(?P<vunit>mm\^3|cm\^3|m\^3)\s*-->\s*'
    r'center of "mass":\s*\(\s*(?P<x>[-\d.eE+]+),\s*(?P<y>[-\d.eE+]+),\s*(?P<z>[-\d.eE+]+)\)\s*(?P<cunit>mm|cm|m),\s*'
    r'moments of inetria:\s*\(\s*(?P<i1>[-\d.eE+]+),\s*(?P<i2>[-\d.eE+]+),\s*(?P<i3>[-\d.eE+]+)\)\s*(?P<iunit>.*)$')

VOL_F = {"mm^3": 1.0, "cm^3": 1e3, "m^3": 1e9}
LEN_F = {"mm": 1.0, "cm": 10.0, "m": 1000.0}


def inertiaFactor(unit):
    unit = unit.strip()
    k = 0
    m = re.match(r"x10\^(\d+)\s+(.*)", unit)
    if m:
        k, unit = int(m.group(1)), m.group(2)
    base = 1e5 if unit.startswith("cm") else 1.0
    return base * 10 ** k


def parseCOG(text):
    rows = {}
    for line in text.splitlines():
        m = COG_RE.match(line.strip())
        if not m:
            continue
        g = m.groupdict()
        lf, inf = LEN_F[g["cunit"]], inertiaFactor(g["iunit"])
        rows[g["name"].strip()] = {
            "vol": float(g["vol"]) * VOL_F[g["vunit"]],
            "cm": [float(g[c]) * lf for c in "xyz"],
            "I": [float(g[c]) * inf for c in ("i1", "i2", "i3")],
        }
    return rows


def close(a, b, rel, abs_=0.0):
    return abs(a - b) <= max(abs_, rel * max(abs(a), abs(b)))


def compareCOG(a, b):
    """Return list of human-readable differences (empty = match)."""
    diffs = []
    for name in sorted(set(a) | set(b)):
        if name not in a or name not in b:
            diffs.append(f"row '{name}' only in {'first' if name in a else 'second'}")
            continue
        ra, rb = a[name], b[name]
        if not close(ra["vol"], rb["vol"], VOL_REL_TOL):
            diffs.append(f"{name}: volume {ra['vol']:.6g} vs {rb['vol']:.6g} mm^3")
        for i, ax in enumerate("xyz"):
            if not close(ra["cm"][i], rb["cm"][i], CM_REL_TOL, CM_ABS_TOL):
                diffs.append(f"{name}: cm.{ax} {ra['cm'][i]:.4g} vs {rb['cm'][i]:.4g} mm")
        for i in range(3):
            if not close(ra["I"][i], rb["I"][i], INERTIA_REL_TOL):
                diffs.append(f"{name}: I{i + 1} {ra['I'][i]:.4g} vs {rb['I'][i]:.4g} mm^5")
    return diffs


def gdmlSummary(path):
    """Structural fingerprint of a GDML file: what Geant4 would see."""
    try:
        root = ET.parse(path).getroot()
    except Exception as e:
        return {"parse_error": str(e)[:200]}
    s = {}
    sec = {c.tag: c for c in root}
    if "solids" in sec:
        for sol in sec["solids"]:
            s["solid:" + sol.tag] = s.get("solid:" + sol.tag, 0) + 1
            if sol.tag in ("tessellated",):
                s["facets"] = s.get("facets", 0) + len(list(sol))
    if "materials" in sec:
        s["materials"] = sum(1 for c in sec["materials"] if c.tag == "material")
    if "structure" in sec:
        for ref in sec["structure"].iter("materialref"):
            k = "materialref:" + ref.get("ref", "?")
            s[k] = s.get(k, 0) + 1
        for c in sec["structure"]:
            s["struct:" + c.tag] = s.get("struct:" + c.tag, 0) + 1
        s["physvol"] = sum(1 for _ in sec["structure"].iter("physvol"))
    if "define" in sec:
        s["defines"] = len(list(sec["define"]))
    return s


def diffSummary(a, b):
    out = []
    for k in sorted(set(a) | set(b)):
        if k in ("defines", "materials"):
            continue   # naming / G4_* definition policy; materialref:* carries the assignment
        if a.get(k) != b.get(k):
            out.append(f"{k}: {a.get(k, 0)} -> {b.get(k, 0)}")
    return out


def loadVersion(vdir):
    recs = {}
    with open(os.path.join(vdir, "results.jsonl"), encoding="utf-8") as f:
        for line in f:
            try:
                r = json.loads(line)
            except ValueError:
                continue
            recs[r["key"]] = r     # last one wins on resume
    env = {}
    envPath = os.path.join(vdir, "env.json")
    if os.path.isfile(envPath):
        with open(envPath, encoding="utf-8") as f:
            env = json.load(f)
    g4 = {}
    g4Path = os.path.join(vdir, "geant4.jsonl")
    if os.path.isfile(g4Path):
        with open(g4Path, encoding="utf-8") as f:
            for line in f:
                try:
                    r = json.loads(line)
                    g4[r["gdml"]] = r
                except ValueError:
                    pass
    return {"dir": vdir, "recs": recs, "env": env, "g4": g4, "done": os.path.isfile(os.path.join(vdir, "DONE"))}


def baselineFor(fcName):
    stem = fcName[:-6]
    cog = os.path.join(TC, "FC_COGs", stem + ".txt")
    cogRows = None
    if os.path.isfile(cog) and os.path.getsize(cog) > 0:
        with open(cog, encoding="utf-8", errors="replace") as f:
            cogRows = parseCOG(f.read())
    gd = sorted(glob.glob(os.path.join(TC, "FC_exported_gdml", glob.escape(stem) + "-*.gdml")))
    return cogRows, (gd[0] if gd else None)


def firstErrorLine(rec):
    if rec.get("missing_addons"):
        msgs = rec.get("missing_addon_messages") or []
        return ("[open] needs addon(s) not importable here: " + ", ".join(rec["missing_addons"])
                + (f" ({msgs[0][:150]})" if msgs else ""))
    if rec.get("status") == "CONSOLE_ERRORS" and rec.get("console_errors"):
        n = len(rec["console_errors"])
        return f"[console] {rec['console_errors'][0][:200]}" + (f" (+{n - 1} more)" if n > 1 else "")
    for st in ("open", "export", "cog", "import", "roundtrip"):
        e = rec.get(st, {}).get("error")
        if e:
            lines = [l for l in e.strip().splitlines() if l.strip()]
            return f"[{st}] {lines[-1][:220]}" if lines else f"[{st}] error"
    if rec.get("stage"):
        return f"[{rec['stage']}] {rec.get('status')} exit {rec.get('exit_code')}"
    return ""


def main(runDir):
    versions = {}
    for d in sorted(glob.glob(os.path.join(runDir, "*", "results.jsonl"))):
        vdir = os.path.dirname(d)
        versions[os.path.basename(vdir)] = loadVersion(vdir)
    if not versions:
        print("no results found under", runDir)
        return 1
    labels = list(versions)
    keys = sorted(set().union(*(v["recs"].keys() for v in versions.values())))

    report = {"run_dir": runDir, "versions": {}, "files": {}}
    for lab, v in versions.items():
        counts = {}
        for r in v["recs"].values():
            counts[r["status"]] = counts.get(r["status"], 0) + 1
        report["versions"][lab] = {"env": v["env"], "complete": v["done"], "counts": counts}

    for key in keys:
        kind, name = key.split(":", 1)
        entry = {"kind": kind, "file": name, "by_version": {}, "flags": []}
        baseCog, baseGdml = baselineFor(name) if kind == "fc" else (None, None)
        baseSum = gdmlSummary(baseGdml) if baseGdml else None
        cogs, sums = {}, {}
        for lab, v in versions.items():
            r = v["recs"].get(key)
            if r is None:
                entry["by_version"][lab] = {"status": "NOT_RUN"}
                continue
            e = {"status": r["status"], "error": firstErrorLine(r)}
            if r.get("open", {}).get("invalid_after_recompute"):
                e["invalid_after_recompute"] = r["open"]["invalid_after_recompute"]
            if r.get("dialogs"):
                e["dialogs"] = [f"{d.get('title')}: {d.get('text')}" for d in r["dialogs"]]
            logMsgs = []
            for st in ("export", "roundtrip"):
                logMsgs += ((r.get(st) or {}).get("log") or {}).get("messages") or []
            if logMsgs:
                e["export_log"] = logMsgs
            cogText = r.get("cog", {}).get("text")
            if cogText:
                cogs[lab] = parseCOG(cogText)
                if baseCog is not None:
                    e["cog_vs_baseline"] = compareCOG(baseCog, cogs[lab])
            gpath = r.get("export", {}).get("path") or r.get("roundtrip", {}).get("path")
            if gpath and r.get("export", r.get("roundtrip", {})).get("ok"):
                sums[lab] = gdmlSummary(os.path.join(v["dir"], gpath))
                if baseSum is not None:
                    e["gdml_vs_baseline"] = diffSummary(baseSum, sums[lab])
                g4 = v["g4"].get(os.path.basename(gpath))
                if g4:
                    e["geant4"] = g4.get("status")
            entry["by_version"][lab] = e

        statuses = {lab: e["status"] for lab, e in entry["by_version"].items()}
        if len(set(statuses.values())) > 1:
            entry["flags"].append("status differs between versions")
        if kind == "fc" and baseCog is not None:
            for lab, st in statuses.items():
                if not st.startswith("OK") and st != "NOT_RUN":
                    entry["flags"].append(f"REGRESSION in {lab}: baseline had COG output, now {st}")
        for lab, e in entry["by_version"].items():
            if e.get("cog_vs_baseline"):
                entry["flags"].append(f"COG differs from baseline in {lab}")
            if e.get("gdml_vs_baseline"):
                entry["flags"].append(f"GDML structure differs from baseline in {lab}")
            if e.get("invalid_after_recompute"):
                entry["flags"].append(f"invalid objects after recompute in {lab}")
            if e.get("export_log"):
                n = len(e["export_log"])
                entry["flags"].append(f"export log in {lab}: {e['export_log'][0][:200]}"
                                      + (f" (+{n - 1} more)" if n > 1 else ""))
            if e.get("geant4") and e["geant4"] != "OK":
                entry["flags"].append(f"Geant4 {e['geant4']} in {lab}")
        if len(cogs) > 1:
            first = labels[0] if labels[0] in cogs else next(iter(cogs))
            for lab in cogs:
                if lab != first:
                    d = compareCOG(cogs[first], cogs[lab])
                    if d:
                        entry["flags"].append(f"COG differs {first} vs {lab}")
                        entry["by_version"][lab]["cog_vs_" + first] = d
        if len(sums) > 1:
            first = next(iter(sums))
            for lab in sums:
                if lab != first:
                    d = diffSummary(sums[first], sums[lab])
                    if d:
                        entry["flags"].append(f"GDML structure differs {first} vs {lab}")
                        entry["by_version"][lab]["gdml_vs_" + first] = d
        report["files"][key] = entry

    with open(os.path.join(runDir, "report.json"), "w", encoding="utf-8") as f:
        json.dump(report, f, indent=1)
    md = renderMarkdown(report, labels)
    with open(os.path.join(runDir, "report.md"), "w", encoding="utf-8") as f:
        f.write(md)
    print(md.split("\n## Files needing attention")[0])
    return 0


def renderMarkdown(rep, labels):
    L = [f"# GDML autotest report\n", f"Run: `{rep['run_dir']}`\n", "## Versions\n"]
    L.append("| | " + " | ".join(labels) + " |")
    L.append("|---" * (len(labels) + 1) + "|")
    for field in ("freecad_version", "python", "occ", "gdml_commit", "gdml_xml_backend", "gdml_path", "machine"):
        L.append(f"| {field} | " + " | ".join(str(rep["versions"][l]["env"].get(field, "?")) for l in labels) + " |")
    L.append("| run complete | " + " | ".join("yes" if rep["versions"][l]["complete"] else "**no**" for l in labels) + " |")
    allStatus = sorted(set().union(*(rep["versions"][l]["counts"].keys() for l in labels)))
    for st in allStatus:
        L.append(f"| {st} | " + " | ".join(str(rep["versions"][l]["counts"].get(st, 0)) for l in labels) + " |")

    attention = [(k, e) for k, e in rep["files"].items()
                 if e["flags"] or any(not v["status"].startswith("OK") for v in e["by_version"].values())]
    regress = [(k, e) for k, e in attention if any(f.startswith("REGRESSION") for f in e["flags"])]
    L.append(f"\n**{len(rep['files'])} files tested, {len(attention)} need attention, {len(regress)} regressions vs baseline.**\n")

    L.append("## Status by file\n")
    L.append("| file | " + " | ".join(labels) + " |")
    L.append("|---" * (len(labels) + 1) + "|")
    for k, e in rep["files"].items():
        cells = []
        for l in labels:
            st = e["by_version"].get(l, {}).get("status", "NOT_RUN")
            cells.append(st if st.startswith("OK") else f"**{st}**")
        L.append(f"| {e['file']} | " + " | ".join(cells) + " |")

    L.append("\n## Files needing attention\n")
    for k, e in attention:
        L.append(f"### {e['file']}\n")
        for fl in dict.fromkeys(e["flags"]):
            L.append(f"- {fl}")
        for l in labels:
            v = e["by_version"].get(l)
            if not v:
                continue
            bits = []
            if v.get("error"):
                bits.append(f"error `{v['error']}`")
            for fld in [x for x in v if x.startswith(("cog_vs", "gdml_vs"))] + ["invalid_after_recompute", "dialogs"]:
                if v.get(fld):
                    items = v[fld][:6]
                    more = f" (+{len(v[fld]) - 6} more)" if len(v[fld]) > 6 else ""
                    bits.append(f"{fld}: " + "; ".join(str(i) for i in items) + more)
            if bits:
                L.append(f"- **{l}** ({v['status']}): " + " — ".join(bits))
        L.append("")
    return "\n".join(L) + "\n"


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(2)
    sys.exit(main(os.path.abspath(sys.argv[1])))
