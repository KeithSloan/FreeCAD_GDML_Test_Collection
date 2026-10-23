#!/usr/bin/env python3
"""Compare FreeCAD centres of mass with Geant4's (load_gdml_color -printCM).

    autotest/geant4_cog.py <run_dir> [--vol-tol 0.01] [--cm-tol 0.01] [--inertia-tol 0.02]

For every version directory in <run_dir> that has both cogs/ (FreeCAD,
calcCenterOfMass.py) and geant4_out/ (Geant4, from run_geant4.sh), pairs each
test file's FreeCAD COG with the Geant4 output of its exported GDML and
compares the "Total" lines:

  volume   relative difference           <= --vol-tol      (default 1 %)
  cm       distance                      <= --cm-tol x Rg  (default 1 % of the
           radius of gyration), never tighter than the printed precision
  inertia  each diagonal term (about the origin), relative to the largest
                                         <= --inertia-tol  (default 2 %)

Geant4's figures are Monte Carlo estimates (10^6 points per volume by default),
so small differences are expected; the tolerances allow for that.

Per-volume lines are paired by name where both sides print the same name;
those pairs are reported for information but do not decide the status.

Writes <version>/geant4_cog.jsonl and <run_dir>/geant4_cog.md.
Status per file: MATCH, DIFF, NO_GEANT4 (no parsable Geant4 moments - see
geant4.jsonl for why), NO_FC (no FreeCAD total).
"""
import argparse
import glob
import json
import math
import os
import re
import sys

MARK = "***** geometry moments ****"

# Both sides print, per volume:
#   <name>  Volume:  V unit^3 --> center of "mass": (x, y, z) unit, moments of
#   inertia: (Ixx, Iyy, Izz) [x10^k] mm^5|cm^5 [Entity type: G4Box]
# (FreeCAD spells it "inetria").
NUM = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?|[-+]?nan|[-+]?inf"
LINE = re.compile(
    r"^(?P<name>\S+)\s+Volume:\s*(?P<v>" + NUM + r")\s*(?P<vu>m\^3|cm\^3|mm\^3)\s*"
    r"-->\s*center of \"mass\":\s*\(\s*(?P<x>" + NUM + r")\s*,\s*(?P<y>" + NUM + r")\s*,\s*(?P<z>" + NUM + r")\s*\)\s*"
    r"(?P<cu>m|cm|mm)\s*,\s*moments of in[a-z]+:\s*\(\s*(?P<ixx>" + NUM + r")\s*,\s*(?P<iyy>" + NUM + r")\s*,\s*(?P<izz>" + NUM + r")\s*\)\s*"
    r"(?:x10\^(?P<k>-?\d+)\s*)?(?P<iu>cm\^5|mm\^5)"
    r"(?:\s+Entity type:\s*(?P<etype>\S+))?",
    re.IGNORECASE)

LEN = {"mm": 1.0, "cm": 10.0, "m": 1000.0}
VOL = {"mm^3": 1.0, "cm^3": 1e3, "m^3": 1e9}
INER = {"mm^5": 1.0, "cm^5": 1e5}


def _decimals(s):
    s = s.lower()
    if "e" in s:
        mant, exp = s.split("e")
        d = len(mant.split(".")[1]) if "." in mant else 0
        return d - int(exp)
    return len(s.split(".")[1]) if "." in s else 0


def parse_line(line):
    """One volume line -> dict in mm / mm^3 / mm^5, or None."""
    m = LINE.search(line.strip())
    if not m:
        return None
    g = m.groupdict()
    lu = LEN[g["cu"].lower()]
    scale = INER[g["iu"].lower()] * (10.0 ** int(g["k"] or 0))
    vals = [float(g[c]) for c in ("v", "x", "y", "z", "ixx", "iyy", "izz")]
    if any(math.isnan(v) or math.isinf(v) for v in vals):
        return {"name": g["name"], "bad": True}
    # resolution of the printed centre of mass (half a unit in the last digit)
    res = max(0.5 * 10 ** -_decimals(g[c]) for c in ("x", "y", "z")) * lu
    return {
        "name": g["name"],
        "volume": vals[0] * VOL[g["vu"].lower()],
        "cm": [vals[1] * lu, vals[2] * lu, vals[3] * lu],
        "cm_res": res,
        "inertia": [vals[4] * scale, vals[5] * scale, vals[6] * scale],
        "etype": g.get("etype"),
    }


def parse_text(text, between_marks):
    """All volume lines in text -> list of dicts.  For Geant4 output only the
    block between the two 'geometry moments' markers is used."""
    lines = text.splitlines()
    if between_marks:
        idx = [i for i, l in enumerate(lines) if MARK in l]
        if len(idx) < 2:
            return []
        lines = lines[idx[0] + 1:idx[1]]
    out = []
    for l in lines:
        p = parse_line(l)
        if p:
            out.append(p)
    return out


def total(entries):
    for e in entries:
        if e["name"].lower() == "total":
            return e
    return None


def compare(fc, g4, vol_tol, cm_tol, in_tol):
    """Compare two parsed entries -> (ok, details dict)."""
    d = {}
    problems = []
    if fc.get("bad") or g4.get("bad"):
        return False, {"problems": ["nan/inf in output"]}
    dv = abs(g4["volume"] - fc["volume"]) / max(abs(fc["volume"]), 1e-12)
    d["volume_fc_mm3"] = fc["volume"]
    d["volume_g4_mm3"] = g4["volume"]
    d["volume_rel"] = dv
    if dv > vol_tol:
        problems.append(f"volume {fmt(fc['volume'])} vs {fmt(g4['volume'])} mm^3 ({100*dv:.2f} %)")
    # radius of gyration about the centre of mass, from the trace about the
    # origin: tr(I_origin) = 2*sum(r^2 dV); about cm subtract 2*V*|cm|^2
    V = fc["volume"]
    tr = sum(fc["inertia"])
    r2 = tr / (2 * V) - sum(c * c for c in fc["cm"]) if V > 0 else 0
    rg = math.sqrt(r2) if r2 > 0 else (V ** (1 / 3) if V > 0 else 1.0)
    dcm = math.dist(fc["cm"], g4["cm"])
    tol = max(cm_tol * rg, fc["cm_res"] + g4["cm_res"])
    d["cm_fc_mm"] = fc["cm"]
    d["cm_g4_mm"] = g4["cm"]
    d["cm_dist_mm"] = dcm
    d["cm_tol_mm"] = tol
    d["rg_mm"] = rg
    if dcm > tol:
        problems.append("cm (" + ", ".join(f"{c:.3g}" for c in fc["cm"]) + ") vs ("
                        + ", ".join(f"{c:.3g}" for c in g4["cm"]) + f") mm, off by {dcm:.3g} mm (tol {tol:.2g})")
    imax = max(max(abs(i) for i in fc["inertia"]), 1e-30)
    di = [abs(a - b) / imax for a, b in zip(fc["inertia"], g4["inertia"])]
    d["inertia_fc_mm5"] = fc["inertia"]
    d["inertia_g4_mm5"] = g4["inertia"]
    d["inertia_rel"] = max(di)
    if max(di) > in_tol:
        problems.append("inertia (" + ", ".join(fmt(i) for i in fc["inertia"]) + ") vs ("
                        + ", ".join(fmt(i) for i in g4["inertia"]) + f") mm^5 ({100*max(di):.1f} %)")
    d["problems"] = problems
    return not problems, d


def fmt(x):
    return f"{x:.4g}"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("run_dir")
    ap.add_argument("--vol-tol", type=float, default=0.01)
    ap.add_argument("--cm-tol", type=float, default=0.01)
    ap.add_argument("--inertia-tol", type=float, default=0.02)
    a = ap.parse_args()

    md = ["# FreeCAD vs Geant4 centre of mass", "",
          f"Run: `{os.path.abspath(a.run_dir)}`", "",
          f"Tolerances: volume {100*a.vol_tol:g} %, centre of mass {100*a.cm_tol:g} % of radius of gyration, "
          f"inertia {100*a.inertia_tol:g} % of largest term.  Geant4 values are Monte Carlo estimates.", ""]
    any_version = False
    for vdir in sorted(glob.glob(os.path.join(a.run_dir, "*", ""))):
        res_path = os.path.join(vdir, "results.jsonl")
        g4dir = os.path.join(vdir, "geant4_out")
        if not (os.path.isfile(res_path) and os.path.isdir(g4dir)):
            continue
        any_version = True
        ver = os.path.basename(os.path.normpath(vdir))
        g4status = {}
        p = os.path.join(vdir, "geant4.jsonl")
        if os.path.isfile(p):
            for l in open(p, encoding="utf-8"):
                if l.strip():
                    r = json.loads(l)
                    g4status[r["gdml"]] = r.get("status")
        rows = []
        for l in open(res_path, encoding="utf-8"):
            if not l.strip():
                continue
            r = json.loads(l)
            exp = (r.get("export") or {}).get("path")
            cog = (r.get("cog") or {}).get("path")
            if not exp or r.get("kind") != "fc":
                continue
            gname = os.path.basename(exp)
            row = {"file": r["file"], "gdml": gname, "geant4_load": g4status.get(gname)}
            fc_entries = []
            if cog and os.path.isfile(os.path.join(vdir, cog)):
                fc_entries = parse_text(open(os.path.join(vdir, cog), encoding="utf-8", errors="replace").read(), False)
            outp = os.path.join(g4dir, gname + ".out")
            g4_entries = []
            if os.path.isfile(outp):
                g4_entries = parse_text(open(outp, encoding="utf-8", errors="replace").read(), True)
            fct, g4t = total(fc_entries), total(g4_entries)
            if not g4t:
                row["status"] = "NO_GEANT4"
            elif not fct:
                row["status"] = "NO_FC"
            else:
                ok, det = compare(fct, g4t, a.vol_tol, a.cm_tol, a.inertia_tol)
                row.update(det)
                row["status"] = "MATCH" if ok else "DIFF"
                # informational per-volume pairs
                g4by = {e["name"]: e for e in g4_entries if e["name"].lower() != "total"}
                pairs = []
                for e in fc_entries:
                    if e["name"].lower() == "total" or e["name"] not in g4by:
                        continue
                    okv, dv = compare(e, g4by[e["name"]], a.vol_tol, a.cm_tol, a.inertia_tol)
                    pairs.append({"name": e["name"], "ok": okv, "problems": dv.get("problems", [])})
                row["volumes"] = pairs
                row["n_fc_volumes"] = len(fc_entries) - 1
                row["n_g4_volumes"] = len(g4_entries) - 1
            rows.append(row)
        with open(os.path.join(vdir, "geant4_cog.jsonl"), "w", encoding="utf-8") as f:
            for row in rows:
                f.write(json.dumps(row) + "\n")
        counts = {}
        for row in rows:
            counts[row["status"]] = counts.get(row["status"], 0) + 1
        md += [f"## {ver}", "", ", ".join(f"{k} {v}" for k, v in sorted(counts.items())), "",
               "| file | status | volume | cm off / tol (mm) | inertia | G4 load |",
               "|---|---|---|---|---|---|"]
        for row in sorted(rows, key=lambda r: (r["status"] == "MATCH", r["file"])):
            if "volume_rel" in row:
                vol = f"{100*row['volume_rel']:.2f} %"
                cm = f"{row['cm_dist_mm']:.3g} / {row['cm_tol_mm']:.2g}"
                ine = f"{100*row['inertia_rel']:.1f} %"
            else:
                vol = cm = ine = ""
            md.append(f"| {row['file']} | {row['status']} | {vol} | {cm} | {ine} | {row.get('geant4_load') or ''} |")
        md.append("")
        diffs = [r for r in rows if r["status"] == "DIFF"]
        if diffs:
            md += ["### Differences", ""]
            for row in diffs:
                md.append(f"- **{row['file']}** (`{row['gdml']}`, {row['n_fc_volumes']} FC / {row['n_g4_volumes']} G4 volumes)")
                for pr in row["problems"]:
                    md.append(f"  - {pr}")
                bad = [p for p in row.get("volumes", []) if not p["ok"]]
                for p in bad[:8]:
                    md.append(f"  - volume `{p['name']}`: " + "; ".join(p["problems"]))
            md.append("")
    if not any_version:
        print("geant4_cog: no version directory with results.jsonl and geant4_out/", file=sys.stderr)
        return 1
    with open(os.path.join(a.run_dir, "geant4_cog.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(md) + "\n")
    print("\n".join(md))
    return 0


if __name__ == "__main__":
    sys.exit(main())
