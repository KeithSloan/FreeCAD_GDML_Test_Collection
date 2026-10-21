"""Record the file FreeCAD was working on when it crashed or was killed.

    python3 record_abort.py <version_out_dir> <reason> <exit_code>

Reads current.json (heartbeat written by the runner) and appends a result line
so the relaunched runner skips that file. Exit status 0 = a file was recorded
(or had already been recorded) so a relaunch will make progress; 1 = nothing
to record (the runner never reached a file), so relaunching would loop.
"""
import json
import os
import sys


def main():
    out, reason, rc = sys.argv[1], sys.argv[2], sys.argv[3]
    cur = os.path.join(out, "current.json")
    if not os.path.isfile(cur):
        return 1
    with open(cur, encoding="utf-8") as f:
        hb = json.load(f)
    name = hb.get("file")
    if not name:
        os.remove(cur)
        return 1
    kind = "gdml" if name.lower().endswith(".gdml") else "fc"
    key = f"{kind}:{name}"

    results = os.path.join(out, "results.jsonl")
    seen = False
    if os.path.isfile(results):
        with open(results, encoding="utf-8") as f:
            for line in f:
                try:
                    if json.loads(line).get("key") == key:
                        seen = True
                        break
                except ValueError:
                    pass
    if not seen:
        rec = {"key": key, "kind": kind, "file": name, "status": reason,
               "stage": hb.get("stage"), "exit_code": int(rc) if rc.lstrip("-").isdigit() else rc}
        with open(results, "a", encoding="utf-8") as f:
            f.write(json.dumps(rec) + "\n")
    os.remove(cur)
    print(f"recorded {reason} for {name} (stage {hb.get('stage')}, exit {rc})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
