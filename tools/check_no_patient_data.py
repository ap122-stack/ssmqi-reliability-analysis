"""Block commits that contain participant-level or private material.

Usage:  python tools/check_no_patient_data.py            (checks staged files in the current git repo)
        python tools/check_no_patient_data.py PATH ...   (checks the given files or folders)
Exit code 1 if anything is blocked.
"""
import os, re, subprocess, sys

BLOCK_EXT = {".dcm", ".mat", ".dpapi", ".tif", ".tiff", ".png", ".jpg", ".jpeg", ".nvp", ".nvpx", ".qdpx", ".mp4", ".m4a", ".wav"}
BLOCK_NAME = re.compile(r"(^|[\\/])(Admin_Private|Patient Results|Fresh Ultrasound Data|NORC|Blinded_Images|Trace_Results)([\\/]|$)|Randomization_Key|Participant_Clinical|Manifest|\.env($|\.)|secrets|credentials", re.I)
BLOCK_TEXT = [
    (re.compile(r"\bSWE0?\d{2,3}[_A-Z]"), "participant-style case code"),
    (re.compile(r"\bPT_\d{3}\b"), "patient code"),
    (re.compile(r"\b(Source_Image|Case_ID)\b.*\.dcm", re.I), "image file name"),
    (re.compile(r"AKIA[0-9A-Z]{16}|ghp_[A-Za-z0-9]{30,}|-----BEGIN [A-Z ]*PRIVATE KEY-----"), "credential"),
]
MAX_BYTES = 5 * 1024 * 1024

def files_from_args(args):
    if args:
        out = []
        for a in args:
            if os.path.isdir(a):
                for root, dirs, fs in os.walk(a):
                    dirs[:] = [d for d in dirs if d != ".git"]
                    out += [os.path.join(root, f) for f in fs]
            else:
                out.append(a)
        return out
    r = subprocess.run(["git", "diff", "--cached", "--name-only", "--diff-filter=ACM"], capture_output=True, text=True)
    return [l for l in r.stdout.splitlines() if l]

def main():
    bad = []
    for f in files_from_args(sys.argv[1:]):
        if os.path.basename(f) == "check_no_patient_data.py": continue
        ext = os.path.splitext(f)[1].lower()
        if ext in BLOCK_EXT: bad.append((f, "blocked file type " + ext)); continue
        if BLOCK_NAME.search(f): bad.append((f, "blocked path or name")); continue
        try:
            if os.path.getsize(f) > MAX_BYTES: bad.append((f, "larger than 5 MB")); continue
        except OSError:
            continue
        if ext in {".docx", ".pdf", ".xlsx", ".zip"}:
            continue
        try:
            txt = open(f, encoding="utf-8", errors="ignore").read()
        except Exception:
            continue
        for rx, why in BLOCK_TEXT:
            m = rx.search(txt)
            if m: bad.append((f, why + ": " + m.group(0)[:40])); break
    if bad:
        print("BLOCKED (%d):" % len(bad))
        for f, why in bad: print("  %s  <- %s" % (f, why))
        sys.exit(1)
    print("OK: no blocked files")

if __name__ == "__main__":
    main()
