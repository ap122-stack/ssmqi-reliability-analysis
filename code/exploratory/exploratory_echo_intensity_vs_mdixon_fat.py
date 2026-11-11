"""EXPLORATORY (not preregistered; written after the final registered run).
Question: is a participant's mean relative echo intensity (rater-averaged, primary inter-rater pairs) associated with the
supraspinatus MRI mDixon fat fraction? Per plane, Spearman rho across participants, percentile bootstrap over participants.
Also prints aggregate cohort descriptives. Reads the private clinical file; writes aggregate rows only (no participant rows).
Usage: python -I exploratory_echo_intensity_vs_mdixon_fat.py
"""
import csv, os, random, statistics as st
ROOT = r"C:\Users\Owner\OneDrive\Rotator Cuff MVI\Reliability Study\Inter_Intra_Rater_Study"
FINAL = os.path.join(ROOT, "Analysis", "Registered", "FINAL_20261106_190031")
clin = {r["Participant"]: r for r in csv.DictReader(open(os.path.join(ROOT, "Admin_Private", "Participant_Clinical_Minimal.csv"), encoding="utf-8-sig"))}
pairs = [r for r in csv.DictReader(open(os.path.join(FINAL, "Primary_AllTracings", "Pair_Level_Metrics.csv"), encoding="utf-8-sig"))
         if r["Comparison"].startswith("Inter-rater, primary")]

def ranks(v):
    order = sorted(range(len(v)), key=lambda i: v[i]); r = [0.0] * len(v); i = 0
    while i < len(v):
        j = i
        while j + 1 < len(v) and v[order[j + 1]] == v[order[i]]: j += 1
        for k in range(i, j + 1): r[order[k]] = (i + j) / 2 + 1
        i = j + 1
    return r
def spearman(x, y):
    if len(x) < 4: return float("nan")
    rx, ry = ranks(x), ranks(y); mx, my = st.mean(rx), st.mean(ry)
    sx = sum((a - mx) ** 2 for a in rx); sy = sum((b - my) ** 2 for b in ry)
    return sum((a - mx) * (b - my) for a, b in zip(rx, ry)) / (sx * sy) ** 0.5 if sx > 0 and sy > 0 else float("nan")

random.seed(20261112)
rows = []
for plane in ("Long axis", "Short axis", "Both planes"):
    per = {}
    for r in pairs:
        if plane != "Both planes" and r["Plane"] != plane: continue
        p = r["Participant"]
        if p in clin and clin[p]["mDixon_SSP_pct"] != "":
            per.setdefault(p, []).append((float(r["EI_First"]) + float(r["EI_Second"])) / 2)
    ids = sorted(per)
    ei = [st.mean(per[p]) for p in ids]; fat = [float(clin[p]["mDixon_SSP_pct"]) for p in ids]
    rho = spearman(ei, fat); bs = []
    for _ in range(2000):
        pick = [random.randrange(len(ids)) for _ in ids]
        b = spearman([ei[i] for i in pick], [fat[i] for i in pick])
        if b == b: bs.append(b)
    bs.sort(); lo, hi = bs[int(0.025 * len(bs))], bs[int(0.975 * len(bs)) - 1]
    rows.append([plane, len(ids), round(rho, 3), round(lo, 3), round(hi, 3)])
    print(plane, "n participants", len(ids), "rho %.3f (95%% CI %.3f to %.3f)" % (rho, lo, hi))
outdir = os.path.join(FINAL, "Exploratory"); os.makedirs(outdir, exist_ok=True)
with open(os.path.join(outdir, "EI_vs_mDixon_SSP_Spearman.csv"), "w", newline="") as f:
    w = csv.writer(f); w.writerow(["Plane", "N_participants", "Spearman_rho", "CI_low", "CI_high"]); w.writerows(rows)
# aggregate descriptives for the 17 participants in the image set
ages = [float(r["Age"]) for r in clin.values()]; fats = [float(r["mDixon_SSP_pct"]) for r in clin.values() if r["mDixon_SSP_pct"] != ""]
print("participants", len(clin), "| female", sum(r["Sex"] == "F" for r in clin.values()), "male", sum(r["Sex"] == "M" for r in clin.values()),
      "| controls", sum(r["Group"] == "control" for r in clin.values()), "cases", sum(r["Group"] == "case" for r in clin.values()))
print("age mean %.1f SD %.1f median %.0f range %.0f-%.0f" % (st.mean(ages), st.stdev(ages), st.median(ages), min(ages), max(ages)))
print("mDixon SSP n=%d median %.1f IQR-ish min %.1f max %.1f mean %.1f" % (len(fats), st.median(fats), min(fats), max(fats), st.mean(fats)))
