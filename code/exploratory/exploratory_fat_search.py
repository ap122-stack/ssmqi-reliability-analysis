"""EXPLORATORY (not preregistered; written after the final registered run).
Part A (declared before running): does the share of Rater 2 images marked untraceable rise with the participant's
supraspinatus MRI mDixon fat fraction? Spearman across participants, participant bootstrap interval.
Part B (hypothesis-generating): does any of 15 declared grayscale/texture features track fat fraction better than
echo intensity? 15 features x 3 plane groups = 45 tests; the multiplicity is handled by a permutation test on the
MAXIMUM |Spearman rho| across all 45 (fat fraction shuffled across participants), so a feature is called out only if it
beats what the best of 45 would reach by chance. Reads the private clinical file; writes aggregate rows only.
Usage: python -I exploratory_fat_search.py
"""
import csv, os, random, statistics as st
ROOT = r"C:\Users\Owner\OneDrive\Rotator Cuff MVI\Reliability Study\Inter_Intra_Rater_Study"
FINAL = os.path.join(ROOT, "Analysis", "Registered", "FINAL_20261106_190031")
OUT = os.path.join(FINAL, "Exploratory")
rd = lambda p: list(csv.DictReader(open(p, encoding="utf-8-sig")))
clin = {r["Participant"]: float(r["mDixon_SSP_pct"]) for r in rd(os.path.join(ROOT, "Admin_Private", "Participant_Clinical_Minimal.csv")) if r["mDixon_SSP_pct"] != ""}

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
def boot_ci(x, y, B=2000):
    bs = []
    for _ in range(B):
        pick = [random.randrange(len(x)) for _ in x]
        b = spearman([x[i] for i in pick], [y[i] for i in pick])
        if b == b: bs.append(b)
    bs.sort(); return bs[int(0.025 * len(bs))], bs[int(0.975 * len(bs)) - 1]

random.seed(20261112)
# ---------- Part A: untraceable share vs fat
def untrace(rater, sess):
    folder = "Rater_0%d" % rater
    cl = rd(os.path.join(ROOT, folder, "Session_%d" % sess, "Session_%d_Case_List.csv" % sess))
    key = {k["Case_ID"]: k["Patient_Group_Inferred"] for k in rd(os.path.join(ROOT, "Admin_Private", "%s_Session_%d_Randomization_Key.csv" % (folder, sess)))}
    d = {}
    for c in cl:
        p = key[c["Case_ID"]]; d.setdefault(p, []).append(0 if c["Status"].strip().lower() == "complete" else 1)
    return d
rowsA = []
for label, parts in (("Rater 2 Session 2", [untrace(2, 2)]), ("Rater 2 both sessions", [untrace(2, 1), untrace(2, 2)])):
    share = {}
    for d in parts:
        for p, v in d.items(): share.setdefault(p, []).extend(v)
    ids = sorted(p for p in share if p in clin)
    x = [clin[p] for p in ids]; y = [sum(share[p]) / len(share[p]) for p in ids]
    rho = spearman(x, y); lo, hi = boot_ci(x, y)
    rowsA.append([label, len(ids), round(rho, 3), round(lo, 3), round(hi, 3)])
    print("A:", label, "n", len(ids), "rho %.3f (%.3f to %.3f)" % (rho, lo, hi))
with open(os.path.join(OUT, "Untraceable_share_vs_mDixon_SSP.csv"), "w", newline="") as f:
    w = csv.writer(f); w.writerow(["Comparison", "N_participants", "Spearman_rho", "CI_low", "CI_high"]); w.writerows(rowsA)

# ---------- Part B: feature search with max-statistic permutation
feat = [r for r in rd(os.path.join(FINAL, "Exploratory", "Tracing_Texture_Features.csv")) if r["Complete"].lower() in ("1", "true") and r["Session"] == "1"]
names = [k for k in feat[0] if k not in ("Rater", "Session", "Source_Image", "Participant", "Axis", "Complete")]
tests = []   # (plane, feature, participants, x-values)
for plane in ("Long axis", "Short axis", "Both planes"):
    for nm in names:
        per = {}
        for r in feat:
            if plane != "Both planes" and r["Axis"] != plane: continue
            if r["Participant"] in clin and r[nm] not in ("", "NaN"): per.setdefault(r["Participant"], []).append(float(r[nm]))
        ids = sorted(per); tests.append((plane, nm, ids, [st.mean(per[p]) for p in ids]))
def all_rho(shuffle):
    out = []
    for plane, nm, ids, x in tests:
        y = [clin[p] for p in ids]
        if shuffle: y = [clin[p] for p in shuffle(ids)]
        out.append(spearman(x, y))
    return out
obs = all_rho(None)
allids = sorted(clin); mx = []
for _ in range(4000):
    perm = dict(zip(allids, random.sample(allids, len(allids))))
    rr = []
    for plane, nm, ids, x in tests:
        rr.append(spearman(x, [clin[perm[p]] for p in ids]))
    mx.append(max(abs(v) for v in rr if v == v))
rowsB = []
for (plane, nm, ids, x), r in zip(tests, obs):
    adj = (1 + sum(m >= abs(r) for m in mx)) / (len(mx) + 1) if r == r else float("nan")
    rowsB.append([plane, nm, len(ids), round(r, 3), round(adj, 4)])
rowsB.sort(key=lambda z: -abs(z[3]))
with open(os.path.join(OUT, "Texture_features_vs_mDixon_SSP_search.csv"), "w", newline="") as f:
    w = csv.writer(f); w.writerow(["Plane", "Feature", "N_participants", "Spearman_rho", "FWER_adjusted_p_maxstat"]); w.writerows(rowsB)
print("B: 45 tests; chance-level max|rho| 95th percentile = %.3f" % sorted(mx)[int(0.95 * len(mx))])
print("top 10 by |rho|:"); [print("  ", r) for r in rowsB[:10]]
print("EI_mean rows:", [r for r in rowsB if r[1] == "EI_mean"])
