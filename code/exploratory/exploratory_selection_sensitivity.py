"""EXPLORATORY (not preregistered) selection-sensitivity summaries from the FINAL registered output.

Reads only Pair_Level_Metrics.csv (written by the locked final run) and the session case lists / keys.
It recomputes nothing in the registered pipeline and changes no file. Output is printed and written
to Exploratory_Selection_Sensitivity.csv next to the FINAL tables. Label every number from this file
"exploratory" in the manuscript.

  1. share of usable pairs with Dice >= 0.70
  2. worst-case mean Dice: a pair with an untraceable outline counts as Dice 0 (denominator = all matched images)
  3. common subset: images accepted in all four sessions (n, mean Dice for the primary comparison)

Usage: python -I exploratory_selection_sensitivity.py
"""
import csv, os, random, statistics
ROOT = r"C:\Users\Owner\OneDrive\Rotator Cuff MVI\Reliability Study\Inter_Intra_Rater_Study"
FINAL = os.path.join(ROOT, "Analysis", "Registered", "FINAL_20261106_190031", "Primary_AllTracings")
IMAGES = {"Long axis": 57, "Short axis": 23}   # matched images per plane (Table 6)
random.seed(20261107)

pairs = list(csv.DictReader(open(os.path.join(FINAL, "Pair_Level_Metrics.csv"), encoding="utf-8")))
comps = []
for r in pairs:
    if r["Comparison"] not in comps: comps.append(r["Comparison"])

def boot_ci(vals_by_part, fn, B=2000):
    parts = list(vals_by_part)
    out = []
    for _ in range(B):
        pick = [random.choice(parts) for _ in parts]
        v = [x for p in pick for x in vals_by_part[p]]
        out.append(fn(v))
    out.sort()
    return out[int(0.025 * B)], out[int(0.975 * B)]

rows = []
for c in comps:
    for pl in ("Long axis", "Short axis"):
        g = [r for r in pairs if r["Comparison"] == c and r["Plane"] == pl]
        if not g: continue
        d = [float(r["Dice"]) for r in g]
        n, N = len(d), IMAGES[pl]
        share = sum(x >= 0.70 for x in d) / n
        worst = sum(d) / N
        by = {}
        for r in g: by.setdefault(r["Participant"], []).append(float(r["Dice"]))
        lo, hi = boot_ci({p: v for p, v in by.items()}, lambda v: sum(x >= 0.70 for x in v) / len(v))
        rows.append([c, pl, n, N, round(100 * share, 1), round(100 * lo, 1), round(100 * hi, 1), round(statistics.mean(d), 3), round(worst, 3)])

# common subset: accepted in all four sessions
def accepted(rater, sess):
    folder = f"Rater_0{rater}"
    cl = list(csv.DictReader(open(os.path.join(ROOT, folder, f"Session_{sess}", f"Session_{sess}_Case_List.csv"), encoding="utf-8-sig")))
    key = {k["Case_ID"]: k["Source_Image"] for k in csv.DictReader(open(os.path.join(ROOT, "Admin_Private", f"{folder}_Session_{sess}_Randomization_Key.csv"), encoding="utf-8-sig"))}
    return {key[c["Case_ID"]] for c in cl if c["Status"].strip().lower() == "complete"}
common = accepted(1, 1) & accepted(1, 2) & accepted(2, 1) & accepted(2, 2)
print("images accepted in all four sessions:", len(common))
prim = [r for r in pairs if r["Comparison"].startswith("Inter-rater, primary") and r["Source_Image"] in common]
for pl in ("Long axis", "Short axis"):
    g = [float(r["Dice"]) for r in prim if r["Plane"] == pl]
    if g: print(f"common subset, primary inter-rater, {pl}: n={len(g)} mean Dice={statistics.mean(g):.3f}")

out = os.path.join(FINAL, "Exploratory_Selection_Sensitivity.csv")
with open(out, "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["Comparison", "Plane", "Usable_pairs", "Matched_images", "Pct_pairs_Dice_ge_0.70", "Pct_CI_low", "Pct_CI_high", "Mean_Dice_usable", "Worst_case_mean_Dice"])
    w.writerows(rows)
for r in rows: print(r)
print("images in common subset:", len(common), "| wrote", out)
