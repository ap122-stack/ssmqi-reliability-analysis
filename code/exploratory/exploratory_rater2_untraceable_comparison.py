"""EXPLORATORY: where do Rater 2's untraceable calls fall across sessions, and how did Rater 1 treat the same images?
Reads only case lists and keys; prints counts, no tracing values. Usage: python -I exploratory_rater2_untraceable_comparison.py
"""
import csv, collections, os
ROOT = r"C:\Users\Owner\OneDrive\Rotator Cuff MVI\Reliability Study\Inter_Intra_Rater_Study"
def load(r, s):
    cl = list(csv.DictReader(open(os.path.join(ROOT, f"Rater_0{r}", f"Session_{s}", f"Session_{s}_Case_List.csv"), encoding="utf-8-sig")))
    key = {k["Case_ID"]: k["Source_Image"] for k in csv.DictReader(open(os.path.join(ROOT, "Admin_Private", f"Rater_0{r}_Session_{s}_Randomization_Key.csv"), encoding="utf-8-sig"))}
    return {key[c["Case_ID"]]: (c["Status"].strip().lower(), c["Axis"]) for c in cl}
r2s1, r2s2, r1s1, r1s2 = load(2, 1), load(2, 2), load(1, 1), load(1, 2)
u2 = [i for i, (s, a) in r2s2.items() if s != "complete"]; u1 = [i for i, (s, a) in r2s1.items() if s != "complete"]
print("R2 S2 untraceable:", len(u2), "| R2 S1 untraceable:", len(u1))
print("of the R2 S2 untraceable, status in R2 S1:", dict(collections.Counter(r2s1[i][0] for i in u2)))
print("of the R2 S1 untraceable, status in R2 S2:", dict(collections.Counter(r2s2[i][0] for i in u1)))
print("untraceable in both R2 sessions:", len(set(u1) & set(u2)))
print("plane of the R2 S2 untraceable:", dict(collections.Counter(r2s2[i][1] for i in u2)), "| all images:", dict(collections.Counter(a for s, a in r2s2.values())))
print("Rater 1 on the R2 S2 untraceable (S1, S2):", dict(collections.Counter(r1s1[i][0] for i in u2)), dict(collections.Counter(r1s2[i][0] for i in u2)))
