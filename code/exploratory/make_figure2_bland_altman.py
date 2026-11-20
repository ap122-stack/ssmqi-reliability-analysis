"""Redraw manuscript Figure 2 (primary inter-rater Bland-Altman) from Pair_Level_Metrics.csv.
Usage: python -I make_figure2_bland_altman.py PAIR_CSV OUT_PNG
Same statistics as the registered analysis (percent difference, R2 minus R1, over pair mean).
Only the bias label position differs from the earlier export."""
import csv, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

pair_csv, out_png = sys.argv[1], sys.argv[2]
rows = list(csv.DictReader(open(pair_csv, encoding="utf-8")))
BLUE = "#1f5fae"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 13})
fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.9), dpi=200)
for ax, plane, letter in zip(axes, ["Long axis", "Short axis"], "AB"):
    s = [r for r in rows if r["Comparison"].startswith("Inter-rater, primary")
         and r["Plane"] == plane and r["EI_First"] and r["EI_Second"]]
    a = np.array([float(r["EI_First"]) for r in s]); b = np.array([float(r["EI_Second"]) for r in s])
    m = (a + b) / 2; d = 100 * (b - a) / m
    bias, sd = d.mean(), d.std(ddof=1)
    lo, hi = bias - 1.96 * sd, bias + 1.96 * sd
    ax.scatter(m, d, s=34, color=BLUE, zorder=3)
    ax.axhline(bias, color="0.25", lw=1.2, zorder=2)
    ax.axhline(lo, color="0.5", lw=1.1, ls="--"); ax.axhline(hi, color="0.5", lw=1.1, ls="--")
    # label sits in the right margin, outside the plotted data, clear of the y-axis tick labels
    ax.annotate(f"Bias {bias:+.1f}%", xy=(1.0, bias), xycoords=("axes fraction", "data"),
                xytext=(-4, 4), textcoords="offset points", ha="right", va="bottom", fontsize=11, color="0.15")
    ax.set_title(f"{plane} (n = {len(s)} pairs)", fontsize=14)
    ax.set_xlabel("Mean echo intensity (0 to 255)")
    ax.set_ylabel("Difference (% of pair mean)")
    ax.text(0.03, 0.95, letter, transform=ax.transAxes, fontsize=24, fontweight="bold", va="top")
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    print(plane, len(s), round(bias, 2), round(lo, 2), round(hi, 2))
fig.tight_layout()
fig.savefig(out_png, dpi=200, facecolor="white")
