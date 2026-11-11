"""Turn a FINAL registered-analysis folder into a mechanical Results draft.

Usage:  python -I build_results_draft_from_final.py <FINAL_folder> [output.md]

It reads Tables 2 to 6 and writes (a) formatted tables and (b) plain factual sentences
with the numbers filled in, plus a pass/fail check of each preregistered threshold
(H1: inter-rater Dice >= 0.70; H2: ICC(A,1) >= 0.75; H3: long-axis disagreement larger
than short-axis). It makes no interpretation and changes no analysis; the author writes
the Discussion. Nothing is hidden: every row of every table is carried through.
"""
import csv, sys, os

PRIMARY = "Inter-rater, primary (R1 S1 vs R2 S1)"
ORDER = [PRIMARY,
         "Inter-rater, secondary (R1 S2 vs R2 S2)",
         "Intra-rater, Rater 1 (S1 vs S2)",
         "Intra-rater, Rater 2 (S1 vs S2)",
         "Sensitivity (R1 S2 vs R2 S1)"]


def load(folder, name):
    with open(os.path.join(folder, name), newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def num(x, d=2):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return "NA"
    if v != v:
        return "NA"
    return f"{v:.{d}f}"


def ci(row, k, d=2):
    return f"{num(row[k], d)} ({num(row[k + '_CI_Low'], d)} to {num(row[k + '_CI_High'], d)})"


def md_table(headers, rows):
    out = ["| " + " | ".join(headers) + " |", "|" + "|".join("---" for _ in headers) + "|"]
    out += ["| " + " | ".join(r) + " |" for r in rows]
    return "\n".join(out)


def sort_rows(rows):
    return sorted(rows, key=lambda r: (ORDER.index(r["Comparison"]) if r["Comparison"] in ORDER else 99, r["Plane"]))


def main():
    folder = sys.argv[1]
    outp = sys.argv[2] if len(sys.argv) > 2 else os.path.join(folder, "Results_Draft_Auto.md")
    t2, t3 = load(folder, "Table2_Outline_Agreement.csv"), load(folder, "Table3_EchoIntensity_Agreement.csv")
    t4, t5 = load(folder, "Table4_Variance_Components_ICC.csv"), load(folder, "Table5_Decision_Study.csv")
    t6 = load(folder, "Table6_Traceability_Agreement.csv")
    L = ["# Results draft (auto-generated from " + os.path.basename(folder.rstrip('/\\')) + ")", "",
         "Mechanical numbers only. Edit wording; do not change any number. Report every row.", ""]

    L += ["## Table 2. Outline agreement", "",
          md_table(["Comparison", "Plane", "Pairs", "Dice, mean (95% CI)", "Jaccard, mean (95% CI)",
                    "HD95 mm, mean (95% CI)", "MSD mm, mean (95% CI)"],
                   [[r["Comparison"], r["Plane"], r["Pairs"],
                     f'{num(r["Dice_Mean"])} ({num(r["Dice_CI_Low"])} to {num(r["Dice_CI_High"])})',
                     f'{num(r["Jaccard_Mean"])} ({num(r["Jaccard_CI_Low"])} to {num(r["Jaccard_CI_High"])})',
                     f'{num(r["HD95_mm_Mean"])} ({num(r["HD95_mm_CI_Low"])} to {num(r["HD95_mm_CI_High"])})',
                     f'{num(r["MSD_mm_Mean"])} ({num(r["MSD_mm_CI_Low"])} to {num(r["MSD_mm_CI_High"])})']
                    for r in sort_rows(t2)]), ""]

    L += ["## Table 3. Echo-intensity agreement", "",
          md_table(["Comparison", "Plane", "Pairs", "Bias % (95% CI)", "Lower LoA %", "Upper LoA %",
                    "Prop. bias p", "SEM (95% CI)", "MDC95 (95% CI)", "CV % (95% CI)"],
                   [[r["Comparison"], r["Plane"], r["Pairs"],
                     f'{num(r["Bias_Percent"])} ({num(r["Bias_CI_Low"])} to {num(r["Bias_CI_High"])})',
                     f'{num(r["LoA_Lower"])} ({num(r["LoA_Lower_CI_Low"])} to {num(r["LoA_Lower_CI_High"])})',
                     f'{num(r["LoA_Upper"])} ({num(r["LoA_Upper_CI_Low"])} to {num(r["LoA_Upper_CI_High"])})',
                     num(r["ProportionalBias_p"], 3),
                     f'{num(r["SEM"])} ({num(r["SEM_CI_Low"])} to {num(r["SEM_CI_High"])})',
                     f'{num(r["MDC95"])} ({num(r["MDC95_CI_Low"])} to {num(r["MDC95_CI_High"])})',
                     f'{num(r["CV_Percent"])} ({num(r["CV_CI_Low"])} to {num(r["CV_CI_High"])})']
                    for r in sort_rows(t3)]), ""]

    L += ["## Table 4. Variance components and ICC(A,1)", "",
          md_table(["Plane", "Tracings", "Participants", "Participant var", "Image var", "Rater var",
                    "Residual var", "ICC(A,1) (95% CI)"],
                   [[r["Plane"], r["Tracings"], r["Participants"], num(r["Participant_Var"], 1), num(r["Image_Var"], 1),
                     num(r["Rater_Var"], 1), num(r["Residual_Var"], 1),
                     f'{num(r["ICC_A1"])} ({num(r["ICC_CI_Low"])} to {num(r["ICC_CI_High"])})'] for r in t4]), ""]

    L += ["## Table 5. Decision study (projected ICC)", "",
          md_table(["Plane", "Design", "Projected ICC (95% CI)"],
                   [[r["Plane"], r["Design"], f'{num(r["Projected_ICC"])} ({num(r["CI_Low"])} to {num(r["CI_High"])})']
                    for r in t5]), ""]

    L += ["## Table 6. Traceability agreement", "",
          md_table(["Comparison", "Plane", "Images", "% agreement (95% CI)", "Kappa (95% CI)", "PABAK (95% CI)"],
                   [[r["Comparison"], r["Plane"], r["Images"],
                     f'{num(r["PercentAgreement"], 1)} ({num(r["PA_CI_Low"], 1)} to {num(r["PA_CI_High"], 1)})',
                     f'{num(r["Kappa"])} ({num(r["Kappa_CI_Low"])} to {num(r["Kappa_CI_High"])})',
                     f'{num(r["PABAK"])} ({num(r["PABAK_CI_Low"])} to {num(r["PABAK_CI_High"])})']
                    for r in sort_rows(t6)]), ""]

    # Preregistered threshold checks (mechanical)
    L += ["## Preregistered threshold checks", ""]
    prim = {r["Plane"]: r for r in t2 if r["Comparison"] == PRIMARY}
    for plane, r in prim.items():
        lo, m = float(r["Dice_CI_Low"]), float(r["Dice_Mean"])
        L.append(f"- H1 ({plane}): primary inter-rater Dice {num(m)} (95% CI {num(r['Dice_CI_Low'])} to "
                 f"{num(r['Dice_CI_High'])}); point estimate {'meets' if m >= 0.70 else 'is below'} 0.70; "
                 f"lower CI bound {'is at or above' if lo >= 0.70 else 'is below'} 0.70.")
    for r in t4:
        lo, m = float(r["ICC_CI_Low"]), float(r["ICC_A1"])
        L.append(f"- H2 ({r['Plane']}): ICC(A,1) {num(m)} (95% CI {num(r['ICC_CI_Low'])} to {num(r['ICC_CI_High'])}); "
                 f"point estimate {'meets' if m >= 0.75 else 'is below'} 0.75; lower CI bound "
                 f"{'is at or above' if lo >= 0.75 else 'is below'} 0.75.")
    if "Long axis" in prim and "Short axis" in prim:
        a, b = float(prim["Long axis"]["Dice_Mean"]), float(prim["Short axis"]["Dice_Mean"])
        L.append(f"- H3 (exploratory): primary inter-rater Dice was {'lower' if a < b else 'higher'} in long-axis "
                 f"({num(a)}) than in short-axis ({num(b)}) images; difference {num(a - b)}.")
    L.append("")

    with open(outp, "w", encoding="utf-8") as f:
        f.write("\n".join(L))
    print("wrote", outp)


if __name__ == "__main__":
    main()
