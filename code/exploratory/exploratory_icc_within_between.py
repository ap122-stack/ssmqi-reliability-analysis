"""EXPLORATORY (post hoc; not preregistered; written after the final registered run).

The registered H2 was tested with one pooled ICC(A,1). This script refits the SAME model as the frozen
code (ssmqi_stats_core.fit_components: EI ~ Rater + (1|Participant) + (1|Image), REML, with
sigma2_rater = delta^2/2 and ICC(A,1) = (s2_participant + s2_image) / (that + s2_rater + s2_residual))
on subsets, using numpy only (no MATLAB):
  pooled check  : all accepted tracings per plane; must reproduce Table 4 (ICC 0.80 long, 0.95 short)
  within rater  : one rater, Session 1 vs Session 2 (the session plays the role of "rater" in the model)
  between raters: both raters within one session (primary = Session 1, secondary = Session 2)
Echo intensity: Exploratory\\Tracing_Texture_Features.csv, column EI_mean (same displayed image, artifact
mask and 0-255 rescale as the registered measure). 95% CIs: bootstrap resampling whole participants.
Usage: python -I exploratory_icc_within_between.py SEED NBOOT OUTCSV
"""
import csv, os, sys, math
import numpy as np

FINAL = r"C:\Users\Owner\OneDrive\Rotator Cuff MVI\Reliability Study\Inter_Intra_Rater_Study\Analysis\Registered\FINAL_20261106_190031\Exploratory"

def nelder_mead(f, x0, step=1.0, iters=400, tol=1e-9):
    n = len(x0); pts = [np.array(x0, float)] + [np.array(x0, float) + step * np.eye(n)[i] for i in range(n)]
    vals = [f(p) for p in pts]
    for _ in range(iters):
        o = np.argsort(vals); pts = [pts[i] for i in o]; vals = [vals[i] for i in o]
        if abs(vals[-1] - vals[0]) < tol: break
        c = np.mean(pts[:-1], axis=0)
        xr = c + (c - pts[-1]); fr = f(xr)
        if fr < vals[0]:
            xe = c + 2 * (c - pts[-1]); fe = f(xe)
            if fe < fr: pts[-1], vals[-1] = xe, fe
            else: pts[-1], vals[-1] = xr, fr
        elif fr < vals[-2]: pts[-1], vals[-1] = xr, fr
        else:
            xc = c + 0.5 * (pts[-1] - c); fc = f(xc)
            if fc < vals[-1]: pts[-1], vals[-1] = xc, fc
            else:
                pts = [pts[0]] + [pts[0] + 0.5 * (p - pts[0]) for p in pts[1:]]; vals = [vals[0]] + [f(p) for p in pts[1:]]
    o = int(np.argmin(vals)); return pts[o], vals[o]

def fit(y, fixed_codes, part_codes, img_codes, x0=None):
    """REML for y = X b + Zp u + Zi v + e. fixed_codes: 0/1 per row (two levels)."""
    n = len(y)
    X = np.column_stack([np.ones(n), fixed_codes.astype(float)]) if len(set(fixed_codes)) > 1 else np.ones((n, 1))
    Zp = np.zeros((n, part_codes.max() + 1)); Zp[np.arange(n), part_codes] = 1
    Zi = np.zeros((n, img_codes.max() + 1)); Zi[np.arange(n), img_codes] = 1
    A, B = Zp @ Zp.T, Zi @ Zi.T; I = np.eye(n)
    def nll(th):
        sp, si, se = np.exp(th)
        V = sp * A + si * B + se * I
        try: L = np.linalg.cholesky(V)
        except np.linalg.LinAlgError: return 1e12
        ld = 2 * np.log(np.diag(L)).sum()
        Vi_X = np.linalg.solve(V, X); XtViX = X.T @ Vi_X
        sign, ldx = np.linalg.slogdet(XtViX)
        if sign <= 0: return 1e12
        b = np.linalg.solve(XtViX, X.T @ np.linalg.solve(V, y)); r = y - X @ b
        return ld + ldx + r @ np.linalg.solve(V, r)
    v = max(np.var(y), 1e-6)
    th, _ = nelder_mead(nll, np.log([0.5 * v, 0.2 * v, 0.3 * v]) if x0 is None else x0, step=1.2 if x0 is None else 0.3, iters=400 if x0 is None else 120)
    sp, si, se = np.exp(th)
    V = sp * A + si * B + se * I; Vi_X = np.linalg.solve(V, X)
    b = np.linalg.solve(X.T @ Vi_X, X.T @ np.linalg.solve(V, y))
    delta = float(b[1]) if X.shape[1] > 1 else 0.0
    s2r = delta ** 2 / 2; num = sp + si
    return {"theta": th, "s2_part": sp, "s2_img": si, "s2_fixed": s2r, "s2_resid": se, "icc": num / (num + s2r + se)}

def codes(labels):
    u = {k: i for i, k in enumerate(dict.fromkeys(labels))}; return np.array([u[k] for k in labels])

def main(seed, nboot, outcsv):
    rows = list(csv.DictReader(open(os.path.join(FINAL, "Tracing_Texture_Features.csv"), encoding="utf-8-sig")))
    rows = [r for r in rows if r["Complete"] in ("1", "True", "true") and r["EI_mean"] not in ("", "NaN")]
    planes = ["Long axis", "Short axis"]
    specs = [("Pooled (registered model)", lambda r: True, "rater"),
             ("Within Rater 1 (S1 vs S2)", lambda r: r["Rater"] == "1", "session"),
             ("Within Rater 2 (S1 vs S2)", lambda r: r["Rater"] == "2", "session"),
             ("Between raters, Session 1", lambda r: r["Session"] == "1", "rater"),
             ("Between raters, Session 2", lambda r: r["Session"] == "2", "rater")]
    rng = np.random.default_rng(seed); out = []
    for pl in planes:
        for name, pred, fx in specs:
            sub = [r for r in rows if r["Axis"] == pl and pred(r)]
            y = np.array([float(r["EI_mean"]) for r in sub])
            f = codes([r["Rater"] if fx == "rater" else r["Session"] for r in sub])
            pc = codes([r["Participant"] for r in sub]); ic = codes([r["Source_Image"] for r in sub])
            est = fit(y, f, pc, ic)
            parts = sorted(set(r["Participant"] for r in sub)); idx = {p: [i for i, r in enumerate(sub) if r["Participant"] == p] for p in parts}
            bs = []
            for _ in range(nboot):
                pick = rng.integers(0, len(parts), len(parts)); ii = []; pl_ = []; im_ = []
                for q, k in enumerate(pick):
                    for i in idx[parts[k]]: ii.append(i); pl_.append("b%d" % q); im_.append(sub[i]["Source_Image"] + "_b%d" % q)
                try: bs.append(fit(y[ii], f[ii], codes(pl_), codes(im_), x0=est["theta"])["icc"])
                except Exception: pass
            lo, hi = np.percentile(bs, [2.5, 97.5]) if bs else (float("nan"),) * 2
            out.append([pl, name, "Rater" if fx == "rater" else "Session", len(sub), len(parts), round(est["icc"], 3), round(lo, 3), round(hi, 3),
                        round(est["s2_part"], 1), round(est["s2_img"], 1), round(est["s2_fixed"], 1), round(est["s2_resid"], 1)])
            print(out[-1], flush=True)
    with open(outcsv, "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["Plane", "Subset", "Fixed_Effect", "Tracings", "Participants", "ICC_A1", "CI_Low", "CI_High", "Var_Participant", "Var_Image", "Var_Fixed", "Var_Residual"]); w.writerows(out)

if __name__ == "__main__":
    main(int(sys.argv[1]), int(sys.argv[2]), sys.argv[3])
