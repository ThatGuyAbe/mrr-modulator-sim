"""loaders for the matlab outputs (matlab/step*.m) and small helpers shared by
figures.py and tables.py. every number in the figures and tables comes from these
files. root is the out folder the matlab scripts write (out_smoke for a smoke test).
"""
import os
import numpy as np
import scipy.io

STEPS = {0: "step0_checks", 1: "step1_sweeps", 2: "step2_mechanism", 3: "step3_eyes",
         4: "step4_ring_dynamics", 5: "step5_convergence", 6: "step6_population", 7: "step7_tuning_error"}


def loadmat(path):
    return scipy.io.loadmat(path, squeeze_me=True, struct_as_record=False)


def path(root, step, name):
    return os.path.join(root, STEPS[step], name)


def have(root, step, name):
    return os.path.exists(path(root, step, name))


def sweeps(root):
    """step 1 sweeps and targets. arrays keep the matlab axis order:
    oma_sw[detuning, Q, power, rate], dw_opt[Q, power, rate]."""
    m = loadmat(path(root, 1, "results.mat"))
    d = {k: m[k] for k in ("df_grid", "hw_GHz", "oma_sw", "veo_sw", "ber_sw", "oma_lin", "veo_lin",
                           "dw_opt", "dw_veo", "dw_lin", "oma_pk", "veo_at_pk", "veo_pk", "oma_0", "veo_0")}
    d["Q_values"] = [int(q) for q in np.atleast_1d(m["Q_values"])]
    d["powers"] = [float(p) for p in np.atleast_1d(m["powers"])]
    d["datarates"] = [float(r) for r in np.atleast_1d(m["datarates"])]
    nQ, nP, nR = len(d["Q_values"]), len(d["powers"]), len(d["datarates"])
    for k in ("dw_opt", "dw_veo", "oma_pk", "veo_at_pk", "veo_pk", "oma_0", "veo_0"):
        d[k] = np.asarray(d[k], float).reshape(nQ, nP, nR)
    for k in ("oma_sw", "veo_sw", "ber_sw"):
        d[k] = np.asarray(d[k], float).reshape(len(d["df_grid"]), nQ, nP, nR)
    d["dw_lin"] = np.atleast_1d(d["dw_lin"]).astype(float)
    d["hw_GHz"] = np.atleast_1d(d["hw_GHz"]).astype(float)
    d["info"] = m.get("run_info")
    return d


def target(R1, Q, P, rate, which="dw_opt"):
    return float(R1[which][R1["Q_values"].index(Q), R1["powers"].index(float(P)), R1["datarates"].index(float(rate))])


def pop(root, Q):
    m = loadmat(path(root, 6, f"pop_results_Q{Q}.mat"))
    d = dict(m)
    d["configs"] = np.atleast_2d(m["configs"]).astype(float)
    for k in ("OMA_asf", "OMA_cor", "OMA_opt", "EO_asf", "EO_cor", "EO_opt", "BER_asf", "BER_cor", "BER_opt"):
        d[k] = np.asarray(m[k], float).reshape(int(m["N_rings"]), -1)
    for k in ("LV_asf", "LV_cor", "LV_opt"):
        d[k] = np.asarray(m[k], float).reshape(int(m["N_rings"]), 4, -1)
    d["asf_rings"] = np.atleast_1d(m["asf_rings"]).astype(int) - 1          # 0-based
    d["df_opt"] = np.atleast_1d(m["df_opt"]).astype(float)
    return d


def pop_qs(root):
    """the Q values that have a population file, in order."""
    return [q for q in (3000, 6000, 9000, 12000, 15000) if have(root, 6, f"pop_results_Q{q}.mat")]


def cfg_index(d, P, rate):
    c = d["configs"]
    hit = np.where((np.abs(c[:, 0] - P) < 1e-9) & (np.abs(c[:, 1] - rate) < 1e-9))[0]
    assert hit.size == 1, (P, rate)
    return int(hit[0])


def lber(b):
    return np.log10(np.maximum(b, 1e-30))


def case_arrays(d, ci, case):
    """oma, veo, ber of one case in one configuration. af uses only rings 1..n_asf."""
    key = {"AF": "asf", "ZDC": "cor", "PAC": "opt"}[case]
    rows = d["asf_rings"] if case == "AF" else slice(None)
    return d["OMA_" + key][rows, ci], d["EO_" + key][rows, ci], d["BER_" + key][rows, ci]


def stats(root, Qs=None):
    """population statistics. S[Q][(P, rate)][case] = dict(oma, veo, lber, inv, neg, n, inv_pct)."""
    out = {}
    for f in sorted(os.listdir(os.path.join(root, STEPS[6]))):
        if not (f.startswith("pop_results_Q") and f.endswith(".mat")):
            continue
        Q = int(f[len("pop_results_Q"):-4])
        if Qs and Q not in Qs:
            continue
        d = pop(root, Q)
        out[Q] = {}
        for ci, (P, rate) in enumerate(d["configs"]):
            out[Q][(float(P), float(rate))] = {}
            for case in ("AF", "ZDC", "PAC"):
                o, e, b = case_arrays(d, ci, case)
                assert np.all(np.isfinite(o)) and np.all(np.isfinite(e)), (Q, P, rate, case)
                out[Q][(float(P), float(rate))][case] = dict(
                    oma=float(o.mean()), veo=float(e.mean()), lber=float(lber(b).mean()), inv=int((e < 0).sum()),
                    neg=int((o < 0).sum()), n=int(o.size), inv_pct=100.0 * float((e < 0).mean()),
                    dw=float(d["df_opt"][ci]))
    return dict(sorted(out.items()))


def num(x, fmt=".2f"):
    """signed number for latex with a true minus sign."""
    s = format(x, fmt)
    return f"${s}$" if s.startswith("-") else s


def rate_txt(r):
    return f"{r:g}"
