"""thesis tables of the microring part (latex, booktabs), from the matlab outputs.

  table_dwopt.tex           power-aware detuning and low-power optimum of the nominal ring (step 1)
  table_linear_limit.tex    1 uW linear-limit check (step 5)
  table_convergence.tex     time step, warm-up and data window (step 5)
  table_bistable.tex        bistable bands of the nominal ring at 4 mW (step 4)
  table_population.tex      mean oma, veo, log10(ber) and invalid rings, every configuration (steps 1, 6)
  table_receiver_floor.tex  ranking of zdc and pac for three receiver floors (step 6)
  table_tuning_power.tex    tuning power per ring for zdc and pac (steps 1, 6)

each table is made only when its inputs exist. the numbers are also printed for checking.
run:  python3 tables.py <out folder of matlab> [table folder, default <out>/tables]
"""
import os
import sys
import numpy as np
from scipy.special import erfc
import mrr_data as D
from mrr_data import num


def MINUS(s):
    return s.replace("-", "$-$") if s.startswith("-") else s


def band(i):
    return r"\rowcolor{band}" if i % 2 else ""


def write(out, name, body):
    open(os.path.join(out, name), "w").write(body)
    print(f"--- {name}\n{body}")


def table_dwopt(root, out):
    R = D.sweeps(root)
    Qs, Ps, Rs = R["Q_values"], R["powers"], R["datarates"]
    body = []
    for iq, q in enumerate(Qs):
        for ip, p in enumerate(Ps):
            qc = rf"\multirow{{{len(Ps)}}}{{*}}{{{q}}}" if ip == 0 else ""
            lc = rf"\multirow{{{len(Ps)}}}{{*}}{{{num(R['dw_lin'][iq])}}}" if ip == 0 else ""
            vals = " & ".join(num(R["dw_opt"][iq, ip, ir]) for ir in range(len(Rs)))
            body.append(f"{qc} & {lc} & {p:g} & {vals} \\\\")
        if iq < len(Qs) - 1:
            body.append(r"\midrule")
    ncol = 3 + len(Rs)
    tex = (r"""\begin{table}[!t]
\centering
\caption{Power-Aware Detuning $\Delta\omega_{\mathrm{opt}}/2\pi$ and Low-Power Optimum $\Delta\omega_{\mathrm{lin}}/2\pi$ of the Nominal Ring (GHz)}
\label{tab:dwopt}
\setlength{\tabcolsep}{5pt}
\begin{tabular}{@{}""" + "c" * ncol + r"""@{}}
\toprule
 & & & \multicolumn{""" + str(len(Rs)) + r"""}{c}{$\Delta\omega_{\mathrm{opt}}/2\pi$ at bit rate (Gb/s)} \\
\cmidrule(l){4-""" + str(ncol) + r"""}
$Q$ & $\Delta\omega_{\mathrm{lin}}/2\pi$ & $P_{\mathrm{in}}$ (mW) & """ + " & ".join(D.rate_txt(r) for r in Rs) + r""" \\
\midrule
""" + "\n".join(body) + r"""
\bottomrule
\end{tabular}
\\[3pt]
\parbox{\columnwidth}{\footnotesize $\Delta\omega_{\mathrm{lin}}$ is the detuning of largest OMA of the nominal ring at 0.01~mW and 2.5~Gb/s. A negative detuning places the laser on the red side of the cold resonance.}
\end{table}
""")
    write(out, "table_dwopt.tex", tex)
    print("dw_opt (GHz). dw_lin at 0.01 mW, 2.5 Gb/s. negative = laser red of the cold resonance")
    print("    Q    HW   dw_lin  P  | " + "  ".join(f"{r:>5g} Gb/s" for r in Rs))
    for iq, q in enumerate(Qs):
        for ip, p in enumerate(Ps):
            lead = f"{q:6d} {R['hw_GHz'][iq]:5.2f} {R['dw_lin'][iq]:+6.2f}" if ip == 0 else " " * 19
            print(f"{lead}  {p:g}  | " + "  ".join(f"{R['dw_opt'][iq, ip, ir]:+10.2f}" for ir in range(len(Rs))))


def table_bistable(root, out):
    H = D.loadmat(D.path(root, 4, "hysteresis.mat"))
    HYS, dwc = np.asarray(H["HYS"], float), np.asarray(H["dw_cw"], float)
    HQ = [int(q) for q in np.atleast_1d(H["HQ"])]
    Vc = [float(v) for v in np.atleast_1d(H["V_cw"])]
    tg = {15000: (float(H["d15_25"]), float(H["d15_15"])), 9000: (float(H["d9_25"]), float(H["d9_15"]))}
    rows = []
    for k, q in enumerate((15000, 9000)):
        iq = HQ.index(q)
        for j, v in enumerate(Vc):
            b, r = HYS[:, 0, iq, j, 0] / 4e-3, HYS[:, 0, iq, j, 1] / 4e-3
            bi = np.abs(b - r) > 0.01
            lo, hi = dwc[bi].min(), dwc[bi].max()
            far = min(abs(tg[q][0] - hi), abs(tg[q][1] - hi))
            same = max(abs(b[np.argmin(abs(dwc - x))] - r[np.argmin(abs(dwc - x))]) for x in (0.0,) + tg[q])
            rows.append(f"{band(k)}{q if j == 0 else ''} & {num(v, '.1f')} & {num(lo, '.2f')} to {num(hi, '.2f')} & "
                        f"{hi - lo:.2f} & {np.abs(b - r).max():.3f} & {HYS[:, 3, iq, j, 0].max():.3f} & {far:.2f} \\\\")
            print(f"Q {q} V {v}: largest direction difference at dw = 0 and the two targets {same:.1e} P_in")
    write(out, "table_bistable.tex", "\\begin{tabular}{ccccccc}\n\\toprule\n"
          "$Q$ & $V_{pn}$ (V) & Bistable range (GHz) & Width (GHz) & Largest difference ($P_{\\mathrm{in}}$) & "
          "Peak $\\Delta T$ (K) & Gap to nearest target (GHz) \\\\\n\\midrule\n" + "\n".join(rows)
          + "\n\\bottomrule\n\\end{tabular}\n")
    if D.have(root, 4, "selfpulse.mat"):
        SP = list(np.atleast_1d(D.loadmat(D.path(root, 4, "selfpulse.mat"))["SP"]))
        print(f"self-pulsation: largest Pt pp {max(float(s.Pt_pp) for s in SP) / 4e-3:.1e} P_in, "
              f"largest DT pp {1e3 * max(float(s.DT_pp) for s in SP):.1e} mK over {len(SP)} cw runs")


def table_checks(root, out):
    m = D.loadmat(D.path(root, 5, "results.mat"))
    Qv = [int(q) for q in m["Qv"]]
    rows = [f"{band(i)}{q} & {m['hw_ex'][i]:.4f} & {m['hw_sim'][i]:.4f} & {num(1e3 * (m['hw_sim'][i] - m['hw_ex'][i]), '.2f')} & "
            f"{m['dmax'][i]:.1e} \\\\" for i, q in enumerate(Qv)]
    write(out, "table_linear_limit.tex", "\\begin{tabular}{ccccc}\n\\toprule\n$Q$ & Exact half-width (GHz) & Simulated (GHz) & "
          "Difference (MHz) & Largest $|T_{\\mathrm{sim}}-T_{\\mathrm{exact}}|/P_{\\mathrm{in}}$ \\\\\n\\midrule\n"
          + "\n".join(rows) + "\n\\bottomrule\n\\end{tabular}\n")
    names = {"dt 0.5 ps": "Time step 0.5~ps", "dt 0.25 ps": "Time step 0.25~ps",
             "warm-up doubled": "Warm-up doubled (500~ns ramp, 500~ns hold)", "2x symbols": "Data window 2000 symbols"}
    vn = [str(x) for x in np.atleast_1d(m["var_name"])]
    wst = np.asarray(m["worst"], float)
    assert np.all(np.asarray(m["base_chk"]) == 0)
    pc = lambda x: "$<0.001$" if x < 0.001 else f"{x:.3f}"      # noqa: E731
    rows = [f"{band(i)}{names[n]} & {pc(wst[k, 0])} & {pc(wst[k, 1])} \\\\" for i, (k, n) in
            enumerate((k, n) for k, n in enumerate(vn) if n in names)]
    write(out, "table_convergence.tex", "\\begin{tabular}{lcc}\n\\toprule\nVariant & Largest OMA change (\\%) & Largest VEO change (\\%) \\\\\n"
          "\\midrule\n" + "\n".join(rows) + "\n\\bottomrule\n\\end{tabular}\n")


def table_population(root, out):
    QS = D.pop_qs(root)
    R1 = D.sweeps(root)
    POP = {q: D.pop(root, q) for q in QS}

    # mean oma, veo, log10(ber) and invalid rings of every configuration
    rows, n_all, n_win = [], 0, 0
    for k, q in enumerate(QS):
        d = POP[q]
        first = True
        for ci, (P, r) in enumerate(d["configs"]):
            res = {c: D.case_arrays(d, ci, c) for c in ("AF", "ZDC", "PAC")}
            mo = {c: res[c][0].mean() for c in res}
            mv = {c: res[c][1].mean() for c in res}
            mb = {c: D.lber(res[c][2]).mean() for c in res}
            inv = {c: int(np.sum(res[c][1] < 0)) for c in res}
            win = mo["PAC"] > mo["ZDC"] and mb["PAC"] < mb["ZDC"] and inv["PAC"] <= inv["ZDC"]
            n_all += 1
            n_win += win
            dw = D.target(R1, q, P, r)
            mark = "" if win else "$^{*}$"
            rows.append(f"{band(k)}{q if first else ''} & {P:g} & {r:g} & {num(dw, '.2f')} & {num(mo['AF'], '.3f')} & "
                        f"{mo['ZDC']:.3f} & {mo['PAC']:.3f} & {mv['ZDC']:.3f} & {mv['PAC']:.3f} & "
                        f"{num(mb['AF'], '.2f')} & {num(mb['ZDC'], '.2f')} & {num(mb['PAC'], '.2f')}{mark} & "
                        f"{inv['AF']} & {inv['ZDC']} & {inv['PAC']} \\\\")
            first = False
    print(f"PAC higher mean OMA, lower mean BER and no more invalid rings in {n_win} of {n_all} configurations")
    write(out, "table_population.tex", "\\begin{tabular}{ccccccccccccccc}\n\\toprule\n"
          " & & & & \\multicolumn{3}{c}{Mean OMA} & \\multicolumn{2}{c}{Mean VEO} & \\multicolumn{3}{c}{Mean $\\log_{10}(\\mathrm{BER})$}"
          " & \\multicolumn{3}{c}{Invalid rings} \\\\\n"
          "\\cmidrule(lr){5-7}\\cmidrule(lr){8-9}\\cmidrule(lr){10-12}\\cmidrule(l){13-15}\n"
          "$Q$ & $P_{\\mathrm{in}}$ & Rate & $\\Delta\\omega_{\\mathrm{opt}}/2\\pi$ & AF & ZDC & PAC & ZDC & PAC & AF & ZDC & PAC & "
          "AF & ZDC & PAC \\\\\n & (mW) & (Gb/s) & (GHz) & & & & & & & & & {\\footnotesize of 500} & {\\footnotesize of 2000} & "
          "{\\footnotesize of 2000} \\\\\n\\midrule\n" + "\n".join(rows) + "\n\\bottomrule\n\\end{tabular}\n")

    # receiver floor. the ber is recomputed from the saved level means and spreads
    def mlber(lv, rx):
        mu1, mu0, s1, s0 = lv[:, 0], lv[:, 1], lv[:, 2], lv[:, 3]
        qf = (mu1 - mu0) / (np.sqrt(s1 ** 2 + rx ** 2) + np.sqrt(s0 ** 2 + rx ** 2))
        return np.log10(np.maximum(0.5 * erfc(qf / np.sqrt(2)), 1e-30))
    gaps, chk = {rx: [] for rx in (0.02, 0.05, 0.10)}, 0.0
    for q in QS:
        d = POP[q]
        for ci, (P, r) in enumerate(d["configs"]):
            for rx in gaps:
                z, o = mlber(d["LV_cor"][:, :, ci], rx), mlber(d["LV_opt"][:, :, ci], rx)
                if rx == 0.05:
                    chk = max(chk, np.abs(z - D.lber(d["BER_cor"][:, ci])).max(), np.abs(o - D.lber(d["BER_opt"][:, ci])).max())
                gaps[rx].append((z.mean() - o.mean(), q, P, r))
    print(f"receiver floor check: recomputed BER at 5 % equals the saved BER to {chk:.1e} decades")
    rows = []
    for i, rx in enumerate(gaps):
        g = gaps[rx]
        lo, hi = min(g), max(g)
        where = lambda x: f"$Q={x[1]}$, {x[2]:g}~mW, {x[3]:g}~Gb/s"      # noqa: E731
        fl = '.3f' if abs(lo[0]) < 0.005 else '.2f'
        rows.append(f"{band(i)}{100 * rx:.0f} & {sum(1 for x in g if x[0] > 0)} & {num(lo[0], fl)} & {where(lo)} & "
                    f"{num(hi[0], '.2f')} & {where(hi)} \\\\")
    write(out, "table_receiver_floor.tex", "\\begin{tabular}{cccccc}\n\\toprule\n"
          f"$\\sigma_{{\\mathrm{{rx}}}}$ (\\% of $P_{{\\mathrm{{in}}}}$) & PAC lower (of {len(gaps[0.05])}) & "
          "\\multicolumn{2}{c}{Smallest gap (decades)} & \\multicolumn{2}{c}{Largest gap (decades)} \\\\\n"
          "\\cmidrule(lr){3-4}\\cmidrule(l){5-6}\n & & gap & where & gap & where \\\\\n\\midrule\n"
          + "\n".join(rows) + "\n\\bottomrule\n\\end{tabular}\n")

    # tuning power. d_i = as-fabricated laser detuning of ring i, the same draw in every
    # population file. the tuner shifts the resonance by d_i - D, D = 0 (zdc) or dw_opt (pac).
    # bidirectional tuner: |d_i - D| / eta. heater only, one design offset for every ring:
    # pac costs D / eta more than zdc per ring.
    dfab = np.asarray(POP[QS[0]]["df_fab_GHz"], float).ravel()
    for q in QS[1:]:
        assert np.array_equal(dfab, np.asarray(POP[q]["df_fab_GHz"], float).ravel())
    c, lam = 2.99792458e8, 1553.49e-9
    gpn = c / lam ** 2 * 1e-9 / 1e9
    ETA = [("Wang 2020", "0.15 nm/mW", 0.15 * gpn), ("Jayatilleka 2015", "0.25 nm/mW", 0.25 * gpn),
           ("Sun 2016", "3.8 $\\mu$W/GHz", 1 / 3.8e-3)]
    s0 = np.abs(dfab)
    Ds = [float(R1["dw_opt"][iq, ip, ir]) for iq in range(len(R1["Q_values"])) for ip in range(len(R1["powers"]))
          for ir in range(len(R1["datarates"]))]
    dmean = [np.abs(dfab - dd).mean() - s0.mean() for dd in Ds]          # GHz, change of the mean shift
    pct = [100 * x / s0.mean() for x in dmean]
    p95 = [100 * (np.percentile(np.abs(dfab - dd), 95) - np.percentile(s0, 95)) / np.percentile(s0, 95) for dd in Ds]
    print(f"tuning power: {dfab.size} rings, mean |d| {s0.mean():.2f} GHz, 95th pct {np.percentile(s0, 95):.2f} GHz. "
          f"PAC change of the mean power {min(pct):+.2f} to {max(pct):+.2f} %, of the 95th pct {min(p95):+.2f} to {max(p95):+.2f} %. "
          f"targets {min(Ds):+.2f} to {max(Ds):+.2f} GHz")
    rows = []
    for i, (src, eff, e) in enumerate(ETA):
        rows.append(f"{band(i)}{src} & {eff} & {e:.1f} & {s0.mean() / e:.2f} & {np.percentile(s0, 95) / e:.2f} & "
                    f"{num(1e3 * min(dmean) / e, '.2f')} to {num(1e3 * max(dmean) / e, '.2f')} & "
                    f"{num(1e3 * min(Ds) / e, '.0f')} to {num(1e3 * max(Ds) / e, '.0f')} \\\\")
    write(out, "table_tuning_power.tex", "\\begin{tabular}{llccccc}\n\\toprule\n"
          " & & & \\multicolumn{2}{c}{ZDC power (mW)} & \\multicolumn{2}{c}{PAC minus ZDC ($\\mu$W)} \\\\\n"
          "\\cmidrule(lr){4-5}\\cmidrule(l){6-7}\n"
          "Source & Heater efficiency & $\\eta_h$ (GHz/mW) & mean & 95th pct. & bidirectional, mean & heater only \\\\\n"
          "\\midrule\n" + "\n".join(rows) + "\n\\bottomrule\n\\end{tabular}\n")


TABLES = [(table_dwopt, [(1, "results.mat")]),
          (table_bistable, [(4, "hysteresis.mat")]),
          (table_checks, [(5, "results.mat")]),
          (table_population, [(1, "results.mat"), (6, "pop_results_Q3000.mat")])]


def main():
    if len(sys.argv) < 2:
        sys.exit("usage: python3 tables.py <out folder of matlab> [table folder]")
    root = os.path.abspath(sys.argv[1])
    out = os.path.abspath(sys.argv[2]) if len(sys.argv) > 2 else os.path.join(root, "tables")
    os.makedirs(out, exist_ok=True)
    for f, needs in TABLES:
        if all(D.have(root, s, n) for s, n in needs):
            f(root, out)
        else:
            print(f"{f.__name__}: skipped, input not found yet")


if __name__ == "__main__":
    main()
