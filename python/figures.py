"""thesis figures of the microring part, from the matlab outputs.

  A1_mechanism           oma and veo of the nominal ring versus detuning at Q 3000, 9000,
                         15000, and the thermal and carrier shifts versus power (steps 1, 2)
  A2_population_metrics  mean oma, mean log10(ber) and invalid rings versus Q (step 6)
  A3_veo                 mean veo versus Q and veo histograms of four configurations (step 6)
  A4_eye                 eye diagrams of the step 3 ring, zdc and pac (step 3)
  B2_traces              through-port power and laser position from the hot resonance,
                         nominal ring, Q 15000, 2.5 Gb/s (step 4)
  bistability            cw transmission at 4 mW in both sweep directions, Q 15000 (steps 1, 4, 5)
  B4_robustness          mean veo and mean log10(ber) versus the tuning error (steps 6, 7)

each figure is made only when its inputs exist. pdf and png in the output folder.
run:  python3 figures.py <out folder of matlab> [figure folder, default <out>/figs]
"""
import os
import sys
import functools
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import mrr_data as D

# style: times-like serif, stix math, inward ticks on all sides, framed legends, one
# colour and marker per input power in every figure
PCOL = {2: "#1f4e9c", 3: "#c0392b", 4: "#000000"}
PMARK = {2: "o", 3: "s", 4: "^"}
INK = "#000000"
GRAY = "#7f7f7f"
FULL_W = 6.0            # thesis text width, inches
COL_W = 4.2             # single-panel figures, inches


def apply_style():
    mpl.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Liberation Serif", "STIXGeneral", "DejaVu Serif"],
        "mathtext.fontset": "stix",
        "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8,
        "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "legend.fontsize": 7,
        "axes.linewidth": 0.6, "axes.edgecolor": INK, "axes.labelcolor": INK,
        "axes.spines.top": True, "axes.spines.right": True,
        "axes.grid": False, "axes.axisbelow": True,
        "axes.titlepad": 3.0, "axes.labelpad": 2.0,
        "xtick.direction": "in", "ytick.direction": "in",
        "xtick.top": True, "ytick.right": True,
        "xtick.major.size": 3.0, "ytick.major.size": 3.0,
        "xtick.minor.size": 1.6, "ytick.minor.size": 1.6,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "xtick.minor.width": 0.5, "ytick.minor.width": 0.5,
        "xtick.color": INK, "ytick.color": INK,
        "xtick.major.pad": 2.5, "ytick.major.pad": 2.5,
        "lines.linewidth": 1.0, "lines.markersize": 3.6, "lines.markeredgewidth": 0.8,
        "legend.frameon": True, "legend.framealpha": 1.0, "legend.edgecolor": INK,
        "legend.fancybox": False, "legend.borderpad": 0.35, "legend.borderaxespad": 0.4,
        "legend.handlelength": 2.2, "legend.handletextpad": 0.5, "legend.labelspacing": 0.25,
        "legend.columnspacing": 1.0,
        "patch.linewidth": 0.6,
        "pdf.fonttype": 42, "ps.fonttype": 42,
        "savefig.dpi": 600,
    })


def styled(f):
    """run a figure function with the style applied, then restore the defaults."""
    @functools.wraps(f)
    def wrapped(*a, **k):
        with plt.rc_context():
            apply_style()
            return f(*a, **k)
    return wrapped


def letter(ax, s, x=0.03, y=0.95, ha="left", va="top", **kw):
    """panel label inside the axes, upper left by default."""
    return ax.text(x, y, f"({s})", transform=ax.transAxes, ha=ha, va=va, color=INK, **kw)


def legend_frame(leg):
    fr = leg.get_frame()
    fr.set_linewidth(0.5)
    fr.set_edgecolor(INK)
    return leg


def save(fig, out, name, png_dpi=600, pdf_dpi=600):
    for ext, dpi in (("pdf", pdf_dpi), ("png", png_dpi)):
        fig.savefig(os.path.join(out, f"{name}.{ext}"), dpi=dpi)
    plt.close(fig)
    print("saved", os.path.join(out, name + ".{pdf,png}"))


@styled
def fig_mechanism(root, out):
    R = D.sweeps(root)
    T = D.loadmat(D.path(root, 2, "results.mat"))
    TR = np.atleast_1d(T["TRs"])
    jobs = np.atleast_2d(T["jobs"])            # Q, P (mW), dw, case (1 = zero, 2 = target)
    Pm = T["P"]
    f_zero = Pm.dW0_dVpn * abs(Pm.V_high) / (2 * np.pi * 1e9)     # zero-bit resonance, GHz red of the cold one
    f_one = Pm.dW0_dVpn * abs(Pm.V_low) / (2 * np.pi * 1e9)
    r25 = R["datarates"].index(2.5)
    df = R["df_grid"]
    plt.rcParams["hatch.linewidth"] = 0.4
    pcol = {2.0: PCOL[2], 3.0: PCOL[3], 4.0: PCOL[4], 0.01: GRAY}
    pmrk = {2.0: "o", 3.0: "s", 4.0: "^", 0.01: "D"}
    BAND = "#ececec"
    C_TH, C_FC, C_SUM = "#c0392b", "#1f4e9c", "#000000"      # thermal (red shift), carrier (blue shift), sum

    def curves(q):
        iq = R["Q_values"].index(q)
        out_ = [(0.01, R["oma_lin"][:, iq], R["veo_lin"][:, iq])]
        for ip, p in enumerate(R["powers"]):
            out_.append((p, R["oma_sw"][:, iq, ip, r25], R["veo_sw"][:, iq, ip, r25]))
        return out_

    def panel(ax, q, xlim, ylim, label):
        ax.axvspan(-f_one, -f_zero, color=BAND, zorder=0, lw=0)
        ax.axvline(-(f_one + f_zero) / 2, color=GRAY, lw=0.6, ls=(0, (1, 1.5)), zorder=0.5)
        ax.axvline(0, color=INK, lw=0.5, zorder=1)
        ax.axhline(0, color=INK, lw=0.5, zorder=1)
        for p, o, v in curves(q):
            lw = 0.9 if p < 0.1 else 1.0
            ax.plot(df, o, "-", color=pcol[p], lw=lw, zorder=2 if p < 0.1 else 3)
            ax.plot(df, v, "--", color=pcol[p], lw=lw - 0.1, dashes=(3.5, 2), zorder=2 if p < 0.1 else 3)
            k = int(np.argmax(o))
            ax.plot(df[k], o[k], pmrk[p], color=pcol[p], ms=3.4 if p < 0.1 else 4.0, mfc="white", mew=0.9, zorder=6)
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
        ax.set_xlabel(r"$\Delta\omega/2\pi$ (GHz)")
        ax.set_ylabel("Signed OMA, VEO")
        letter(ax, label)
        ax.text(0.97, 0.95, f"$Q$ = {q}", transform=ax.transAxes, ha="right", va="top")

    def yrange(qs, xlim):
        m = (df >= xlim[0]) & (df <= xlim[1])
        v = np.concatenate([np.concatenate([o[m], w[m]]) for q in qs for _, o, w in curves(q)])
        lo, hi = v.min(), v.max()
        pad = 0.08 * (hi - lo)
        return lo - pad, hi + 0.12 * (hi - lo)

    # shared range for Q 9000 and 15000: from the reddest 4 mW target to the bluest low-power one
    hi_q = [q for q in (9000, 15000) if q in R["Q_values"]]
    tg = [D.target(R, q, p, 2.5) for q in hi_q for p in R["powers"]] + [R["dw_lin"][R["Q_values"].index(q)] for q in hi_q]
    hw = {q: R["hw_GHz"][R["Q_values"].index(q)] for q in R["Q_values"]}
    xl = (5 * np.floor((min(tg) - 1.6 * max(hw[q] for q in hi_q)) / 5),
          5 * np.ceil((max(tg) + 2.2 * max(hw[q] for q in hi_q)) / 5))
    fig = plt.figure(figsize=(FULL_W, 4.2))
    gs = fig.add_gridspec(2, 2, hspace=0.32, wspace=0.28)
    ax3, ax9, ax15, axd = (fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, 0]),
                           fig.add_subplot(gs[1, 1]))
    x3 = (-40, 70)
    panel(ax3, 3000, x3, yrange([3000], x3), "a")
    y_hi = yrange(hi_q, xl)
    panel(ax9, 9000, xl, y_hi, "b")
    panel(ax15, 15000, xl, y_hi, "c")

    # (d) shifts versus input power at Q 9000
    q = 9000
    Ps = sorted(set(jobs[jobs[:, 0] == q, 1]))
    allv = []
    for case, ls, mfc in ((2, "-", None), (1, "--", "white")):
        th, fc, lo, hi = [], [], [], []
        for p in Ps:
            j = int(np.where((jobs[:, 0] == q) & (np.abs(jobs[:, 1] - p) < 1e-9) & (jobs[:, 3] == case))[0][0])
            S = TR[j]
            th.append(S.thermal_GHz); fc.append(S.carrier_GHz)
            lo.append(S.carrier_bitmean_min); hi.append(S.carrier_bitmean_max)
        th, fc = np.array(th), np.array(fc)
        allv += list(th) + list(fc) + lo + hi + list(th + fc)
        if case == 2:
            axd.fill_between(Ps, lo, hi, color=C_FC, alpha=0.15, lw=0, zorder=1)
        else:
            axd.fill_between(Ps, lo, hi, facecolor="none", edgecolor=C_FC, alpha=0.35, hatch="////", lw=0, zorder=1)
        for y, col, mk in ((th, C_TH, "o"), (fc, C_FC, "s")):
            axd.plot(Ps, y, ls, color=col, lw=1.0, marker=mk, ms=3.6, mfc=col if mfc is None else mfc, mec=col,
                     mew=0.8, zorder=3, clip_on=False)
        axd.plot(Ps, th + fc, ls, color=C_SUM, lw=0.6, zorder=4, clip_on=False)
    axd.axhline(0, color=INK, lw=0.5, zorder=1)
    span = max(allv) - min(allv)
    axd.set_ylim(min(min(allv), 0) - 0.12 * span, max(allv) + 0.55 * span)
    axd.set_xlim(-0.15, max(Ps) + 0.15)
    axd.set_xticks(range(0, int(max(Ps)) + 1))
    axd.set_xlabel(r"$P_{\mathrm{in}}$ (mW)")
    axd.set_ylabel("Resonance shift (GHz)")
    axd.tick_params(right=False)
    sec = axd.secondary_yaxis("right", functions=(lambda x: x / hw[q], lambda x: x * hw[q]))
    sec.set_ylabel("Resonance shift (half-widths)")
    sec.tick_params(direction="in", length=3.0, width=0.6)
    letter(axd, "d")
    axd.text(0.97, 0.95, f"$Q$ = {q}, 2.5 Gb/s", transform=axd.transAxes, ha="right", va="top")
    hb = [Line2D([], [], color=C_TH, marker="o", ms=3.4, lw=1.0, label="thermal"),
          Line2D([], [], color=C_FC, marker="s", ms=3.4, lw=1.0, label="carrier"),
          Line2D([], [], color=C_SUM, lw=0.6, label="sum"),
          Line2D([], [], color=INK, ls="-", lw=1.0, label=r"at $\Delta\omega_{\mathrm{opt}}$"),
          Line2D([], [], color=INK, ls="--", lw=1.0, label=r"at $\Delta\omega=0$"),
          Patch(facecolor=C_FC, alpha=0.25, label="carrier, bit to bit")]
    legend_frame(axd.legend(handles=hb, loc="upper left", bbox_to_anchor=(0.0, 0.88), ncol=2))

    hl = [Line2D([], [], color=pcol[p], lw=1.0, marker=pmrk[p], ms=3.6, mfc="white", mew=0.9,
                 label=("0.01 mW" if p < 0.1 else f"{p:g} mW")) for p in (0.01, 2.0, 3.0, 4.0)]
    hl += [Line2D([], [], color=INK, ls="-", lw=1.0, label="OMA"),
           Line2D([], [], color=INK, ls="--", dashes=(3.5, 2), lw=0.9, label="VEO"),
           Patch(facecolor=BAND, edgecolor="none", label="drive swing"),
           Line2D([], [], color=GRAY, lw=0.6, ls=(0, (1, 1.5)), label="bias, $-$2.5 V")]
    legend_frame(ax3.legend(handles=hl, loc="lower right", ncol=2, handlelength=1.8))
    fig.subplots_adjust(left=0.075, right=0.93, top=0.98, bottom=0.085)
    save(fig, out, "A1_mechanism")
    print("  x range of (b), (c):", xl)


def _pop_axes_data(root):
    S = D.stats(root)
    QS = list(S)
    PS = sorted({p for q in QS for (p, r) in S[q]})
    RS = sorted({r for q in QS for (p, r) in S[q]})
    return S, QS, PS, RS


def _power_handles(PS, PMS):
    handles = [Line2D([], [], color=PCOL[int(p)], lw=1.0, marker=PMARK[int(p)], ms=PMS[int(p)], mfc=PCOL[int(p)],
                      mec=PCOL[int(p)], label=f"{p:g} mW") for p in PS]
    handles += [Line2D([], [], color=INK, lw=1.0, ls="--", dashes=(4, 2), marker="o", ms=3.8, mfc="white", mec=INK,
                       mew=0.9, label="ZDC"),
                Line2D([], [], color=INK, lw=1.0, ls="-", marker="o", ms=3.8, mfc=INK, mec=INK, label="PAC")]
    return handles


@styled
def fig_population(root, out):
    # zdc dashed with open markers, pac solid with filled markers. af is left out.
    S, QS, PS, RS = _pop_axes_data(root)
    PMS = {2: 3.8, 3: 3.4, 4: 4.2}

    def val(q, p, r, case, key):
        return S[q][(p, r)][case][key]

    ROWS = [("Mean OMA", "oma"), ("Mean log$_{10}$(BER)", "lber")]
    if any(S[q][k][c]["inv"] for q in S for k in S[q] for c in ("ZDC", "PAC")):     # only when there is something to show
        ROWS.append(("Invalid rings (%)", "inv_pct"))

    def lims(key):
        v = [val(q, p, r, c, key) for q in QS for p in PS for r in RS for c in ("ZDC", "PAC")]
        lo, hi = min(v), max(v)
        if key == "inv_pct":
            return (-4, 104) if hi > 20 else (-0.05 * max(hi, 1) - 0.2, max(hi, 1) * 1.15 + 0.2)
        pad = 0.08 * (hi - lo if hi > lo else 1)
        lo = min(lo - pad, 0) if key == "oma" else lo - pad
        hi = hi + pad if key == "oma" else min(hi + pad, 0.5)
        return lo, hi

    fig, axes = plt.subplots(len(ROWS), len(RS), figsize=(FULL_W, 1.55 * len(ROWS) + 0.2), sharex=True,
                             sharey="row", squeeze=False)
    for ri, (ylab, key) in enumerate(ROWS):
        yl = lims(key)
        for ci, r in enumerate(RS):
            ax = axes[ri, ci]
            if key == "oma":
                ax.axhline(0, color=INK, lw=0.5, zorder=1)
            for p in PS:
                pi = int(round(p))
                ax.plot(QS, [val(q, p, r, "PAC", key) for q in QS], "-", color=PCOL[pi], lw=1.0, marker=PMARK[pi],
                        ms=PMS[pi], mfc=PCOL[pi], mec=PCOL[pi], zorder=3)
                ax.plot(QS, [val(q, p, r, "ZDC", key) for q in QS], "--", color=PCOL[pi], lw=1.0, dashes=(4, 2),
                        marker=PMARK[pi], ms=PMS[pi], mfc="white", mec=PCOL[pi], mew=0.9, zorder=4)
            ax.set_ylim(*yl)
            ax.set_xlim(min(QS) - 700, max(QS) + 700)
            ax.set_xticks(QS)
            if key == "lber":
                letter(ax, "abcdefghi"[len(RS) * ri + ci], y=0.05, va="bottom")
            else:
                letter(ax, "abcdefghi"[len(RS) * ri + ci])
            if ri == 0:
                ax.set_title(f"{r:g} Gb/s")
            if ci == 0:
                ax.set_ylabel(ylab)
            if ri == len(ROWS) - 1:
                ax.set_xlabel("$Q$")
    legend_frame(fig.legend(handles=_power_handles(PS, PMS), loc="upper center", ncol=5, bbox_to_anchor=(0.5, 0.995)))
    H = fig.get_size_inches()[1]
    fig.subplots_adjust(left=0.085, right=0.985, top=1 - 0.46 / H, bottom=0.36 / H, hspace=0.10, wspace=0.08)
    save(fig, out, "A2_population_metrics")


@styled
def fig_veo(root, out):
    # (a), (b) mean veo versus Q. (c)-(f) veo histograms, pac gray fill, zdc dashed outline.
    HIST_CASES = [(9000, 2, 2.5), (9000, 4, 2.5), (9000, 4, 15.0), (15000, 4, 2.5)]
    S, QS, PS, RS = _pop_axes_data(root)
    SHADE = "#ececec"
    PMS = {2: 3.8, 3: 3.4, 4: 4.2}
    FILL, EDGE = "#9a9a9a", "#4d4d4d"

    fig = plt.figure(figsize=(FULL_W, 4.0))
    gs = fig.add_gridspec(2, 4, height_ratios=[1.0, 0.9], hspace=0.42, wspace=0.10)
    top = [fig.add_subplot(gs[0, 0:2])]
    top.append(fig.add_subplot(gs[0, 2:4], sharey=top[0]))
    hist_axes = [fig.add_subplot(gs[1, k]) for k in range(4)]
    for k in range(1, 4):
        hist_axes[k].sharex(hist_axes[0])
        hist_axes[k].sharey(hist_axes[0])

    allm = [S[q][(p, r)][c]["veo"] for q in QS for p in PS for r in RS for c in ("ZDC", "PAC")]
    ylo, yhi = min(allm) - 0.06, max(allm) + 0.08
    if min(allm) < 0.1:                 # show the invalid region only when the data come near it
        ylo = min(ylo, -0.05)
    for ax, r, label in zip(top, RS[:2], "ab"):
        ax.axhspan(-1, 0, color=SHADE, zorder=0, lw=0)
        ax.axhline(0, color=INK, lw=0.5, zorder=1)
        for p in PS:
            pi = int(round(p))
            ax.plot(QS, [S[q][(p, r)]["PAC"]["veo"] for q in QS], "-", color=PCOL[pi], lw=1.0, marker=PMARK[pi],
                    ms=PMS[pi], mfc=PCOL[pi], mec=PCOL[pi], zorder=3)
            ax.plot(QS, [S[q][(p, r)]["ZDC"]["veo"] for q in QS], "--", color=PCOL[pi], lw=1.0, dashes=(4, 2),
                    marker=PMARK[pi], ms=PMS[pi], mfc="white", mec=PCOL[pi], mew=0.9, zorder=4)
        ax.set_xlim(min(QS) - 700, max(QS) + 700)
        ax.set_xticks(QS)
        ax.set_ylim(ylo, yhi)
        ax.set_xlabel("$Q$")
        ax.set_title(f"{r:g} Gb/s")
        letter(ax, label)
    top[0].set_ylabel("Mean VEO")
    plt.setp(top[1].get_yticklabels(), visible=False)
    handles = _power_handles(PS, PMS)
    if ylo < 0:
        handles.append(Patch(facecolor=SHADE, edgecolor="none", label="invalid, VEO $<$ 0"))
    legend_frame(fig.legend(handles=handles, loc="upper center", ncol=len(handles), bbox_to_anchor=(0.5, 0.995)))

    data = {}
    for q, p, r in HIST_CASES:
        d = D.pop(root, q)
        ci = D.cfg_index(d, p, r)
        data[(q, p, r)] = (d["EO_cor"][:, ci], d["EO_opt"][:, ci])
    allv = np.concatenate([np.concatenate(v) for v in data.values()])
    bw = 0.01
    lo = np.floor(((min(allv.min(), 0) if allv.min() < 0.05 else allv.min()) - 0.02) / bw) * bw
    hi = np.ceil((allv.max() + 0.02) / bw) * bw
    edges = np.arange(lo, hi + bw / 2, bw)
    ymax = 0
    for ax, key, label in zip(hist_axes, HIST_CASES, "cdef"):
        zdc, pac = data[key]
        n_r = zdc.size
        ax.axvspan(lo - 1, 0, color=SHADE, zorder=0, lw=0)
        ax.axvline(0, color=INK, lw=0.5, zorder=1)
        for vals, filled in ((pac, True), (zdc, False)):
            n, _ = np.histogram(vals, bins=edges)
            pct = 100.0 * n / n_r
            ymax = max(ymax, pct.max())
            nz = np.nonzero(n)[0]
            a, b = nz[0], nz[-1] + 1
            if filled:
                ax.stairs(pct[a:b], edges[a:b + 1], fill=True, color=FILL, alpha=0.6, linewidth=0, zorder=2)
                ax.stairs(pct[a:b], edges[a:b + 1], color=EDGE, linewidth=0.6, baseline=None, zorder=2)
            else:
                ax.stairs(pct[a:b], edges[a:b + 1], color=INK, linewidth=0.9, linestyle=(0, (3, 1.5)), baseline=None,
                          zorder=3)
        q, p, r = key
        ax.set_title(f"$Q$ = {q}, {p:g} mW, {r:g} Gb/s", fontsize=7.5)
        letter(ax, label, x=0.97, ha="right")
        ax.set_xlabel("VEO")
    hist_axes[0].set_ylabel("Rings per bin (%)")
    hist_axes[0].set_xlim(lo, hi)
    hist_axes[0].set_ylim(0, ymax * 1.12)
    for ax in hist_axes[1:]:
        plt.setp(ax.get_yticklabels(), visible=False)
    legend_frame(hist_axes[0].legend(
        handles=[Patch(facecolor=FILL, alpha=0.6, edgecolor=EDGE, linewidth=0.6, label="PAC"),
                 Patch(facecolor="none", edgecolor=INK, linestyle=(0, (3, 1.5)), linewidth=0.9, label="ZDC")],
        loc="upper left", handlelength=1.4))
    fig.subplots_adjust(left=0.075, right=0.985, top=0.885, bottom=0.095)
    save(fig, out, "A3_veo")
    for key, (z, a) in data.items():
        print(f"  {key}: ZDC VEO {z.min():+.3f} to {z.max():+.3f} ({(z < 0).sum()} below 0), "
              f"PAC {a.min():+.3f} to {a.max():+.3f} ({(a < 0).sum()} below 0)")


@styled
def fig_eyes(root, out):
    # each panel overlays the 250 four-symbol groups of the data window. oma and veo are
    # recomputed from the waveforms first and must equal the matlab values.
    E = D.loadmat(D.path(root, 3, "results.mat"))
    bits = np.asarray(E["prbs_bits"]).astype(int)
    N_SYM, GROUP = int(E["N_sym"]), 4
    names = [str(x) for x in np.atleast_1d(E["case_name"])]

    def case(name, p, rate):
        for c in range(len(names)):
            if names[c] == name and E["case_P_mW"][c] == p and abs(E["case_rate"][c] - rate) < 1e-9:
                return c
        raise KeyError((name, p, rate))

    def center(c):
        pt = np.asarray(E["Pt_all"][c], float).ravel()
        m = int(E["Mbit_all"][c])
        ctr = int(np.floor(m / 2 + 0.5))                 # matlab round(m/2)
        return pt, m, ctr, pt[np.arange(N_SYM) * m + ctr - 1]

    for c in range(len(names)):
        pt, m, ctr, cs = center(c)
        pin = E["case_P_mW"][c] * 1e-3
        oma = (cs[bits == 1].mean() - cs[bits == 0].mean()) / pin
        veo = (cs[bits == 1].min() - cs[bits == 0].max()) / pin
        assert abs(oma - E["oma"][c]) < 1e-12 and abs(veo - E["veo"][c]) < 1e-12, c

    TRACE = "#1f4e9c"
    cols = [(2, 2.5), (4, 2.5), (4, 15.0)]
    rows = ["ZDC", "PAC"]
    ytop = 1.0
    for nm in rows:
        for P, rate in cols:
            pt, m, ctr, cs = center(case(nm.lower(), P, rate))
            ytop = max(ytop, pt.max() / (P * 1e-3) * 1.02)
    fig, axes = plt.subplots(2, 3, figsize=(FULL_W, 3.55), sharex=True, sharey=True)
    for ri, nm in enumerate(rows):
        for ci, (P, rate) in enumerate(cols):
            ax = axes[ri, ci]
            c = case(nm.lower(), P, rate)
            pt, m, ctr, cs = center(c)
            pin = P * 1e-3
            t = np.arange(GROUP * m) / m
            traces = (pt / pin).reshape(N_SYM // GROUP, GROUP * m)
            ax.add_collection(LineCollection([np.column_stack([t, tr]) for tr in traces], colors=TRACE,
                                             linewidths=0.3, alpha=0.10, rasterized=True))
            x0 = (ctr - 1) / m
            ax.axvline(x0, color=GRAY, lw=0.6, ls=(0, (3, 2)), zorder=3)
            top_, bot = cs[bits == 1].min() / pin, cs[bits == 0].max() / pin
            veo = float(E["veo"][c])
            if veo > 0:
                ax.plot([x0, x0], [bot, top_], color=INK, lw=1.4, zorder=5, solid_capstyle="butt")
                for y in (bot, top_):
                    ax.plot([x0 - 0.06, x0 + 0.06], [y, y], color=INK, lw=0.9, zorder=5)
            ax.text(0.97, 0.95, f"VEO {veo:.3f}".replace("-", "−"), transform=ax.transAxes, ha="right", va="top",
                    fontsize=7.5, color=INK, bbox=dict(facecolor="white", edgecolor="none", pad=1.0))
            letter(ax, "abcdef"[3 * ri + ci], bbox=dict(facecolor="white", edgecolor="none", pad=1.0))
            ax.set_xlim(0, GROUP)
            ax.set_ylim(0, ytop)
            ax.set_xticks(range(GROUP + 1))
            if ri == 0:
                ax.set_title(f"{P} mW, {rate:g} Gb/s")
            if ci == 0:
                ax.set_ylabel(r"$P_{\mathrm{thr}}/P_{\mathrm{in}}$")
                ax.text(-0.30, 0.5, nm, transform=ax.transAxes, rotation=90, ha="center", va="center", fontsize=8.5,
                        fontweight="bold")
            if ri == 1:
                ax.set_xlabel("Time (UI)")
    handles = [Line2D([], [], color=TRACE, lw=1.0, label="through-port power, 250 four-symbol groups"),
               Line2D([], [], color=GRAY, lw=0.8, ls=(0, (3, 2)), label="decision instant"),
               Line2D([], [], color=INK, lw=1.4, marker="_", ms=6, mew=0.9, label="vertical eye opening (VEO)")]
    legend_frame(fig.legend(handles=handles, loc="upper center", ncol=3, bbox_to_anchor=(0.5, 0.995),
                            handlelength=2.4, columnspacing=1.6))
    fig.subplots_adjust(left=0.095, right=0.99, top=0.84, bottom=0.12, hspace=0.12, wspace=0.08)
    save(fig, out, "A4_eye")
    for c in range(len(names)):
        print(f"  {names[c]} {E['case_P_mW'][c]:g} mW {E['case_rate'][c]:g} Gb/s dw {E['case_dw'][c]:+.2f} GHz: "
              f"OMA {E['oma'][c]:+.3f} VEO {E['veo'][c]:+.3f} log10 BER {np.log10(max(E['ber'][c], 1e-30)):.2f}")


@styled
def fig_traces(root, out):
    # columns: zdc 2 mW (case1), zdc 4 mW (case2), pac 4 mW (case3), data bits 20 to 39.
    # rows: through-port power with the decision samples, and the laser position from the
    # hot resonance in half-widths (red side shaded). every case must match step 1 or 3.
    TRm = D.loadmat(D.path(root, 4, "traces.mat"))
    TR = {t.name: t for t in np.atleast_1d(TRm["TRB"])}
    bits = np.asarray(TRm["prbs_bits"]).astype(int)
    zoom = np.asarray(TRm["zoom_bits"]).astype(int)          # 1-based data bit numbers
    for t in TR.values():
        assert float(t.check) == 0.0, (t.name, t.check)
    RED_BAND = "#f3dede"
    TRACE = PCOL[2]

    def dec_idx(t, n):
        m = int(t.M_bit)
        ctr = max(int(np.floor(m / 2 + 0.5)), 1)        # matlab round(M/2)
        return np.arange(n) * m + ctr - 1

    cases = [("case1", "ZDC, 2 mW"), ("case2", "ZDC, 4 mW"), ("case3", None)]
    zb = bits[zoom[0] - 1: zoom[-1]]
    fig, axes = plt.subplots(2, 3, figsize=(FULL_W, 3.8), sharex=True)
    for ci, (nm, title) in enumerate(cases):
        t = TR[nm]
        Z = np.asarray(t.Z, float)
        m = int(t.M_bit)
        tt = np.arange(Z.shape[0]) / m
        k = dec_idx(t, len(zb))
        if title is None:
            title = "PAC, 4 mW, " + r"$\Delta\omega/2\pi$" + f" = {float(t.dw):.2f} GHz".replace("-", "−")
        ax = axes[0, ci]
        y = Z[:, 0] / float(t.P_in_W)
        ax.plot(tt, y, color=TRACE, lw=0.9)
        ax.plot(tt[k][zb == 1], y[k][zb == 1], "o", ms=2.8, color=INK, zorder=5)
        ax.plot(tt[k][zb == 0], y[k][zb == 0], "o", ms=2.8, mfc="white", mec=INK, mew=0.7, zorder=5)
        ax.set_ylim(0, 1.1)
        ax.set_title(f"{title}\nOMA {float(t.oma):.3f}, VEO {float(t.veo):.3f}")
        letter(ax, "abc"[ci], y=0.97)
        ax = axes[1, ci]
        d = Z[:, 5]
        ax.axhspan(-1.0, 0, color=RED_BAND, lw=0, zorder=0)
        ax.axhline(0, color=INK, lw=0.5)
        ax.plot(tt, d, color=TRACE, lw=0.9)
        ax.plot(tt[k][zb == 1], d[k][zb == 1], "o", ms=2.8, color=INK, zorder=5)
        ax.plot(tt[k][zb == 0], d[k][zb == 0], "o", ms=2.8, mfc="white", mec=INK, mew=0.7, zorder=5)
        ax.set_ylim(-0.6, 3.0)
        ax.set_xlim(0, len(zb))
        ax.set_xticks(range(0, len(zb) + 1, 5))
        ax.set_xlabel(f"Time (UI), data bits {zoom[0]} to {zoom[-1]}")
        letter(ax, "def"[ci], y=0.03, va="bottom", x=0.015)
    axes[0, 0].set_ylabel(r"$P_{\mathrm{thr}}/P_{\mathrm{in}}$")
    axes[1, 0].set_ylabel(r"$\delta_\omega$ (half-widths)")
    for ax in axes[:, 1:].ravel():
        ax.tick_params(labelleft=False)
    handles = [Line2D([], [], ls="", marker="o", ms=3.5, color=INK, label="decision sample, one bit"),
               Line2D([], [], ls="", marker="o", ms=3.5, mfc="white", mec=INK, label="decision sample, zero bit"),
               Patch(fc=RED_BAND, ec="none", label="laser red of the hot resonance")]
    legend_frame(fig.legend(handles=handles, loc="upper center", ncol=5, bbox_to_anchor=(0.5, 1.0), fontsize=7))
    fig.subplots_adjust(left=0.075, right=0.99, top=0.80, bottom=0.12, hspace=0.12, wspace=0.07)
    save(fig, out, "B2_traces", png_dpi=400, pdf_dpi=None)


@styled
def fig_bistability(root, out):
    # nominal ring at Q 15000, zero-bit voltage (-0.5 V), 4 mW, swept blue to red (solid)
    # and red to blue (dashed), with the 1 uW response of the cold ring. the numbers for
    # Q 9000 are printed as well.
    V_ZERO = -0.5                        # zero-bit junction voltage (V)
    QS = (15000, 9000)
    RATE = 2.5
    H = D.loadmat(D.path(root, 4, "hysteresis.mat"))
    p_mW = float(H["P_cw_mW"])
    assert p_mW == 4.0, p_mW
    hq = [int(q) for q in np.atleast_1d(H["HQ"])]
    vc = [float(v) for v in np.atleast_1d(H["V_cw"])]
    hys = np.asarray(H["HYS"], float)
    dw = np.asarray(H["dw_cw"], float)
    curves = {q: (hys[:, 0, hq.index(q), vc.index(V_ZERO), 0] / (p_mW * 1e-3),       # blue to red
                  hys[:, 0, hq.index(q), vc.index(V_ZERO), 1] / (p_mW * 1e-3)) for q in QS}
    dWdV = float(H["dW0_dVpn"])
    R1 = D.sweeps(root)
    tgt = {q: D.target(R1, q, p_mW, RATE) for q in QS}

    # cold ring: the 1 uW response at the bias (-2.5 V) moved to the zero-bit voltage. the
    # linear response depends on dw - (dw0/dv) v only, so the -2.5 V curve moves 5 GHz
    # toward positive detuning.
    R5 = D.loadmat(D.path(root, 5, "results.mat"))
    dfg = np.asarray(R5["dfg"], float)
    qv = [int(q) for q in np.atleast_1d(R5["Qv"])]
    shift = dWdV * (V_ZERO - (-2.5)) / (2 * np.pi * 1e9)
    cold = {q: (dfg + shift, np.asarray(R5["T_sim"], float)[:, qv.index(q)]) for q in QS}
    for q in QS:                                         # the cold zero-bit dip sits at (dw0/dv)*(-0.5 V)/2pi
        x, t = cold[q]
        assert abs(x[np.argmin(t)] - dWdV * V_ZERO / (2 * np.pi * 1e9)) < 0.26, (q, x[np.argmin(t)])

    COL, SHADE = PCOL[4], "#bdbdbd"
    XL = (-40, 20)
    report = []
    for q in QS:
        b2r, r2b = curves[q]
        bi = np.abs(b2r - r2b) > 0.01
        d = tgt[q]
        k = int(np.argmin(np.abs(dw - d)))
        assert abs(dw[k] - d) < 1e-9
        lo, hi = (dw[bi].min(), dw[bi].max()) if bi.any() else (np.nan, np.nan)
        report.append(f"  Q {q}, {p_mW:g} mW, {V_ZERO} V: bistable {lo:+.2f} to {hi:+.2f} GHz (largest difference "
                      f"{np.abs(b2r - r2b).max():.3f} P_in), PAC target {d:+.2f} GHz ({d - hi:+.2f} GHz from the blue "
                      f"edge), blue-to-red {b2r[k]:.3f}, red-to-blue {r2b[k]:.3f} P_in")

    fig, ax = plt.subplots(figsize=(COL_W, 2.95))
    q = 15000
    x, t = cold[q]
    m = (x >= XL[0]) & (x <= XL[1])
    ax.plot(x[m], t[m], color=GRAY, lw=0.9, zorder=2)
    ax.axvline(0, color=INK, lw=0.5, zorder=1)
    b2r, r2b = curves[q]
    bi = np.abs(b2r - r2b) > 0.01
    ax.fill_between(dw, np.minimum(b2r, r2b), np.maximum(b2r, r2b), where=bi, color=SHADE, alpha=0.55,
                    lw=0, zorder=1)
    ax.plot(dw, b2r, "-", color=COL, lw=1.15, zorder=4)
    ax.plot(dw, r2b, "--", color=COL, lw=1.05, dashes=(3.5, 1.8), zorder=4)
    d = tgt[q]
    k = int(np.argmin(np.abs(dw - d)))
    ax.plot([d, d], [0, 1.06], ":", color=COL, lw=0.9, zorder=3)
    assert abs(b2r[k] - r2b[k]) < 0.01                  # both sweep directions reach the same state at the target
    ax.plot(d, b2r[k], "o", ms=4.0, color=COL, mfc=COL, mec=COL, zorder=7)
    ax.set_xlim(*XL)
    ax.set_ylim(0, 1.06)
    ax.tick_params(axis="x", pad=5.0)                # keeps the -40 label clear of the 0.0 label at the corner
    ax.set_xlabel(r"$\Delta\omega/2\pi$ (GHz)")
    ax.set_ylabel(r"$P_{\mathrm{thr}}/P_{\mathrm{in}}$")
    ax.set_title(f"$Q$ = {q}, 4 mW, zero-bit voltage ($-$0.5 V)", fontsize=7.5, pad=3)
    ax.text(1.0, 1.0, "ZDC", ha="left", va="center", fontsize=7, zorder=9)
    handles = [Line2D([], [], color=GRAY, lw=0.9, label="1 $\\mu$W (cold)"),
               Line2D([], [], color=INK, lw=1.1, ls="-", label="blue to red"),
               Line2D([], [], color=INK, lw=1.0, ls="--", dashes=(3.5, 1.8), label="red to blue"),
               Patch(facecolor="#d0d0d0", edgecolor="none", label="bistable"),
               Line2D([], [], color=INK, lw=0.9, ls=":", marker="o", ms=4, mfc=INK, mec=INK,
                      label="PAC target, 2.5 Gb/s")]
    legend_frame(fig.legend(handles=handles, loc="upper center", ncol=3, bbox_to_anchor=(0.53, 1.0), fontsize=6.5,
                            handlelength=1.8, columnspacing=0.9, handletextpad=0.4, labelspacing=0.2))
    fig.subplots_adjust(left=0.15, right=0.97, top=0.845, bottom=0.14)
    save(fig, out, "bistability")
    print("\n".join(report))


@styled
def fig_tuning_error(root, out):
    # mean veo (top) and mean log10(ber) (bottom) over the rings versus sigma_res for the
    # four configurations, zdc dashed with open markers, pac solid with filled markers.
    # sigma_res = 0 must equal the step 6 results exactly.
    M = D.loadmat(D.path(root, 7, "results.mat"))
    VEO, BER = np.asarray(M["VEO"], float), np.asarray(M["BER"], float)
    C, sig = np.atleast_2d(M["C"]), np.asarray(M["sig_all"], float)
    for k, (q, p, r) in enumerate(C):
        d = D.pop(root, int(q))
        ci = D.cfg_index(d, p, r)
        assert np.array_equal(VEO[:, 0, 0, k], d["EO_cor"][:, ci]) and np.array_equal(VEO[:, 0, 1, k], d["EO_opt"][:, ci])
    COL = PCOL[4]
    titles = [f"$Q$ = {int(q)}, {p:g} mW, {r:g} Gb/s" for q, p, r in C]
    fig, axes = plt.subplots(2, 4, figsize=(FULL_W, 3.9), sharex=True)
    for k in range(len(C)):
        for t, (ls, mfc) in enumerate((("--", "white"), ("-", COL))):
            mv = VEO[:, :, t, k].mean(axis=0)
            mb = np.log10(np.maximum(BER[:, :, t, k], 1e-30)).mean(axis=0)
            axes[0, k].plot(sig, mv, ls, color=COL, marker="^", mfc=mfc, mec=COL)
            axes[1, k].plot(sig, mb, ls, color=COL, marker="^", mfc=mfc, mec=COL)
        axes[0, k].set_title(titles[k], fontsize=7.5)
        axes[1, k].set_xlabel(r"$\sigma_{\mathrm{res}}$ (GHz)")
        axes[1, k].set_xticks(sig)
        axes[1, k].set_xlim(-0.3, 2.3)
        for row, lt in ((0, "abcd"[k]), (1, "efgh"[k])):
            ax = axes[row, k]
            ylo, yhi = ax.get_ylim()
            ax.set_ylim(ylo - 0.2 * (yhi - ylo), yhi)      # room for the panel letter under the data
            letter(ax, lt, y=0.04, va="bottom", x=0.05)
    axes[0, 0].set_ylabel("Mean VEO")
    axes[1, 0].set_ylabel(r"Mean $\log_{10}(\mathrm{BER})$")
    handles = [Line2D([], [], color=COL, ls="--", marker="^", mfc="white", mec=COL, label="ZDC"),
               Line2D([], [], color=COL, ls="-", marker="^", mfc=COL, mec=COL, label="PAC")]
    legend_frame(fig.legend(handles=handles, loc="upper center", ncol=2, bbox_to_anchor=(0.5, 1.0)))
    fig.subplots_adjust(left=0.08, right=0.99, top=0.88, bottom=0.11, hspace=0.12, wspace=0.38)
    save(fig, out, "B4_robustness", png_dpi=400, pdf_dpi=None)


FIGS = [(fig_mechanism, [(1, "results.mat"), (2, "results.mat")]),
        (fig_population, [(6, "pop_results_Q3000.mat")]),
        (fig_veo, [(6, "pop_results_Q9000.mat"), (6, "pop_results_Q15000.mat")]),
        (fig_eyes, [(3, "results.mat")]),
        (fig_traces, [(4, "traces.mat")]),
        (fig_bistability, [(1, "results.mat"), (4, "hysteresis.mat"), (5, "results.mat")]),
        (fig_tuning_error, [(7, "results.mat"), (6, "pop_results_Q9000.mat"), (6, "pop_results_Q15000.mat")])]


def main():
    if len(sys.argv) < 2:
        sys.exit("usage: python3 figures.py <out folder of matlab> [figure folder]")
    root = os.path.abspath(sys.argv[1])
    out = os.path.abspath(sys.argv[2]) if len(sys.argv) > 2 else os.path.join(root, "figs")
    os.makedirs(out, exist_ok=True)
    for f, needs in FIGS:
        if all(D.have(root, s, n) for s, n in needs):
            f(root, out)
        else:
            print(f"{f.__name__}: skipped, input not found yet")


if __name__ == "__main__":
    main()
