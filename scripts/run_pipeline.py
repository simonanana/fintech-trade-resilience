#!/usr/bin/env python3
"""Run the full analysis: integrity checks, diagnostics, estimation, inference,
alternative designs, figures and tables.

Usage
-----
    python scripts/run_pipeline.py                         # uses data/processed/panel_eastasia_v4.parquet
    python scripts/run_pipeline.py --panel path/to/panel.parquet --n-ri 5000
    python scripts/run_pipeline.py --skip-ri               # fast run, no randomization inference
    python scripts/run_pipeline.py --panel data/processed/panel_synthetic.parquet --out demo_output

Outputs go to <out>/results/tables and <out>/results/figures.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ftres import (corridors, data, diagnostics, hypotheses, inference, plots,  # noqa: E402
                   report, spec_curve, synth)
from ftres.config import (FE_MAIN, INDEX_VARIANT, INDEX_VARIANTS, N_RI_DEFAULT,  # noqa: E402
                          REPO_ROOT, TREATED_EXPORTER, Paths)
from ftres.estimation import coef, fit  # noqa: E402
from ftres.style import apply_style, savefig  # noqa: E402

#: Values recorded in the paper, used only when --skip-ri is set.
PAPER_RI = {"b_obs": -1.4627, "sd": 1.1237, "p": 0.1833, "n_draws": 300,
            "draws": np.array([]), "mc_se": 0.0223}


def step(msg: str) -> None:
    print(f"\n[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def main(argv=None) -> dict:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--panel", type=Path, default=None, help="processed panel (.parquet)")
    ap.add_argument("--out", type=Path, default=REPO_ROOT, help="output root (default: repo root)")
    ap.add_argument("--n-ri", type=int, default=N_RI_DEFAULT, help="randomization-inference draws")
    ap.add_argument("--skip-ri", action="store_true", help="skip randomization inference")
    ap.add_argument("--no-figures", action="store_true")
    ap.add_argument("--variant", default=INDEX_VARIANT, help="fintech index standardisation")
    ap.add_argument("--partner-findex", type=Path, default=None,
                    help="CSV (iso3, year, findex_digital_pay) to run Route B")
    args = ap.parse_args(argv)

    logging.basicConfig(level=logging.WARNING, format="  [%(levelname)s] %(message)s")
    # pyfixest reports separation and dropped collinear terms as UserWarnings;
    # separation is expected here (all zeros sit in never-trading dyads).
    warnings.filterwarnings("ignore", category=UserWarning, module="pyfixest")
    paths = Paths(root=args.out).ensure()
    panel_path = args.panel or Paths().panel
    tab, figdir = paths.tables, paths.figures
    figs = not args.no_figures
    if figs:
        apply_style()
    R: dict = {}

    # ------------------------------------------------------------------ data
    step(f"load and prepare {panel_path}")
    D = data.prepare_panel(data.load_panel(panel_path), variant=args.variant)
    spl = data.splice_check(D)
    spl.to_csv(tab / "tab0_splice_check.csv")
    bal_tab, balanced, failing = data.balance_check(D)
    bal_tab.to_csv(tab / "tab0_balance_check.csv")
    S = data.build_samples(D, balanced)
    main_ = S.main
    print(data.sample_summary(S).to_string(index=False))
    print(f"  exporter-years failing the coverage rule: {failing}")
    R.update(n_exporters=main_.iso3_o.nunique(), n_obs=len(main_), n_dyads=main_.pair.nunique())

    # ----------------------------------------------------------- diagnostics
    step("identification diagnostics")
    var = diagnostics.variation_by_exporter(S.all10)
    var.to_csv(tab / "tabA1_variation_by_exporter.csv", index=False)
    print(diagnostics.wave1_changes(S.all10).to_string(index=False))
    stepwise = diagnostics.stepwise_deletion(S.all10)
    stepwise.to_csv(tab / "tabA2_core_battery_by_sample.csv", index=False)
    print(stepwise.round(3).to_string(index=False))

    ne = diagnostics.effective_units(main_)
    R["N_eff"], R["neff_shares"] = ne["N_eff"], ne["shares"].round(4).to_dict()
    print(f"  N_eff = {ne['N_eff']:.2f} (nominal {ne['n_units']} exporters, {ne['n_clusters']} dyads)")
    pd.DataFrame({"share": ne["shares"]}).to_csv(tab / "tabA5_neff_shares.csv")

    scale = diagnostics.moderator_scale(main_)
    R["sd_moderator"], R["mean_moderator"] = scale["sd"], scale["mean"]

    b1, _, t1 = coef(fit(f"trade_usd ~ ln_tariff | {FE_MAIN}", main_), "ln_tariff")
    m_int = fit(hypotheses.F_H1, main_)
    b1i, b3 = coef(m_int, "ln_tariff")[0], coef(m_int, "tXfin")[0]
    R.update(beta1_no_interaction=b1, t_beta1=t1, beta1_interacted=b1i, beta3=b3,
             turning_point=diagnostics.turning_point(b1i, b3))
    dec = diagnostics.elasticity_decomposition(b1i, b3, scale["values"], b1)
    dec.to_csv(tab / "tabA6_elasticity_decomposition.csv")
    fin = scale["values"]
    R["index_treated"] = float(fin.get(TREATED_EXPORTER, np.nan))
    R["implied_treated"] = b1i + b3 * R["index_treated"]
    R["implied_mean"] = b1i + b3 * scale["mean"]
    beyond = fin[fin < R["turning_point"]] if b3 < 0 else fin[fin > R["turning_point"]]
    R["units_beyond_turning_point"] = ", ".join(beyond.index)

    diagnostics.splice_level_shift(D).to_csv(tab / "tab0_splice_test.csv", index=False)

    # ------------------------------------------------------------ estimation
    step("gravity ladder, dose-response, hypotheses")
    hypotheses.sample_definition_effects(S).to_csv(tab / "tabB0_sample_definition_effects.csv", index=False)
    hypotheses.gravity_ladder(main_).to_csv(tab / "tabB1_gravity_ladder.csv", index=False)
    bins = pd.concat([hypotheses.dose_response(S.all10, "all 10"),
                      hypotheses.dose_response(main_, "9 (no HKG)")], ignore_index=True)
    bins.to_csv(tab / "tabA3_tariff_bins.csv", index=False)
    hyp = hypotheses.hypothesis_table(main_)
    hyp.to_csv(tab / "tabB2_hypotheses_with_fragility.csv", index=False)
    print(hyp[["hypothesis", "coef", "t_pair", "t_exp", "t_jack", "coef_drop_CHN"]]
          .round(3).to_string(index=False))
    R["hypotheses"] = hyp

    # ------------------------------------------------------------- inference
    if args.skip_ri:
        step("randomization inference SKIPPED -- power table uses the paper's recorded null sd")
        ri_all = PAPER_RI
    else:
        step(f"randomization inference, {args.n_ri} draws (x3)")
        ri_all = inference.ri_interaction(S.all10, n=args.n_ri, progress=True)
        ri_h1 = inference.ri_interaction(main_, n=args.n_ri)
        ri_h3 = inference.ri_triple(main_, n=args.n_ri)
        R["ri_h1_p"], R["ri_h3_p"] = ri_h1["p"], ri_h3["p"]
        json.dump({k: {kk: vv for kk, vv in v.items() if kk != "draws"}
                   for k, v in (("all10_H1", ri_all), ("main_H1", ri_h1), ("main_H3", ri_h3))},
                  open(tab / "tabB3_ri.json", "w"), indent=2, default=float)
    R.update(ri_sd=ri_all["sd"], ri_n_draws=ri_all["n_draws"])
    power = inference.power_table(ri_all["sd"], abs(b1), R["n_exporters"], scale["sd"])
    power.to_csv(tab / "tabA4_power.csv", index=False)
    print(power.round(3).to_string(index=False))
    R["power"] = power

    # ----------------------------------------------------- synthetic control
    step("synthetic control")
    donors = [e for e in S.exporters_main if e != TREATED_EXPORTER]
    sc_input = S.all10[S.all10.balanced | (S.all10.iso3_o != TREATED_EXPORTER)]
    sc = synth.synthetic_control(sc_input, donors=donors)
    placebos = {}
    if sc:
        sc["path"].to_csv(tab / "tabB4_synthetic_control_china.csv")
        print(sc["weights"][sc["weights"] > .01].round(3).to_string())
        print(f"  pre-RMSPE {sc['rmspe_pre']:.4f}")
        placebos = synth.in_space_placebos(S.all10, sc, sc["weights"].index, S.exporters_main)
    R["synth"] = sc

    # ------------------------------------------------- alternative designs
    step("Route A: payment corridors")
    C = corridors.add_corridor(main_)
    cor = corridors.corridor_effects(C)
    cor.to_csv(tab / "tabC1_corridor.csv", index=False)
    rr = corridors.reallocation_robustness(C)
    rr.to_csv(tab / "tabC2_corridor_robustness.csv", index=False)
    ev = corridors.event_study(C)
    ev.to_csv(tab / "tabC3_corridor_event_study.csv", index=False)
    print(cor.round(3).to_string(index=False))
    print(ev.round(3).to_string(index=False))

    if args.partner_findex and args.partner_findex.exists():
        step("Route B: bilateral complementarity")
        rb, ne_b = corridors.bilateral_complementarity(main_, pd.read_csv(args.partner_findex))
        rb["N_eff_tXmin"] = ne_b["N_eff"]
        rb.to_csv(tab / "tabC4_bilateral_fintech.csv", index=False)
        print(rb.round(4).to_string(index=False))

    # ----------------------------------------------------- specification curve
    step("specification curve")
    spec_all = spec_curve.specification_curve(S.all10)
    spec_all.to_csv(tab / "tabD1_specification_curve.csv", index=False)
    spec = spec_curve.identified(spec_all)
    R["spec_summary"] = spec_curve.summarise(spec_all)
    print(R["spec_summary"])

    # --------------------------------------------------------------- figures
    if figs:
        step("figures")
        savefig(plots.identification_anatomy(D, S.all10, var, stepwise), figdir,
                "figA1_identification_anatomy")
        savefig(plots.chn_hkg_contrast(D, S.all10), figdir, "figA2_chn_hkg_expose")
        savefig(plots.dose_response(bins), figdir, "figA3_tariff_bins")
        if ri_all["draws"].size:
            savefig(plots.power_gap(ri_all, power, abs(b1), scale["sd"], R["n_exporters"]),
                    figdir, "figA4_power_gap")
        loo = {lab: hypotheses.leave_one_out(main_, f, t) for lab, f, t in (
            ("H1 interaction", hypotheses.F_H1, "tXfin"),
            ("H2 payments", hypotheses.F_H2, "tXfin_pay"),
            ("H2 digital credit", hypotheses.F_H2, "tXfin_cred"),
            ("H3 triple", hypotheses.F_H3, "tXfinXinst"))}
        savefig(plots.fragility(loo), figdir, "figB1_fragility_all_hypotheses")
        if sc:
            savefig(plots.synthetic_control(sc, placebos), figdir, "figB2_synthetic_control_china")
        savefig(plots.corridor(ev, rr), figdir, "figC1_corridor")
        savefig(plots.specification_curve(spec, float(power.MDE.iloc[0]),
                                          list(spec_curve.SAMPLE_RULES),
                                          [v for v in INDEX_VARIANTS if v in set(spec.variant)],
                                          list(spec_curve.FE_STRUCTURES)),
                figdir, "figD1_specification_curve")

    # --------------------------------------------------------------- outputs
    numbers = {k: v for k, v in R.items()
               if isinstance(v, (int, float, np.floating, np.integer, dict)) and k != "synth"}
    numbers["power"] = power.to_dict("records")
    (tab / "numbers.json").write_text(json.dumps(numbers, indent=2, default=float))
    report.write_report(R, paths.root / "results" / "REPORT.md")
    step(f"done -> {paths.root / 'results'}")
    return R


if __name__ == "__main__":
    main()
