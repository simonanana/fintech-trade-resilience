#!/usr/bin/env python3
"""Write a synthetic panel with the production schema for demos and tests.

    python scripts/make_synthetic_panel.py            # -> data/processed/panel_synthetic.parquet
    python scripts/run_pipeline.py --panel data/processed/panel_synthetic.parquet --out demo_output --skip-ri

The synthetic data reproduce the structure of the problem (one heavily targeted
exporter), not the paper's estimates.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ftres.config import Paths  # noqa: E402
from ftres.synthetic import make_panel  # noqa: E402

if __name__ == "__main__":
    out = Paths().data_processed / "panel_synthetic.parquet"
    out.parent.mkdir(parents=True, exist_ok=True)
    df = make_panel()
    df.to_parquet(out, index=False)
    print(f"wrote {out} {df.shape}")
