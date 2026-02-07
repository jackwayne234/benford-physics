#!/usr/bin/env python3
"""Fig 1: epsilon(d) fingerprints for the 4 fundamental distributions."""
import json, sys, os
sys.path.insert(0, os.path.dirname(__file__))
from plotly.subplots import make_subplots
import plotly.graph_objects as go
from plotly_style import apply_dark_theme, DIST_COLORS, save_figure

BASE = os.path.join(os.path.dirname(__file__), "..")
OUT = os.path.join(BASE, "results/round_trip/figures/01_fingerprints.html")

DISTS = [
    ("Bose-Einstein",    "results/individual/bose_einstein_numerical.json"),
    ("Fermi-Dirac",      "results/individual/fermi_dirac_numerical.json"),
    ("Maxwell-Boltzmann", "results/individual/maxwell_boltzmann_numerical.json"),
    ("Planck",           "results/individual/planck_radiation_spectrum.json"),
]

fig = make_subplots(
    rows=2, cols=2,
    subplot_titles=["placeholder"] * 4,
    horizontal_spacing=0.1, vertical_spacing=0.12,
)

for i, (name, path) in enumerate(DISTS):
    row, col = divmod(i, 2)
    row += 1; col += 1
    with open(os.path.join(BASE, path)) as f:
        data = json.load(f)
    eps = data["per_digit_deviation"]
    delta_b = data["delta_b"]
    digits = list(range(1, 10))
    vals = [eps.get(str(d), eps.get(d, 0)) for d in digits]
    colors = [DIST_COLORS[name] if v >= 0 else DIST_COLORS[name] for v in vals]

    # Sign pattern string
    signs = "".join("+" if v >= 0 else "\u2212" for v in vals)

    fig.add_trace(go.Bar(
        x=digits, y=vals,
        marker_color=DIST_COLORS[name],
        marker_line_width=0,
        name=name,
        hovertemplate=(
            f"<b>{name}</b><br>"
            "Digit %{x}<br>"
            "\u03b5(d) = %{y:.5f}<extra></extra>"
        ),
        showlegend=True,
    ), row=row, col=col)

    # Update subplot title
    fig.layout.annotations[i].text = (
        f"<b>{name}</b>  (\u03b4_B = {delta_b:.4f})<br>"
        f"<span style='font-size:11px;color:#8b949e'>{signs}</span>"
    )

apply_dark_theme(fig,
    height=700, width=1100,
    title=dict(
        text=(
            "<b>\u03b5(d) Fingerprints \u2014 The Four Fundamental Distributions</b><br>"
            "<span style='font-size:13px;color:#8b949e'>Click legend to toggle | Hover for values</span>"
        ),
        x=0.5, xanchor="center", font=dict(size=20),
    ),
)

for i in range(4):
    row, col = divmod(i, 2)
    fig.update_xaxes(title_text="First Digit d", row=row+1, col=col+1, dtick=1)
    fig.update_yaxes(title_text="\u03b5(d)", row=row+1, col=col+1)
    fig.add_hline(y=0, line_dash="dot", line_color="#3a3a5a", line_width=1, row=row+1, col=col+1)

save_figure(fig, OUT)
print("Fig 1 done.")
