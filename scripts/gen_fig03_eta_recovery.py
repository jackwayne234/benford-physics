#!/usr/bin/env python3
"""Fig 3: Eta recovery — delta_B vs alpha."""
import json, sys, os
sys.path.insert(0, os.path.dirname(__file__))
import plotly.graph_objects as go
from plotly_style import apply_dark_theme, add_conform_zone, save_figure

BASE = os.path.join(os.path.dirname(__file__), "..")
OUT = os.path.join(BASE, "results/round_trip/figures/03_eta_recovery.html")

with open(os.path.join(BASE, "results/round_trip/eta_recovery.json")) as f:
    data = json.load(f)

sweep = data["sweep"]
xs = [s["alpha"] for s in sweep]
ys = [s["delta_b"] for s in sweep]

fig = go.Figure()

# Main curve
fig.add_trace(go.Scatter(
    x=xs, y=ys, mode="lines+markers",
    line=dict(color="#f78166", width=2.5),
    marker=dict(size=9),
    name="\u03b4_B(\u03b1)",
    hovertemplate=(
        "<b>\u03b1 = %{x:.1f}</b><br>"
        "\u03b4_B = %{y:.6f}<extra></extra>"
    ),
))

# Shaded area
fig.add_trace(go.Scatter(
    x=xs + xs[::-1], y=ys + [0]*len(ys),
    fill="toself", fillcolor="rgba(247,129,102,0.08)",
    line=dict(width=0), showlegend=False, hoverinfo="skip",
))

# Mark alpha=0 (BE)
be = sweep[0]
fig.add_trace(go.Scatter(
    x=[0], y=[be["delta_b"]], mode="markers",
    marker=dict(symbol="square", size=14, color="#58a6ff", line=dict(color="white", width=2)),
    name=f"\u03b1=0: BE (\u03b4_B={be['delta_b']:.4f})",
    hovertemplate="<b>Bose-Einstein (\u03b1=0)</b><br>\u03b4_B = %{y:.6f}<extra></extra>",
))

# Mark alpha=1 (FD)
fd = next(s for s in sweep if s["alpha"] == 1.0)
fig.add_trace(go.Scatter(
    x=[1.0], y=[fd["delta_b"]], mode="markers",
    marker=dict(symbol="diamond", size=14, color="#f78166", line=dict(color="white", width=2)),
    name=f"\u03b1=1: FD (\u03b4_B={fd['delta_b']:.4f})",
    hovertemplate="<b>Fermi-Dirac (\u03b1=1)</b><br>\u03b4_B = %{y:.6f}<extra></extra>",
))

# Inversion annotation
inv = data["inversion"]
fig.add_annotation(
    x=0.5, y=inv["target_delta_b"] - 0.001,
    text=(
        f"<b>\u03b4_B = {inv['target_delta_b']:.4f} \u2192 \u03b1 = {inv['recovered_alpha']:.3f}</b><br>"
        f"\u03b7(1) = ln(2) = {inv['eta_1_value']:.4f}"
    ),
    showarrow=False, font=dict(color="#f78166", size=13),
)

# Vertical line at alpha=1
fig.add_vline(x=1.0, line_dash="dot", line_color="#8b949e", line_width=1, opacity=0.5)

# Horizontal line at FD delta_b
fig.add_hline(y=fd["delta_b"], line_dash="dot", line_color="#f78166", line_width=1, opacity=0.3)

apply_dark_theme(fig,
    height=600, width=1000,
    title=dict(
        text=(
            "<b>Eta Recovery: n(x) = 1/(e<sup>x</sup>\u22121) \u2212 \u03b1\u00b72/(e<sup>2x</sup>\u22121)</b><br>"
            "<span style='font-size:13px;color:#8b949e'>\u03b4_B = 0.0117 inverts to \u03b1 = 1 (exact), confirming \u03b7(1) = ln(2) | Hover for details</span>"
        ),
        x=0.5, xanchor="center", font=dict(size=18),
    ),
    xaxis=dict(title="\u03b1  (Dirichlet eta modulation strength)", dtick=0.2),
    yaxis=dict(title="\u03b4_B  (Euclidean deviation)"),
)

save_figure(fig, OUT)
print("Fig 3 done.")
