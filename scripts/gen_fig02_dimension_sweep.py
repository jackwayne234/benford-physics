#!/usr/bin/env python3
"""Fig 2: Dimension sweep — delta_B vs exponent n."""
import json, sys, os
sys.path.insert(0, os.path.dirname(__file__))
import plotly.graph_objects as go
from plotly_style import apply_dark_theme, add_conform_zone, save_figure

BASE = os.path.join(os.path.dirname(__file__), "..")
OUT = os.path.join(BASE, "results/round_trip/figures/02_dimension_sweep.html")

with open(os.path.join(BASE, "results/round_trip/dimension_sweep.json")) as f:
    data = json.load(f)

sweep = data["sweep"]
xs = [s["exponent"] for s in sweep]
ys = [s["delta_b"] for s in sweep]
verdicts = [s["verdict"] for s in sweep]

fig = go.Figure()

# Main curve
fig.add_trace(go.Scatter(
    x=xs, y=ys, mode="lines+markers",
    line=dict(color="#d2a8ff", width=2.5),
    marker=dict(size=9),
    name="\u03b4_B(n)",
    hovertemplate=(
        "<b>Exponent n = %{x}</b><br>"
        "\u03b4_B = %{y:.4f}<br>"
        "<extra></extra>"
    ),
))

# Mark n=0 (BE)
fig.add_trace(go.Scatter(
    x=[0], y=[sweep[0]["delta_b"]], mode="markers",
    marker=dict(symbol="square", size=14, color="#58a6ff", line=dict(color="white", width=2)),
    name=f"n=0: BE (\u03b4_B={sweep[0]['delta_b']:.4f})",
    hovertemplate="<b>Bose-Einstein (n=0)</b><br>\u03b4_B = %{y:.4f}<extra></extra>",
))

# Mark n=3 (Planck)
n3 = next(s for s in sweep if s["exponent"] == 3)
fig.add_trace(go.Scatter(
    x=[3], y=[n3["delta_b"]], mode="markers",
    marker=dict(symbol="diamond", size=14, color="#ffd700", line=dict(color="white", width=2)),
    name=f"n=3: Planck (\u03b4_B={n3['delta_b']:.4f})",
    hovertemplate="<b>Planck (n=3)</b><br>\u03b4_B = %{y:.4f}<extra></extra>",
))

# Inversion arrow: horizontal from delta_b to n=3
inv = data["inversion"]
fig.add_annotation(
    x=3, y=inv["planck_measured_delta_b"],
    ax=0.5, ay=inv["planck_measured_delta_b"],
    xref="x", yref="y", axref="x", ayref="y",
    showarrow=True, arrowhead=2, arrowsize=1.5,
    arrowcolor="#ffd700", arrowwidth=2,
)
fig.add_annotation(
    x=1.5, y=inv["planck_measured_delta_b"] + 0.002,
    text=f"<b>\u03b4_B = {inv['planck_measured_delta_b']:.4f} \u2192 n = {inv['recovered_exponent_exact']:.3f}</b>",
    showarrow=False, font=dict(color="#ffd700", size=14),
)

# Shaded area under curve
fig.add_trace(go.Scatter(
    x=xs + xs[::-1],
    y=ys + [0]*len(ys),
    fill="toself", fillcolor="rgba(210,168,255,0.08)",
    line=dict(width=0), showlegend=False, hoverinfo="skip",
))

add_conform_zone(fig)

apply_dark_theme(fig,
    height=600, width=1000,
    title=dict(
        text=(
            "<b>Dimension Sweep: B(x) = x<sup>n</sup> / (e<sup>x</sup> \u2212 1)</b><br>"
            "<span style='font-size:13px;color:#8b949e'>\u03b4_B = 0.028 inverts to n = 3 (exact) | Hover for details</span>"
        ),
        x=0.5, xanchor="center", font=dict(size=19),
    ),
    xaxis=dict(title="Exponent n  (spatial dimensions)", dtick=0.5),
    yaxis=dict(title="\u03b4_B  (Euclidean deviation)"),
)

save_figure(fig, OUT)
print("Fig 2 done.")
