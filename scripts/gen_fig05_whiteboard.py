#!/usr/bin/env python3
"""Fig 5: The Whiteboard — Exotic physics existence filter."""
import json, sys, os
sys.path.insert(0, os.path.dirname(__file__))
import plotly.graph_objects as go
from plotly_style import apply_dark_theme, save_figure

BASE = os.path.join(os.path.dirname(__file__), "..")
OUT = os.path.join(BASE, "results/round_trip/figures/05_whiteboard.html")

with open(os.path.join(BASE, "results/round_trip/whiteboard.json")) as f:
    data = json.load(f)

# Separate existing (computable) vs undefined
exist = [c for c in data["candidates"] if c["status"] == "EXISTS" and c["delta_b"] is not None]
undefined = [c for c in data["candidates"] if c["status"] == "UNDEFINED"]

# Sort by delta_b ascending (best conformance first)
exist.sort(key=lambda c: c["delta_b"])

names = [c["name"] for c in exist]
deltas = [c["delta_b"] for c in exist]
categories = [c.get("category", "") for c in exist]
verdicts = [c.get("verdict", "") for c in exist]

# Color by category
CAT_COLORS = {
    "Fractional Statistics": "#58a6ff",
    "Exotic Mass":          "#ff4444",
    "Black Holes":          "#ffd700",
    "Dark Matter":          "#7ee787",
    "Dark Energy":          "#ff69b4",
    "BSM Particles":        "#d2a8ff",
    "Relativistic QFT":     "#00d4ff",
    "Quantum Gravity":      "#f78166",
}
colors = [CAT_COLORS.get(c, "#888") for c in categories]

fig = go.Figure()

fig.add_trace(go.Bar(
    y=names, x=deltas,
    orientation="h",
    marker_color=colors,
    hovertemplate=(
        "<b>%{y}</b><br>"
        "\u03b4_B = %{x:.4f}<br>"
        "<extra></extra>"
    ),
    text=[f" {d:.4f}" for d in deltas],
    textposition="outside",
    textfont=dict(size=11, color="#c9d1d9"),
))

# Threshold lines
fig.add_vline(x=0.025, line_dash="dash", line_color="#7ee787", line_width=2, opacity=0.6,
              annotation_text="CONFORMS threshold", annotation_position="top",
              annotation_font=dict(color="#7ee787", size=11))
fig.add_vline(x=0.1, line_dash="dash", line_color="#ff4444", line_width=2, opacity=0.6,
              annotation_text="DEVIATES threshold", annotation_position="top",
              annotation_font=dict(color="#ff4444", size=11))

# Add undefined candidates as annotation
if undefined:
    undef_text = "<b>UNDEFINED (cannot exist):</b><br>" + "<br>".join(
        f"\u2022 {c['name']}: {c.get('reason', 'N/A')}" for c in undefined
    )
    fig.add_annotation(
        x=0.5, y=-0.08, xref="paper", yref="paper",
        text=undef_text, showarrow=False,
        font=dict(color="#ff6666", size=10),
        align="left", bgcolor="rgba(255,68,68,0.08)",
        bordercolor="#ff4444", borderwidth=1, borderpad=8,
    )

# Category legend as color key
for cat, color in CAT_COLORS.items():
    if cat in categories:
        fig.add_trace(go.Scatter(
            x=[None], y=[None], mode="markers",
            marker=dict(size=10, color=color),
            name=cat, showlegend=True,
        ))

apply_dark_theme(fig,
    height=max(600, len(names) * 32 + 200),
    width=1000,
    title=dict(
        text=(
            "<b>The Whiteboard: Exotic Physics Existence Filter</b><br>"
            "<span style='font-size:13px;color:#8b949e'>All candidates ranked by \u03b4_B | Hover for details</span>"
        ),
        x=0.5, xanchor="center", font=dict(size=19),
    ),
    xaxis=dict(title="\u03b4_B  (Euclidean deviation from Benford)"),
    yaxis=dict(autorange="reversed"),
    margin=dict(l=250, r=80, t=100, b=120),
)

save_figure(fig, OUT)
print("Fig 5 done.")
