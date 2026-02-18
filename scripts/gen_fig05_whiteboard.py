#!/usr/bin/env python3
"""Fig 5: The Whiteboard — Exotic physics existence filter.
Clean black-background version with values inside bars."""
import json, sys, os
sys.path.insert(0, os.path.dirname(__file__))
import plotly.graph_objects as go
from plotly_style import save_figure

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

# Build text labels: value inside bar for long bars, outside for short ones
# Use white text inside colored bars, bright text outside on black
max_delta = max(deltas)
text_labels = []
text_positions = []
text_colors = []
for d in deltas:
    text_labels.append(f"  {d:.4f}  ")
    if d > max_delta * 0.15:
        text_positions.append("inside")
        text_colors.append("#ffffff")
    else:
        text_positions.append("outside")
        text_colors.append("#ffffff")

fig = go.Figure()

fig.add_trace(go.Bar(
    y=names, x=deltas,
    orientation="h",
    marker=dict(
        color=colors,
        line=dict(color="rgba(255,255,255,0.1)", width=1),
    ),
    hovertemplate=(
        "<b>%{y}</b><br>"
        "\u03b4_B = %{x:.4f}<br>"
        "<extra></extra>"
    ),
    text=text_labels,
    textposition=text_positions,
    textfont=dict(size=13, color=text_colors, family="monospace"),
    width=0.7,
))

# Threshold lines
fig.add_vline(x=0.025, line_dash="dash", line_color="#00ff66", line_width=2, opacity=0.7)
fig.add_vline(x=0.1, line_dash="dash", line_color="#ff3333", line_width=2, opacity=0.7)

# Threshold labels as annotations at the top
fig.add_annotation(
    x=0.025, y=1.02, xref="x", yref="paper",
    text="<b>CONFORMS</b>", showarrow=False,
    font=dict(color="#00ff66", size=12, family="monospace"),
    xanchor="center",
)
fig.add_annotation(
    x=0.1, y=1.02, xref="x", yref="paper",
    text="<b>DEVIATES</b>", showarrow=False,
    font=dict(color="#ff3333", size=12, family="monospace"),
    xanchor="center",
)

# Undefined candidates box at bottom
if undefined:
    undef_lines = []
    for c in undefined:
        undef_lines.append(f"\u2022 {c['name']}")
    undef_text = "<b>UNDEFINED (cannot exist):</b>  " + "  |  ".join(undef_lines)
    fig.add_annotation(
        x=0.5, y=-0.12, xref="paper", yref="paper",
        text=undef_text, showarrow=False,
        font=dict(color="#ff6666", size=11, family="monospace"),
        align="center", bgcolor="rgba(255,50,50,0.08)",
        bordercolor="#ff4444", borderwidth=1, borderpad=10,
    )

# Category legend
for cat, color in CAT_COLORS.items():
    if cat in categories:
        fig.add_trace(go.Scatter(
            x=[None], y=[None], mode="markers",
            marker=dict(size=12, color=color, symbol="square"),
            name=cat, showlegend=True,
        ))

# Layout — black background, clean and readable
n_bars = len(names)
fig.update_layout(
    paper_bgcolor="#000000",
    plot_bgcolor="#000000",
    font=dict(family="monospace", size=14, color="#ffffff"),
    height=max(750, n_bars * 44 + 250),
    width=1400,
    title=dict(
        text=(
            "<b>The Whiteboard: Exotic Physics Existence Filter</b><br>"
            "<span style='font-size:14px;color:#888888'>"
            "All candidates ranked by \u03b4_B  |  Hover for details</span>"
        ),
        x=0.5, xanchor="center",
        font=dict(size=20, color="#ffffff"),
    ),
    xaxis=dict(
        title=dict(
            text="\u03b4_B  (Euclidean deviation from Benford)",
            font=dict(size=15, color="#aaaaaa"),
        ),
        gridcolor="#1a1a1a",
        zeroline=False,
        tickfont=dict(color="#aaaaaa", size=13),
    ),
    yaxis=dict(
        autorange="reversed",
        gridcolor="#1a1a1a",
        zeroline=False,
        tickfont=dict(color="#ffffff", size=14),
    ),
    margin=dict(l=300, r=80, t=110, b=160),
    legend=dict(
        font=dict(size=13, color="#222222"),
        bgcolor="rgba(255,255,255,0.92)",
        bordercolor="#aaaaaa",
        borderwidth=1,
        x=0.98, y=0.5, xanchor="right",
        itemclick="toggle",
        itemdoubleclick="toggleothers",
    ),
    bargap=0.25,
)

save_figure(fig, OUT)
print("Fig 5 done.")
