#!/usr/bin/env python3
"""Fig 4: The Planck Wall — QG models through the Big Bang."""
import json, sys, os
sys.path.insert(0, os.path.dirname(__file__))
import plotly.graph_objects as go
from plotly_style import apply_dark_theme, add_conform_zone, MODEL_COLORS, save_figure

BASE = os.path.join(os.path.dirname(__file__), "..")
OUT = os.path.join(BASE, "results/round_trip/figures/04_planck_wall.html")

# Use extended data (10 models, 94 points) if available, else fall back
for fname in ["planck_wall_extended.json", "planck_wall_hires.json", "planck_wall.json"]:
    path = os.path.join(BASE, "results/round_trip", fname)
    if os.path.exists(path):
        with open(path) as f:
            data = json.load(f)
        break

models = data["models_tested"]
results = data["results"]
descs = data.get("model_descriptions", {})

fig = go.Figure()

for model in models:
    entries = results[model]
    xs = [e["T_planck_units"] for e in entries if e.get("computable", True)]
    ys = [e["delta_b"] for e in entries if e.get("computable", True)]
    color = MODEL_COLORS.get(model, "#888")
    desc = descs.get(model, "")

    fig.add_trace(go.Scatter(
        x=xs, y=ys, mode="lines+markers",
        line=dict(color=color, width=2.5),
        marker=dict(size=5),
        name=model,
        hovertemplate=(
            f"<b>{model}</b><br>"
            "T/T_P = %{x:.3f}<br>"
            "\u03b4_B = %{y:.5f}<br>"
            f"<span style='font-size:11px'>{desc}</span>"
            "<extra></extra>"
        ),
    ))

# Big Bang line at T=1
fig.add_vline(x=1.0, line_color="#ff4444", line_width=3)
fig.add_annotation(
    x=1.0, y=0.95, xref="x", yref="y domain",
    text="<b>BIG BANG</b><br>(T = T_Planck)",
    showarrow=False, font=dict(color="#ff4444", size=15),
    xanchor="left",
)

# Direction arrows
fig.add_annotation(
    x=0.003, y=-0.04, xref="x", yref="y",
    text="Today's universe<br>(T << T_P)",
    showarrow=False, font=dict(color="#ffd700", size=11),
)
fig.add_annotation(
    x=30, y=-0.04, xref="x", yref="y",
    text="Before the<br>Big Bang?",
    showarrow=False, font=dict(color="#7ee787", size=11),
)

add_conform_zone(fig)

apply_dark_theme(fig,
    height=750, width=1400,
    title=dict(
        text=(
            "<b>The Planck Wall: Which Quantum Gravity Models Survive the Big Bang?</b><br>"
            "<span style='font-size:13px;color:#8b949e'>"
            "Click legend to show/hide models | Double-click to isolate | Hover for values</span>"
        ),
        x=0.5, xanchor="center", font=dict(size=19),
    ),
    xaxis=dict(title="T / T_Planck  (higher T = further back in time)", type="log"),
    yaxis=dict(title="\u03b4_B  (Euclidean deviation from Benford)"),
)

save_figure(fig, OUT)
print("Fig 4 done.")
