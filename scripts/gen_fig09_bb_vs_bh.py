#!/usr/bin/env python3
"""Fig 9: Big Bang vs Black Hole — side-by-side comparison."""
import json, sys, os
sys.path.insert(0, os.path.dirname(__file__))
from plotly.subplots import make_subplots
import plotly.graph_objects as go
from plotly_style import apply_dark_theme, add_conform_zone, MODEL_COLORS, save_figure

BASE = os.path.join(os.path.dirname(__file__), "..")
OUT = os.path.join(BASE, "results/round_trip/figures/09_bb_vs_bh.html")

# Load both datasets
for fname in ["planck_wall_extended.json", "planck_wall_hires.json", "planck_wall.json"]:
    path = os.path.join(BASE, "results/round_trip", fname)
    if os.path.exists(path):
        with open(path) as f:
            bb_data = json.load(f)
        break

with open(os.path.join(BASE, "results/round_trip/black_hole_wall.json")) as f:
    bh_data = json.load(f)

SPOTLIGHT = ["Hagedorn", "Causal Set", "Standard", "LQG"]

fig = make_subplots(
    rows=1, cols=2,
    subplot_titles=[
        "<b>Big Bang Wall</b><br><span style='color:#ffd700;font-size:12px'>(T sweeps through T_Planck)</span>",
        "<b>Black Hole Wall</b><br><span style='color:#ffd700;font-size:12px'>(r sweeps through event horizon to singularity)</span>",
    ],
    horizontal_spacing=0.08,
)

for model in SPOTLIGHT:
    color = MODEL_COLORS.get(model, "#888")

    # Left: Big Bang
    if model in bb_data["results"]:
        entries = bb_data["results"][model]
        xs = [e["T_planck_units"] for e in entries if e.get("computable", True)]
        ys = [e["delta_b"] for e in entries if e.get("computable", True)]
        fig.add_trace(go.Scatter(
            x=xs, y=ys, mode="lines+markers",
            line=dict(color=color, width=2.5),
            marker=dict(size=5),
            name=model,
            legendgroup=model,
            hovertemplate=f"<b>{model}</b><br>T/T_P = %{{x:.3f}}<br>\u03b4_B = %{{y:.5f}}<extra></extra>",
        ), row=1, col=1)

    # Right: Black Hole
    if model in bh_data["infalling_observer"]:
        entries = bh_data["infalling_observer"][model]
        xs = [e["r_ratio"] for e in entries if e.get("computable", True)]
        ys = [e["delta_b"] for e in entries if e.get("computable", True)]
        fig.add_trace(go.Scatter(
            x=xs, y=ys, mode="lines+markers",
            line=dict(color=color, width=2.5),
            marker=dict(size=5),
            name=model,
            legendgroup=model,
            showlegend=False,
            hovertemplate=f"<b>{model}</b><br>r/r_s = %{{x:.3f}}<br>\u03b4_B = %{{y:.5f}}<extra></extra>",
        ), row=1, col=2)

# Big Bang line on left
fig.add_vline(x=1.0, line_color="#ff4444", line_width=2, row=1, col=1)
fig.add_annotation(x=1.0, y=0.45, text="<b>BIG<br>BANG</b>", showarrow=False,
                   font=dict(color="#ff4444", size=12), xref="x", yref="y")

# Event horizon + singularity on right
fig.add_vline(x=1.0, line_color="#ff4444", line_width=2, row=1, col=2)
fig.add_annotation(x=1.0, y=0.45, text="<b>EVENT<br>HORIZON</b>", showarrow=False,
                   font=dict(color="#ff4444", size=11), xref="x2", yref="y2")
fig.add_vline(x=0.0, line_color="#ffffff", line_width=1.5, row=1, col=2)
fig.add_annotation(x=0.0, y=0.45, text="<b>SINGULARITY</b>", showarrow=False,
                   font=dict(color="#ffffff", size=11), xref="x2", yref="y2")

add_conform_zone(fig, row=1, col=1)
add_conform_zone(fig, row=1, col=2)

fig.update_xaxes(title_text="T / T_Planck", type="log", row=1, col=1)
fig.update_xaxes(title_text="r / r_s  (\u2192 horizon \u2192 singularity)", autorange="reversed", row=1, col=2)
fig.update_yaxes(title_text="\u03b4_B  (Euclidean deviation)", row=1, col=1)
fig.update_yaxes(title_text="\u03b4_B", row=1, col=2)

apply_dark_theme(fig,
    height=700, width=1500,
    title=dict(
        text=(
            "<b>Same Physics, Two Walls: Big Bang vs Black Hole</b><br>"
            "<span style='font-size:13px;color:#8b949e'>"
            "Do QG models behave the same at both singularities? | Click legend to toggle</span>"
        ),
        x=0.5, xanchor="center", font=dict(size=19),
    ),
)

save_figure(fig, OUT)
print("Fig 9 done.")
