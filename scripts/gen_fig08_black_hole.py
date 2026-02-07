#!/usr/bin/env python3
"""Fig 8: Black Hole Wall — 10 QG models through the event horizon and singularity."""
import json, sys, os
sys.path.insert(0, os.path.dirname(__file__))
import plotly.graph_objects as go
from plotly_style import apply_dark_theme, add_conform_zone, MODEL_COLORS, save_figure

BASE = os.path.join(os.path.dirname(__file__), "..")
OUT = os.path.join(BASE, "results/round_trip/figures/08_black_hole_wall.html")

with open(os.path.join(BASE, "results/round_trip/black_hole_wall.json")) as f:
    data = json.load(f)

models = data["models_tested"]
observer = data["infalling_observer"]
descs = {
    "Standard": "GR+QFT, no modification",
    "LQG": "Polymer dispersion, Brillouin zone",
    "GUP": "Minimum length",
    "DSR": "Energy saturation",
    "Hagedorn": "String theory, exponential state growth",
    "Causal Set": "Discrete random spacetime, Gaussian UV suppression",
    "Asym. Safety": "Running spectral dimension 4\u21922",
    "Horava-Lif.": "Anisotropic scaling",
    "Noncommut.": "Minimum area",
    "CDT": "Sharp dimensional reduction 4\u21922",
}

fig = go.Figure()

for model in models:
    entries = observer[model]
    xs = [e["r_ratio"] for e in entries if e.get("computable", True)]
    ys = [e["delta_b"] for e in entries if e.get("computable", True)]
    zones = [e.get("zone", "") for e in entries if e.get("computable", True)]
    color = MODEL_COLORS.get(model, "#888")
    desc = descs.get(model, "")

    fig.add_trace(go.Scatter(
        x=xs, y=ys, mode="lines+markers",
        line=dict(color=color, width=2),
        marker=dict(size=5),
        name=model,
        visible="legendonly" if model not in ["Causal Set", "Hagedorn", "Standard"] else True,
        hovertemplate=(
            f"<b>{model}</b><br>"
            "r/r_s = %{x:.3f}<br>"
            "\u03b4_B = %{y:.5f}<br>"
            f"<span style='font-size:11px'>{desc}</span>"
            "<extra></extra>"
        ),
    ))

# Event Horizon line at r=1
fig.add_vline(x=1.0, line_color="#ff4444", line_width=3)
fig.add_annotation(
    x=1.0, y=0.95, xref="x", yref="y domain",
    text="<b>EVENT<br>HORIZON</b>",
    showarrow=False, font=dict(color="#ff4444", size=14),
    xanchor="left",
)

# Singularity line at r=0
fig.add_vline(x=0.0, line_color="#ffffff", line_width=2)
fig.add_annotation(
    x=0.0, y=0.95, xref="x", yref="y domain",
    text="<b>SINGULARITY</b>",
    showarrow=False, font=dict(color="#ffffff", size=13),
    xanchor="left",
)

# Zone labels
fig.add_annotation(x=5, y=-0.03, text="Outside<br>(approaching BH)",
                   showarrow=False, font=dict(color="#8b949e", size=11))
fig.add_annotation(x=0.5, y=-0.03, text="Inside<br>(falling to singularity)",
                   showarrow=False, font=dict(color="#ff8888", size=11))
fig.add_annotation(x=-0.5, y=-0.03, text="Other Side?<br>(post-singularity bounce)",
                   showarrow=False, font=dict(color="#7ee787", size=11))

add_conform_zone(fig)

apply_dark_theme(fig,
    height=750, width=1500,
    title=dict(
        text=(
            "<b>Black Hole: Through the Horizon, Through the Singularity, Out the Other Side</b><br>"
            "<span style='font-size:13px;color:#8b949e'>"
            "Click legend to show/hide models | Double-click to isolate | Hover for values</span>"
        ),
        x=0.5, xanchor="center", font=dict(size=18),
    ),
    xaxis=dict(
        title="r / r_s  (far away \u2192 horizon \u2192 singularity \u2192 other side)",
        autorange="reversed",
    ),
    yaxis=dict(title="\u03b4_B  (Euclidean deviation from Benford)"),
)

save_figure(fig, OUT)
print("Fig 8 done.")
