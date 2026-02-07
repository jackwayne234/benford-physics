#!/usr/bin/env python3
"""Fig 10: Hawking Radiation vs Causal Set Theory — multi-panel analysis."""
import json, sys, os
sys.path.insert(0, os.path.dirname(__file__))
from plotly.subplots import make_subplots
import plotly.graph_objects as go
from plotly_style import apply_dark_theme, add_conform_zone, save_figure

BASE = os.path.join(os.path.dirname(__file__), "..")
OUT = os.path.join(BASE, "results/round_trip/figures/10_hawking_vs_causal_set.html")

with open(os.path.join(BASE, "results/round_trip/black_hole_wall.json")) as f:
    bh_data = json.load(f)

with open(os.path.join(BASE, "results/round_trip/whiteboard.json")) as f:
    wb_data = json.load(f)

# Get Causal Set journey data
cs_entries = bh_data["infalling_observer"]["Causal Set"]
cs_xs = [e["r_ratio"] for e in cs_entries if e.get("computable", True)]
cs_ys = [e["delta_b"] for e in cs_entries if e.get("computable", True)]
cs_eps = [e["epsilon_d"] for e in cs_entries if e.get("computable", True)]

# Get Hawking candidates
hawking = [c for c in wb_data["candidates"] if "Hawking" in c["name"] and c["delta_b"] is not None]
hawking.sort(key=lambda c: c["delta_b"])

hawking_colors = ["#ff7700", "#ffaa00", "#ffd700", "#ff4400"]

fig = make_subplots(
    rows=2, cols=2,
    subplot_titles=[
        "<b>Causal Set Journey + Hawking Reference Levels</b>",
        "<b>Causal Set Fingerprints at 3 Positions</b>",
        "<b>Full Scale \u2014 Singularity Spike</b>",
        "<b>Hawking Radiation Fingerprints at 4 Greybody Cutoffs</b>",
    ],
    horizontal_spacing=0.1, vertical_spacing=0.15,
    specs=[[{"type": "scatter"}, {"type": "scatter"}],
           [{"type": "scatter"}, {"type": "scatter"}]],
)

# ---- TOP LEFT: Causal Set journey with Hawking reference lines ----
fig.add_trace(go.Scatter(
    x=cs_xs, y=cs_ys, mode="lines+markers",
    line=dict(color="#00d4ff", width=2.5),
    marker=dict(size=5),
    name="Causal Set",
    hovertemplate="<b>Causal Set</b><br>r/r_s = %{x:.3f}<br>\u03b4_B = %{y:.5f}<extra></extra>",
), row=1, col=1)

for i, h in enumerate(hawking):
    omega = h["name"].split("=")[1] if "=" in h["name"] else "?"
    fig.add_hline(
        y=h["delta_b"], line_dash="dash",
        line_color=hawking_colors[i % len(hawking_colors)], line_width=2,
        row=1, col=1,
        annotation_text=f"\u03c9_c={omega}  \u03b4_B={h['delta_b']:.4f}",
        annotation_font=dict(color=hawking_colors[i % len(hawking_colors)], size=10),
        annotation_position="top right",
    )

fig.add_vline(x=1.0, line_color="#ff4444", line_width=2, row=1, col=1)
fig.add_vline(x=0.0, line_color="#ffffff", line_width=1.5, row=1, col=1)
add_conform_zone(fig, row=1, col=1)

# ---- TOP RIGHT: Causal Set fingerprints at 3 positions ----
digits = list(range(1, 10))
positions = []
# Far out (r ~ 5), near horizon (r ~ 1.1), deep inside (r ~ 0.1)
targets = [5.0, 1.1, 0.1]
labels = ["Far out (r=5)", "Near horizon (r\u22481.1)", "Deep inside (r\u22480.1)"]
fp_colors = ["#00d4ff", "#ffd700", "#ff6666"]

for t_val in targets:
    best_idx = min(range(len(cs_xs)), key=lambda j: abs(cs_xs[j] - t_val))
    positions.append(best_idx)

for k, (idx, label, color) in enumerate(zip(positions, labels, fp_colors)):
    if idx < len(cs_eps) and cs_eps[idx]:
        eps = cs_eps[idx]
        vals = [eps.get(str(d), 0) for d in digits]
        fig.add_trace(go.Scatter(
            x=digits, y=vals, mode="lines+markers",
            line=dict(color=color, width=2),
            marker=dict(size=7),
            name=f"CS: {label}  \u03b4_B={cs_ys[idx]:.4f}",
            hovertemplate=f"<b>{label}</b><br>Digit %{{x}}<br>\u03b5(d) = %{{y:.5f}}<extra></extra>",
        ), row=1, col=2)

fig.add_hline(y=0, line_dash="dot", line_color="#3a3a5a", line_width=1, row=1, col=2)

# ---- BOTTOM LEFT: Full scale Causal Set with singularity spike ----
fig.add_trace(go.Scatter(
    x=cs_xs, y=cs_ys, mode="lines+markers",
    line=dict(color="#00d4ff", width=2),
    marker=dict(size=4),
    name="Causal Set (full)",
    showlegend=False,
    hovertemplate="<b>Causal Set</b><br>r/r_s = %{x:.3f}<br>\u03b4_B = %{y:.5f}<extra></extra>",
), row=2, col=1)

fig.add_vline(x=1.0, line_color="#ff4444", line_width=2, row=2, col=1)
fig.add_vline(x=0.0, line_color="#ffffff", line_width=1.5, row=2, col=1)
add_conform_zone(fig, row=2, col=1)

# ---- BOTTOM RIGHT: Hawking fingerprints at 4 greybody cutoffs ----
for i, h in enumerate(hawking):
    eps = h["epsilon_d"]
    vals = [eps.get(str(d), 0) for d in digits]
    omega = h["name"].split("=")[1] if "=" in h["name"] else "?"
    fig.add_trace(go.Scatter(
        x=digits, y=vals, mode="lines+markers",
        line=dict(color=hawking_colors[i % len(hawking_colors)], width=2),
        marker=dict(size=7),
        name=f"Hawking \u03c9_c={omega}  \u03b4_B={h['delta_b']:.4f}",
        hovertemplate=f"<b>Hawking \u03c9_c={omega}</b><br>Digit %{{x}}<br>\u03b5(d) = %{{y:.5f}}<extra></extra>",
    ), row=2, col=2)

fig.add_hline(y=0, line_dash="dot", line_color="#3a3a5a", line_width=1, row=2, col=2)

# Axis labels
fig.update_xaxes(title_text="r / r_s", autorange="reversed", row=1, col=1)
fig.update_xaxes(title_text="First Digit d", dtick=1, row=1, col=2)
fig.update_xaxes(title_text="r / r_s", autorange="reversed", row=2, col=1)
fig.update_xaxes(title_text="First Digit d", dtick=1, row=2, col=2)
fig.update_yaxes(title_text="\u03b4_B", row=1, col=1)
fig.update_yaxes(title_text="\u03b5(d)", row=1, col=2)
fig.update_yaxes(title_text="\u03b4_B", row=2, col=1)
fig.update_yaxes(title_text="\u03b5(d)", row=2, col=2)

apply_dark_theme(fig,
    height=1000, width=1500,
    title=dict(
        text=(
            "<b>Hawking Radiation vs Causal Set Theory</b><br>"
            "<span style='font-size:13px;color:#8b949e'>"
            "Click legend to toggle traces | Hover for precise values</span>"
        ),
        x=0.5, xanchor="center", font=dict(size=19),
    ),
    legend=dict(font=dict(size=11)),
)

save_figure(fig, OUT)
print("Fig 10 done.")
