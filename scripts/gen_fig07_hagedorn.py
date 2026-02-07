#!/usr/bin/env python3
"""Fig 7: Hagedorn spotlight — String theory through the Big Bang."""
import json, sys, os
sys.path.insert(0, os.path.dirname(__file__))
import plotly.graph_objects as go
from plotly_style import apply_dark_theme, add_conform_zone, save_figure

BASE = os.path.join(os.path.dirname(__file__), "..")
OUT = os.path.join(BASE, "results/round_trip/figures/07_hagedorn_spotlight.html")

for fname in ["planck_wall_extended.json", "planck_wall_hires.json", "planck_wall.json"]:
    path = os.path.join(BASE, "results/round_trip", fname)
    if os.path.exists(path):
        with open(path) as f:
            data = json.load(f)
        break

hag = data["results"]["Hagedorn"]
xs = [e["T_planck_units"] for e in hag if e.get("computable", True)]
ys = [e["delta_b"] for e in hag if e.get("computable", True)]

fig = go.Figure()

# Main Hagedorn curve
fig.add_trace(go.Scatter(
    x=xs, y=ys, mode="lines+markers",
    line=dict(color="#ffd700", width=3),
    marker=dict(size=6, color="#ffd700"),
    name="Hagedorn (String Theory)",
    hovertemplate=(
        "<b>Hagedorn</b><br>"
        "T/T_P = %{x:.3f}<br>"
        "\u03b4_B = %{y:.5f}<extra></extra>"
    ),
))

# Find peak chaos
peak_idx = max(range(len(ys)), key=lambda i: ys[i])
fig.add_annotation(
    x=xs[peak_idx], y=ys[peak_idx],
    text=f"<b>Peak chaos</b><br>\u03b4_B = {ys[peak_idx]:.3f}",
    showarrow=True, arrowhead=2, arrowcolor="#ffd700",
    font=dict(color="#ffd700", size=13),
    ax=60, ay=-40,
)

# Post-wall settling (last few points)
post_wall = [(x, y) for x, y in zip(xs, ys) if x > 5]
if post_wall:
    avg_post = sum(y for _, y in post_wall) / len(post_wall)
    fig.add_annotation(
        x=post_wall[-1][0], y=avg_post,
        text=f"<b>Settles to \u03b4_B \u2248 {avg_post:.3f}</b><br>(near-perfect Benford)",
        showarrow=True, arrowhead=2, arrowcolor="#7ee787",
        font=dict(color="#7ee787", size=13),
        ax=-80, ay=-50,
    )

# Big Bang line
fig.add_vline(x=1.0, line_color="#ff4444", line_width=3)
fig.add_annotation(
    x=1.0, y=0.95, xref="x", yref="y domain",
    text="<b>BIG BANG</b><br>(T = T_Planck)",
    showarrow=False, font=dict(color="#ff4444", size=15),
    xanchor="left",
)

# Direction annotations
fig.add_annotation(
    x=0.003, y=0.01, xref="x", yref="y",
    text="<b>\u2190 direction of time</b>",
    showarrow=False, font=dict(color="#8b949e", size=11),
)
fig.add_annotation(
    x=30, y=0.01, xref="x", yref="y",
    text="<b>pre-Big-Bang era \u2192</b>",
    showarrow=False, font=dict(color="#7ee787", size=11),
)

add_conform_zone(fig)

apply_dark_theme(fig,
    height=650, width=1200,
    title=dict(
        text=(
            "<b>Hagedorn (String Theory): Through the Big Bang and Out the Other Side</b><br>"
            "<span style='font-size:13px;color:#8b949e'>"
            "The only model that emerges cleaner after the singularity | Hover for values</span>"
        ),
        x=0.5, xanchor="center", font=dict(size=18),
    ),
    xaxis=dict(title="T / T_Planck  (higher T = further back in time \u2192 through Big Bang \u2192 before)", type="log"),
    yaxis=dict(title="\u03b4_B  (Euclidean deviation)"),
)

save_figure(fig, OUT)
print("Fig 7 done.")
