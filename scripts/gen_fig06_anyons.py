#!/usr/bin/env python3
"""Fig 6: Anyon interpolation — epsilon(d) from BE to FD."""
import json, sys, os
sys.path.insert(0, os.path.dirname(__file__))
from plotly.subplots import make_subplots
import plotly.graph_objects as go
from plotly_style import apply_dark_theme, save_figure

BASE = os.path.join(os.path.dirname(__file__), "..")
OUT = os.path.join(BASE, "results/round_trip/figures/06_anyons.html")

with open(os.path.join(BASE, "results/round_trip/whiteboard.json")) as f:
    data = json.load(f)

# Filter anyon candidates and sort by g value
anyons = [c for c in data["candidates"] if "Anyon" in c["name"] and c["epsilon_d"]]
g_values = {"g=0": 0, "g=0.25": 0.25, "g=0.5": 0.5, "g=0.75": 0.75, "g=1.0": 1.0}
anyons.sort(key=lambda c: next((v for k, v in g_values.items() if k in c["name"]), 99))

# Color gradient blue -> purple -> red
grad_colors = ["#58a6ff", "#8888dd", "#bb77bb", "#cc6699", "#f78166"]
labels = ["boson", "1/4-on", "semion", "3/4-on", "fermion"]

fig = make_subplots(
    rows=1, cols=5,
    subplot_titles=["placeholder"] * 5,
    horizontal_spacing=0.04,
    shared_yaxes=True,
)

digits = list(range(1, 10))

for i, anyon in enumerate(anyons[:5]):
    eps = anyon["epsilon_d"]
    vals = [eps.get(str(d), 0) for d in digits]
    delta_b = anyon["delta_b"]

    # Extract g value from name
    g_str = anyon["name"].split("g=")[1].split(" ")[0] if "g=" in anyon["name"] else "?"

    fig.add_trace(go.Bar(
        x=digits, y=vals,
        marker_color=grad_colors[i],
        marker_line_width=0,
        name=f"g={g_str}",
        hovertemplate=(
            f"<b>g={g_str} ({labels[i]})</b><br>"
            "Digit %{x}<br>"
            "\u03b5(d) = %{y:.5f}<extra></extra>"
        ),
    ), row=1, col=i+1)

    fig.add_hline(y=0, line_dash="dot", line_color="#3a3a5a", line_width=1, row=1, col=i+1)

    fig.layout.annotations[i].text = (
        f"<span style='color:{grad_colors[i]}'><b>g={g_str}</b><br>"
        f"({labels[i]})<br>"
        f"\u03b4_B={delta_b:.4f}</span>"
    )

fig.update_xaxes(dtick=1, title_text="Digit d")
fig.update_yaxes(title_text="\u03b5(d)", col=1)

apply_dark_theme(fig,
    height=500, width=1500,
    title=dict(
        text=(
            "<b>Anyon Interpolation: \u03b5(d) evolves smoothly from BE \u2192 FD</b><br>"
            "<span style='font-size:13px;color:#8b949e'>Click legend to toggle | Hover for values</span>"
        ),
        x=0.5, xanchor="center", font=dict(size=19),
    ),
)

save_figure(fig, OUT)
print("Fig 6 done.")
