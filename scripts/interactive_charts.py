#!/usr/bin/env python3
"""Generate three interactive HTML charts:
  10a — Black Hole:  Causal Set + Hawking radiation (others hidden but toggleable)
  10b — Big Bang:    All 10 QG models visible (toggleable)
  10c — Wormhole:    All 10 QG models + Pure Casimir through a Morris-Thorne throat

Click any legend entry to show/hide that trace.
"""

import json
import os
import sys
import plotly.graph_objects as go

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

FIG_DIR = "results/round_trip/figures"


def load_json(path):
    with open(path) as f:
        return json.load(f)


# ── Shared style ──
DARK_BG = "#0d1117"
PANEL_BG = "#161b22"
GRID_COLOR = "#21262d"
TEXT_COLOR = "#c9d1d9"
MUTED = "#8b949e"

MODEL_COLORS = {
    "Standard":     "#58a6ff",
    "LQG":          "#7ee787",
    "GUP":          "#f78166",
    "DSR":          "#d2a8ff",
    "Hagedorn":     "#ffd700",
    "Causal Set":   "#00d4ff",
    "Asym. Safety": "#ff69b4",
    "Horava-Lif.":  "#ff6347",
    "Noncommut.":   "#98fb98",
    "CDT":          "#dda0dd",
}

HAWK_COLORS = {
    "0.5": "#ffd700",
    "1.0": "#ffaa00",
    "2.0": "#ff7700",
    "5.0": "#ff4400",
}


def dark_layout(title, xaxis_title, yaxis_title, height=700):
    return go.Layout(
        title=dict(text=title, font=dict(size=22, color=TEXT_COLOR),
                   x=0.5, xanchor="center"),
        paper_bgcolor=DARK_BG,
        plot_bgcolor=PANEL_BG,
        font=dict(family="monospace", size=14, color=TEXT_COLOR),
        xaxis=dict(
            title=dict(text=xaxis_title, font=dict(size=16)),
            gridcolor=GRID_COLOR, zeroline=False,
            tickfont=dict(size=14, color=MUTED),
        ),
        yaxis=dict(
            title=dict(text=yaxis_title, font=dict(size=16)),
            gridcolor=GRID_COLOR, zeroline=False,
            tickfont=dict(size=14, color=MUTED),
        ),
        legend=dict(
            font=dict(size=14), bgcolor="rgba(22,27,34,0.8)",
            bordercolor="#30363d", borderwidth=1,
            itemclick="toggle", itemdoubleclick="toggleothers",
        ),
        height=height,
        margin=dict(l=80, r=40, t=100, b=80),
    )


# ═══════════════════════════════════════════════════════════════════════════
# CHART 10a — BLACK HOLE
# ═══════════════════════════════════════════════════════════════════════════

def chart_black_hole():
    bh_data = load_json("results/round_trip/black_hole_wall.json")
    wb_data = load_json("results/round_trip/whiteboard.json")

    # Hawking radiation from whiteboard
    hawking = {}
    for c in wb_data["candidates"]:
        if "Hawking" in c.get("name", "") and c["status"] == "EXISTS":
            hawking[c["name"]] = c

    fig = go.Figure()

    # ── All 10 QG models through the BH ──
    # Default: only Causal Set visible, rest hidden but toggleable
    for model_name in bh_data["models_tested"]:
        entries = bh_data["infalling_observer"][model_name]
        pts = [(e["r_ratio"], e["delta_b"]) for e in entries if e.get("computable")]
        pts.sort(key=lambda p: -p[0])
        if not pts:
            continue
        rs, dbs = zip(*pts)

        visible = model_name == "Causal Set"
        fig.add_trace(go.Scatter(
            x=list(rs), y=list(dbs),
            mode="lines+markers",
            name=model_name,
            line=dict(color=MODEL_COLORS.get(model_name, MUTED), width=3 if visible else 2),
            marker=dict(size=6),
            visible=True if visible else "legendonly",
            hovertemplate=f"<b>{model_name}</b><br>r/r_s = %{{x:.3f}}<br>δ_B = %{{y:.5f}}<extra></extra>",
        ))

    # ── Hawking radiation δ_B as horizontal reference lines ──
    for name, cdata in sorted(hawking.items(), key=lambda x: x[1]["delta_b"]):
        wc = name.split("=")[1]
        color = HAWK_COLORS.get(wc, "#ffd700")
        db = cdata["delta_b"]
        fig.add_trace(go.Scatter(
            x=[10.0, -1.0], y=[db, db],
            mode="lines",
            name=f"Hawking ω_c={wc}  (δ_B={db:.4f})",
            line=dict(color=color, width=2.5, dash="dash"),
            hovertemplate=f"Hawking ω_c={wc}<br>δ_B = {db:.5f}<extra></extra>",
        ))

    # ── Landmark lines (event horizon, singularity) ──
    fig.add_vline(x=1.0, line=dict(color="#ff4444", width=3, dash="solid"),
                  annotation=dict(text="EVENT HORIZON", font=dict(size=16, color="#ff4444"),
                                  yref="paper", y=0.95))
    fig.add_vline(x=0.0, line=dict(color="#ffffff", width=2, dash="solid"),
                  annotation=dict(text="SINGULARITY", font=dict(size=14, color="#ffffff"),
                                  yref="paper", y=0.95, opacity=0.7))

    # Conformance zone
    fig.add_hrect(y0=0, y1=0.025, fillcolor="#7ee787", opacity=0.06,
                  line_width=0)
    fig.add_hline(y=0.025, line=dict(color="#7ee787", width=1, dash="dash"),
                  opacity=0.4)

    layout = dark_layout(
        "Black Hole — Hawking Radiation + Causal Set<br>"
        "<span style='font-size:14px;color:#8b949e'>"
        "Click legend entries to show/hide models  |  Double-click to isolate one</span>",
        "r / r_s  (far away  →  horizon  →  singularity  →  other side)",
        "δ_B  (Euclidean deviation from Benford)",
        height=750,
    )
    layout.xaxis.autorange = "reversed"
    fig.update_layout(layout)

    path = os.path.join(FIG_DIR, "10a_black_hole_interactive.html")
    fig.write_html(path, include_plotlyjs="cdn")
    print(f"  Saved {path}")


# ═══════════════════════════════════════════════════════════════════════════
# CHART 10b — BIG BANG (Planck Wall)
# ═══════════════════════════════════════════════════════════════════════════

def chart_big_bang():
    bb_data = load_json("results/round_trip/planck_wall_extended.json")
    wb_data = load_json("results/round_trip/whiteboard.json")

    # Hawking radiation from whiteboard
    hawking = {}
    for c in wb_data["candidates"]:
        if "Hawking" in c.get("name", "") and c["status"] == "EXISTS":
            hawking[c["name"]] = c

    fig = go.Figure()

    # ── All 10 QG models through the Planck Wall — all visible ──
    for model_name in bb_data["models_tested"]:
        entries = bb_data["results"][model_name]
        pts = [(e["T_planck_units"], e["delta_b"]) for e in entries if e.get("computable")]
        pts.sort()
        if not pts:
            continue
        ts, dbs = zip(*pts)

        fig.add_trace(go.Scatter(
            x=list(ts), y=list(dbs),
            mode="lines+markers",
            name=model_name,
            line=dict(color=MODEL_COLORS.get(model_name, MUTED), width=2.5),
            marker=dict(size=5),
            visible=True,
            hovertemplate=f"<b>{model_name}</b><br>T/T_P = %{{x:.3f}}<br>δ_B = %{{y:.5f}}<extra></extra>",
        ))

    # ── Hawking radiation δ_B as horizontal reference lines (hidden by default) ──
    for name, cdata in sorted(hawking.items(), key=lambda x: x[1]["delta_b"]):
        wc = name.split("=")[1]
        color = HAWK_COLORS.get(wc, "#ffd700")
        db = cdata["delta_b"]
        fig.add_trace(go.Scatter(
            x=[0.001, 100], y=[db, db],
            mode="lines",
            name=f"Hawking ω_c={wc}  (δ_B={db:.4f})",
            line=dict(color=color, width=2.5, dash="dash"),
            visible="legendonly",
            hovertemplate=f"Hawking ω_c={wc}<br>δ_B = {db:.5f}<extra></extra>",
        ))

    # ── Planck temperature landmark ──
    fig.add_vline(x=1.0, line=dict(color="#ff4444", width=3, dash="solid"),
                  annotation=dict(text="T = T_Planck", font=dict(size=16, color="#ff4444"),
                                  yref="paper", y=0.95))

    # Conformance zone
    fig.add_hrect(y0=0, y1=0.025, fillcolor="#7ee787", opacity=0.06,
                  line_width=0)
    fig.add_hline(y=0.025, line=dict(color="#7ee787", width=1, dash="dash"),
                  opacity=0.4)

    layout = dark_layout(
        "Big Bang — All 10 Quantum Gravity Models Through the Planck Wall<br>"
        "<span style='font-size:14px;color:#8b949e'>"
        "Click legend entries to show/hide models  |  Double-click to isolate one</span>",
        "T / T_Planck",
        "δ_B  (Euclidean deviation from Benford)",
        height=750,
    )
    layout.xaxis.type = "log"
    fig.update_layout(layout)

    path = os.path.join(FIG_DIR, "10b_big_bang_interactive.html")
    fig.write_html(path, include_plotlyjs="cdn")
    print(f"  Saved {path}")


# ═══════════════════════════════════════════════════════════════════════════
# CHART 10c — WORMHOLE (Morris-Thorne / Ellis)
# ═══════════════════════════════════════════════════════════════════════════

def chart_wormhole():
    wh_data = load_json("results/round_trip/wormhole_wall.json")
    wb_data = load_json("results/round_trip/whiteboard.json")

    b_0 = wh_data["wormhole"]["b_0"]

    # Hawking radiation from whiteboard
    hawking = {}
    for c in wb_data["candidates"]:
        if "Hawking" in c.get("name", "") and c["status"] == "EXISTS":
            hawking[c["name"]] = c

    fig = go.Figure()

    # ── All 10 QG models through the wormhole — all visible ──
    for model_name in wh_data["models_tested"]:
        entries = wh_data["traversing_observer"][model_name]
        pts = [(e["l_over_b0"], e["delta_b"]) for e in entries if e.get("computable")]
        pts.sort()
        if not pts:
            continue
        ls, dbs = zip(*pts)

        fig.add_trace(go.Scatter(
            x=list(ls), y=list(dbs),
            mode="lines+markers",
            name=model_name,
            line=dict(color=MODEL_COLORS.get(model_name, MUTED), width=2.5),
            marker=dict(size=5),
            visible=True,
            hovertemplate=f"<b>{model_name}</b><br>l/b₀ = %{{x:.1f}}<br>δ_B = %{{y:.5f}}<extra></extra>",
        ))

    # ── Pure Casimir trace ──
    pure_entries = wh_data.get("pure_casimir", [])
    pts = [(e["l_over_b0"], e["delta_b"]) for e in pure_entries if e.get("computable")]
    pts.sort()
    if pts:
        ls, dbs = zip(*pts)
        fig.add_trace(go.Scatter(
            x=list(ls), y=list(dbs),
            mode="lines+markers",
            name="Pure Casimir",
            line=dict(color="#ffffff", width=2.5, dash="dot"),
            marker=dict(size=5, symbol="diamond"),
            visible=True,
            hovertemplate="<b>Pure Casimir</b><br>l/b₀ = %{x:.1f}<br>δ_B = %{y:.5f}<extra></extra>",
        ))

    # ── Hawking radiation δ_B as horizontal reference lines (hidden by default) ──
    x_range = [-100, 100]
    for name, cdata in sorted(hawking.items(), key=lambda x: x[1]["delta_b"]):
        wc = name.split("=")[1]
        color = HAWK_COLORS.get(wc, "#ffd700")
        db = cdata["delta_b"]
        fig.add_trace(go.Scatter(
            x=x_range, y=[db, db],
            mode="lines",
            name=f"Hawking ω_c={wc}  (δ_B={db:.4f})",
            line=dict(color=color, width=2.5, dash="dash"),
            visible="legendonly",
            hovertemplate=f"Hawking ω_c={wc}<br>δ_B = {db:.5f}<extra></extra>",
        ))

    # ── Throat landmark ──
    fig.add_vline(x=0.0, line=dict(color="#ff4444", width=3, dash="solid"),
                  annotation=dict(text="THROAT", font=dict(size=16, color="#ff4444"),
                                  yref="paper", y=0.95))

    # Conformance zone
    fig.add_hrect(y0=0, y1=0.025, fillcolor="#7ee787", opacity=0.06,
                  line_width=0)
    fig.add_hline(y=0.025, line=dict(color="#7ee787", width=1, dash="dash"),
                  opacity=0.4)

    layout = dark_layout(
        "Wormhole — All 10 QG Models Through a Morris-Thorne Throat<br>"
        "<span style='font-size:14px;color:#8b949e'>"
        "Click legend entries to show/hide models  |  Double-click to isolate one</span>",
        "l / b₀  (proper distance from throat)",
        "δ_B  (Euclidean deviation from Benford)",
        height=750,
    )
    fig.update_layout(layout)

    path = os.path.join(FIG_DIR, "10c_wormhole_interactive.html")
    fig.write_html(path, include_plotlyjs="cdn")
    print(f"  Saved {path}")


# ═══════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    print("=" * 60)
    print("Generating interactive charts...")
    print("=" * 60)
    chart_black_hole()
    chart_big_bang()
    chart_wormhole()
    print("\nDone! Open the .html files in a browser.")
