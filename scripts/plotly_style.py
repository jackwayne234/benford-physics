"""Shared dark-theme Plotly style for all interactive figures."""

DARK_THEME = dict(
    paper_bgcolor="#0d1117",
    plot_bgcolor="#161b22",
    font=dict(family="monospace", size=14, color="#c9d1d9"),
    margin=dict(l=90, r=60, t=120, b=90),
)

# Legend style applied LAST so nothing can override it
LEGEND_STYLE = dict(
    font=dict(size=14, color="#000000"),
    bgcolor="rgba(255,255,255,0.95)",
    bordercolor="#888888",
    borderwidth=2,
    itemclick="toggle",
    itemdoubleclick="toggleothers",
)

AXIS_STYLE = dict(
    gridcolor="#21262d",
    zeroline=False,
    tickfont=dict(color="#8b949e", size=13),
    title_font=dict(size=15),
)

# Consistent model colors (matches existing interactive charts)
MODEL_COLORS = {
    "Standard":      "#58a6ff",
    "LQG":           "#7ee787",
    "GUP":           "#f78166",
    "DSR":           "#d2a8ff",
    "Hagedorn":      "#ffd700",
    "Causal Set":    "#00d4ff",
    "Asym. Safety":  "#ff69b4",
    "Horava-Lif.":   "#ff6347",
    "Noncommut.":    "#98fb98",
    "CDT":           "#dda0dd",
}

# Colors for the 4 fundamental distributions
DIST_COLORS = {
    "Bose-Einstein":    "#58a6ff",
    "Fermi-Dirac":      "#f78166",
    "Maxwell-Boltzmann": "#7ee787",
    "Planck":           "#d2a8ff",
}

CONFORM_THRESHOLD = 0.025  # green zone upper bound
DEVIATE_THRESHOLD = 0.1


def apply_dark_theme(fig, **overrides):
    """Apply the dark theme to a plotly figure."""
    layout = {**DARK_THEME, **overrides}
    fig.update_layout(**layout)
    fig.update_xaxes(**AXIS_STYLE)
    fig.update_yaxes(**AXIS_STYLE)
    # Force legend style LAST — black text on white background, always readable
    fig.update_layout(legend=LEGEND_STYLE)
    return fig


def add_conform_zone(fig, row=None, col=None):
    """Add green conformance zone shading."""
    kwargs = {}
    if row is not None:
        kwargs["row"] = row
        kwargs["col"] = col
    fig.add_hrect(
        y0=0, y1=CONFORM_THRESHOLD,
        fillcolor="#7ee787", opacity=0.06,
        line_width=0, **kwargs
    )
    fig.add_hline(
        y=CONFORM_THRESHOLD,
        line_dash="dash", line_color="#7ee787", line_width=1,
        opacity=0.4, **kwargs
    )


def save_figure(fig, path):
    """Save figure as standalone HTML."""
    fig.write_html(path, include_plotlyjs="cdn")
    print(f"  Saved: {path}")
