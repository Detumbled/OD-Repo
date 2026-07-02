import csv
from html import escape
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
NO_VLBI_TRAJECTORY = SCRIPT_DIR / "voyager_od_trajectory_error.csv"
VLBI_TRAJECTORY = SCRIPT_DIR / "voyager_od_trajectory_error_VLBI.csv"
OUTPUT_SVG = SCRIPT_DIR / "voyager_vlbi_comparison.svg"


def read_trajectory(path):
    rows = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            rows.append(
                {
                    "hours": float(row["dt_since_reference_s"]) / 3600.0,
                    "position_error_km": float(row["position_error_norm_km"]),
                    "velocity_error_mm_s": float(row["velocity_error_norm_km_s"]) * 1.0e6,
                }
            )

    if not rows:
        raise RuntimeError(f"No trajectory rows found in {path}")
    return rows


def values(rows, key):
    return [row[key] for row in rows]


def scale(value, low, high, screen_low, screen_high):
    if high == low:
        return 0.5 * (screen_low + screen_high)
    return screen_low + (value - low) * (screen_high - screen_low) / (high - low)


def polyline(rows, x_key, y_key, x_domain, y_domain, panel, color):
    left, top, width, height = panel
    points = []
    for row in rows:
        x = scale(row[x_key], x_domain[0], x_domain[1], left, left + width)
        y = scale(row[y_key], y_domain[0], y_domain[1], top + height, top)
        points.append(f"{x:.2f},{y:.2f}")
    return (
        f'<polyline fill="none" stroke="{color}" stroke-width="2.5" '
        f'points="{" ".join(points)}" />'
    )


def axis_ticks(domain, count):
    low, high = domain
    if count <= 1:
        return [low]
    step = (high - low) / (count - 1)
    return [low + i * step for i in range(count)]


def draw_panel(title, y_label, y_key, no_vlbi, vlbi, panel, x_domain):
    left, top, width, height = panel
    all_y = values(no_vlbi, y_key) + values(vlbi, y_key)
    y_low = min(all_y)
    y_high = max(all_y)
    padding = 0.06 * (y_high - y_low if y_high != y_low else max(abs(y_high), 1.0))
    y_domain = (y_low - padding, y_high + padding)

    parts = [
        f'<text x="{left}" y="{top - 18}" font-size="17" font-weight="700">{escape(title)}</text>',
        f'<rect x="{left}" y="{top}" width="{width}" height="{height}" fill="#ffffff" stroke="#333333" />',
    ]

    for tick in axis_ticks(y_domain, 5):
        y = scale(tick, y_domain[0], y_domain[1], top + height, top)
        parts.append(f'<line x1="{left}" y1="{y:.2f}" x2="{left + width}" y2="{y:.2f}" stroke="#dddddd" />')
        parts.append(
            f'<text x="{left - 10}" y="{y + 4:.2f}" font-size="12" text-anchor="end">'
            f"{tick:.2f}</text>"
        )

    for tick in axis_ticks(x_domain, 7):
        x = scale(tick, x_domain[0], x_domain[1], left, left + width)
        parts.append(f'<line x1="{x:.2f}" y1="{top}" x2="{x:.2f}" y2="{top + height}" stroke="#eeeeee" />')
        parts.append(
            f'<text x="{x:.2f}" y="{top + height + 18}" font-size="12" text-anchor="middle">'
            f"{tick:.0f}</text>"
        )

    parts.append(polyline(no_vlbi, "hours", y_key, x_domain, y_domain, panel, "#1f77b4"))
    parts.append(polyline(vlbi, "hours", y_key, x_domain, y_domain, panel, "#d62728"))

    parts.append(
        f'<text x="{left - 62}" y="{top + height / 2:.2f}" font-size="13" '
        f'text-anchor="middle" transform="rotate(-90 {left - 62},{top + height / 2:.2f})">'
        f"{escape(y_label)}</text>"
    )
    return "\n".join(parts)


def write_svg(no_vlbi, vlbi):
    width = 1120
    height = 820
    panel_width = 900
    panel_height = 255
    left = 150
    panel_one = (left, 145, panel_width, panel_height)
    panel_two = (left, 505, panel_width, panel_height)

    x_domain = (
        min(values(no_vlbi, "hours") + values(vlbi, "hours")),
        max(values(no_vlbi, "hours") + values(vlbi, "hours")),
    )

    no_vlbi_final = no_vlbi[-1]["position_error_km"]
    vlbi_final = vlbi[-1]["position_error_km"]
    improvement = no_vlbi_final - vlbi_final

    parts = [
        '<svg xmlns="http://www.w3.org/2000/svg" width="1120" height="820" viewBox="0 0 1120 820">',
        '<rect width="1120" height="820" fill="#f8f9fb" />',
        '<style>text { font-family: Arial, Helvetica, sans-serif; fill: #202124; }</style>',
        '<text x="560" y="42" text-anchor="middle" font-size="24" font-weight="700">'
        "Voyager OD Trajectory Error: VLBI Comparison</text>",
        '<text x="560" y="72" text-anchor="middle" font-size="15">'
        f"Final position error: {no_vlbi_final:.3f} km without VLBI, "
        f"{vlbi_final:.3f} km with VLBI ({improvement:.3f} km lower)</text>",
        draw_panel(
            "Position Error",
            "Position error norm (km)",
            "position_error_km",
            no_vlbi,
            vlbi,
            panel_one,
            x_domain,
        ),
        draw_panel(
            "Velocity Error",
            "Velocity error norm (mm/s)",
            "velocity_error_mm_s",
            no_vlbi,
            vlbi,
            panel_two,
            x_domain,
        ),
        '<line x1="720" y1="100" x2="750" y2="100" stroke="#1f77b4" stroke-width="3" />',
        '<text x="760" y="105" font-size="14">Range + range-rate</text>',
        '<line x1="720" y1="124" x2="750" y2="124" stroke="#d62728" stroke-width="3" />',
        '<text x="760" y="129" font-size="14">Range + range-rate + VLBI</text>',
        '<text x="600" y="807" text-anchor="middle" font-size="13">Hours since first observation</text>',
        "</svg>",
    ]
    OUTPUT_SVG.write_text("\n".join(parts) + "\n")
    print(f"Wrote {OUTPUT_SVG}")


def main():
    no_vlbi = read_trajectory(NO_VLBI_TRAJECTORY)
    vlbi = read_trajectory(VLBI_TRAJECTORY)
    write_svg(no_vlbi, vlbi)


if __name__ == "__main__":
    main()
