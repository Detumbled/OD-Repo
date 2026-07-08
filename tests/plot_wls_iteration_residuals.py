#!/usr/bin/env python3
"""Plot Voyager WLS residuals and state-error history by iteration."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from html import escape
from pathlib import Path
from typing import Optional

try:
    import matplotlib.pyplot as plt
except ModuleNotFoundError:
    plt = None


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_RESIDUALS = SCRIPT_DIR / "voyager_wls_iteration_residuals_VLBI.csv"
DEFAULT_SUMMARY = SCRIPT_DIR / "voyager_wls_iteration_summary_VLBI.csv"

MEASUREMENT_PLOTS = {
    "range": {
        "title": "Range Residuals by WLS Iteration",
        "ylabel": "Range residual (m)",
        "filename": "wls_range_residuals.png",
    },
    "range_rate": {
        "title": "Range-Rate Residuals by WLS Iteration",
        "ylabel": "Range-rate residual (mm/s)",
        "filename": "wls_range_rate_residuals.png",
    },
    "vlbi": {
        "title": "VLBI Differential-Range Residuals by WLS Iteration",
        "ylabel": "VLBI residual (m)",
        "filename": "wls_vlbi_residuals.png",
    },
}


def read_residuals(path: Path) -> list[dict]:
    rows = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            rows.append(
                {
                    "iteration": int(row["iteration"]),
                    "measurement_type": row["measurement_type"],
                    "station_label": row["station_label"],
                    "hours": float(row["dt_since_reference_s"]) / 3600.0,
                    "residual": float(row["residual_plot_value"]),
                    "sigma": float(row["sigma_plot_value"]),
                    "normalized_residual": float(row["normalized_residual"]),
                }
            )

    if not rows:
        raise RuntimeError(f"No residual rows found in {path}")
    return rows


def read_summary(path: Path) -> list[dict]:
    rows = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            rows.append(
                {
                    "iteration": int(row["iteration"]),
                    "position_error_rms_km": float(row["position_error_rms_km"]),
                    "velocity_error_rms_mm_s": float(row["velocity_error_rms_mm_s"]),
                    "position_error_norm_km": float(row["position_error_norm_km"]),
                    "velocity_error_norm_mm_s": float(row["velocity_error_norm_km_s"]) * 1.0e6,
                }
            )

    if not rows:
        raise RuntimeError(f"No iteration-summary rows found in {path}")
    return rows


def set_window_title(fig, title: str) -> None:
    manager = getattr(fig.canvas, "manager", None)
    if manager is not None and hasattr(manager, "set_window_title"):
        manager.set_window_title(title)


def save_if_requested(fig, output_dir: Optional[Path], filename: str) -> None:
    if output_dir is None:
        return
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / filename
    fig.savefig(output_path, dpi=160, bbox_inches="tight")
    print(f"Wrote {output_path}")


def plot_measurement_residuals(
    residual_rows: list[dict],
    measurement_type: str,
    output_dir: Optional[Path],
    normalized: bool,
) -> None:
    config = MEASUREMENT_PLOTS[measurement_type]
    rows = [row for row in residual_rows if row["measurement_type"] == measurement_type]
    if not rows:
        print(f"Skipping {measurement_type}: no rows")
        return

    rows_by_iteration: dict[int, list[dict]] = defaultdict(list)
    for row in rows:
        rows_by_iteration[row["iteration"]].append(row)

    fig, ax = plt.subplots(figsize=(13, 6))
    set_window_title(fig, config["title"])

    value_key = "normalized_residual" if normalized else "residual"
    ylabel = "Normalized residual (sigma)" if normalized else config["ylabel"]

    for iteration in sorted(rows_by_iteration):
        iteration_rows = rows_by_iteration[iteration]
        ax.scatter(
            [row["hours"] for row in iteration_rows],
            [row[value_key] for row in iteration_rows],
            s=7,
            alpha=0.45,
            label=f"iter {iteration}",
        )

    ax.axhline(0.0, color="black", linewidth=0.9, alpha=0.7)
    ax.set_title(config["title"])
    ax.set_xlabel("Hours since reference epoch")
    ax.set_ylabel(ylabel)
    ax.grid(True, linestyle="--", alpha=0.35)
    ax.legend(title="WLS iteration", ncols=2, fontsize=9)
    fig.tight_layout()

    stem = Path(config["filename"]).stem
    filename = f"{stem}_normalized.png" if normalized else config["filename"]
    save_if_requested(fig, output_dir, filename)


def plot_state_error(summary_rows: list[dict], output_dir: Optional[Path], use_norm: bool) -> None:
    iterations = [row["iteration"] for row in summary_rows]
    if use_norm:
        position_values = [row["position_error_norm_km"] for row in summary_rows]
        velocity_values = [row["velocity_error_norm_mm_s"] for row in summary_rows]
        position_label = "Position error norm (km)"
        velocity_label = "Velocity error norm (mm/s)"
        title_suffix = "Norm"
        filename = "wls_state_error_norms.png"
    else:
        position_values = [row["position_error_rms_km"] for row in summary_rows]
        velocity_values = [row["velocity_error_rms_mm_s"] for row in summary_rows]
        position_label = "Position RMS component error (km)"
        velocity_label = "Velocity RMS component error (mm/s)"
        title_suffix = "Component RMS"
        filename = "wls_state_error_rms.png"

    fig, axes = plt.subplots(2, 1, figsize=(9, 7), sharex=True)
    set_window_title(fig, f"WLS State Error {title_suffix}")

    axes[0].plot(iterations, position_values, marker="o", linewidth=2)
    axes[0].set_ylabel(position_label)
    axes[0].set_title(f"State Error by WLS Iteration ({title_suffix})")
    axes[0].grid(True, linestyle="--", alpha=0.35)

    axes[1].plot(iterations, velocity_values, marker="o", linewidth=2, color="tab:orange")
    axes[1].set_xlabel("WLS iteration")
    axes[1].set_ylabel(velocity_label)
    axes[1].grid(True, linestyle="--", alpha=0.35)
    axes[1].set_xticks(iterations)

    fig.tight_layout()
    save_if_requested(fig, output_dir, filename)


def axis_ticks(low: float, high: float, count: int) -> list[float]:
    if count <= 1 or high == low:
        return [low]
    step = (high - low) / float(count - 1)
    return [low + step * index for index in range(count)]


def padded_domain(values: list[float]) -> tuple[float, float]:
    low = min(values)
    high = max(values)
    span = high - low
    if span == 0.0:
        padding = max(abs(low), 1.0) * 0.1
    else:
        padding = 0.06 * span
    return low - padding, high + padding


def scale(value: float, low: float, high: float, screen_low: float, screen_high: float) -> float:
    if high == low:
        return 0.5 * (screen_low + screen_high)
    return screen_low + (value - low) * (screen_high - screen_low) / (high - low)


def svg_header(width: int, height: int, title: str) -> list[str]:
    return [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">',
        f'<rect width="{width}" height="{height}" fill="#fbfbfd" />',
        '<style>text { font-family: Arial, Helvetica, sans-serif; fill: #202124; } '
        '.axis { stroke: #333333; stroke-width: 1; } '
        '.grid { stroke: #dddddd; stroke-width: 1; }</style>',
        f'<text x="{width / 2:.1f}" y="34" text-anchor="middle" '
        f'font-size="22" font-weight="700">{escape(title)}</text>',
    ]


def write_svg_residual_plot(
    residual_rows: list[dict],
    measurement_type: str,
    output_dir: Path,
    normalized: bool,
) -> None:
    config = MEASUREMENT_PLOTS[measurement_type]
    rows = [row for row in residual_rows if row["measurement_type"] == measurement_type]
    if not rows:
        return

    value_key = "normalized_residual" if normalized else "residual"
    ylabel = "Normalized residual (sigma)" if normalized else config["ylabel"]
    title = config["title"]
    filename = Path(config["filename"]).with_suffix(".svg").name
    if normalized:
        filename = f"{Path(filename).stem}_normalized.svg"

    width = 1200
    height = 650
    left = 95
    top = 70
    plot_width = 920
    plot_height = 500
    colors = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#ff7f0e", "#17becf", "#8c564b"]
    x_values = [row["hours"] for row in rows]
    y_values = [row[value_key] for row in rows]
    x_low, x_high = padded_domain(x_values)
    y_low, y_high = padded_domain(y_values)

    parts = svg_header(width, height, title)
    parts.append(
        f'<rect x="{left}" y="{top}" width="{plot_width}" height="{plot_height}" '
        'fill="#ffffff" stroke="#333333" />'
    )

    for tick in axis_ticks(x_low, x_high, 7):
        x = scale(tick, x_low, x_high, left, left + plot_width)
        parts.append(f'<line class="grid" x1="{x:.2f}" y1="{top}" x2="{x:.2f}" y2="{top + plot_height}" />')
        parts.append(
            f'<text x="{x:.2f}" y="{top + plot_height + 22}" font-size="12" '
            f'text-anchor="middle">{tick:.0f}</text>'
        )

    for tick in axis_ticks(y_low, y_high, 7):
        y = scale(tick, y_low, y_high, top + plot_height, top)
        parts.append(f'<line class="grid" x1="{left}" y1="{y:.2f}" x2="{left + plot_width}" y2="{y:.2f}" />')
        parts.append(
            f'<text x="{left - 10}" y="{y + 4:.2f}" font-size="12" '
            f'text-anchor="end">{tick:.3g}</text>'
        )

    zero_y = scale(0.0, y_low, y_high, top + plot_height, top)
    parts.append(
        f'<line x1="{left}" y1="{zero_y:.2f}" x2="{left + plot_width}" y2="{zero_y:.2f}" '
        'stroke="#111111" stroke-width="1.2" />'
    )

    rows_by_iteration: dict[int, list[dict]] = defaultdict(list)
    for row in rows:
        rows_by_iteration[row["iteration"]].append(row)

    legend_x = left + plot_width + 35
    legend_y = top + 20
    for color_index, iteration in enumerate(sorted(rows_by_iteration)):
        color = colors[color_index % len(colors)]
        for row in rows_by_iteration[iteration]:
            x = scale(row["hours"], x_low, x_high, left, left + plot_width)
            y = scale(row[value_key], y_low, y_high, top + plot_height, top)
            parts.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="1.4" fill="{color}" fill-opacity="0.45" />')
        parts.append(
            f'<circle cx="{legend_x}" cy="{legend_y + 22 * color_index}" r="5" fill="{color}" />'
        )
        parts.append(
            f'<text x="{legend_x + 12}" y="{legend_y + 22 * color_index + 4}" '
            f'font-size="13">iter {iteration}</text>'
        )

    parts.append(
        f'<text x="{left + plot_width / 2:.2f}" y="{height - 26}" '
        'font-size="14" text-anchor="middle">Hours since reference epoch</text>'
    )
    parts.append(
        f'<text x="24" y="{top + plot_height / 2:.2f}" font-size="14" text-anchor="middle" '
        f'transform="rotate(-90 24,{top + plot_height / 2:.2f})">{escape(ylabel)}</text>'
    )
    parts.append("</svg>")

    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / filename
    output_path.write_text("\n".join(parts) + "\n")
    print(f"Wrote {output_path}")


def write_svg_state_error(summary_rows: list[dict], output_dir: Path, use_norm: bool) -> None:
    if use_norm:
        position_key = "position_error_norm_km"
        velocity_key = "velocity_error_norm_mm_s"
        title = "WLS State Error Norms"
        position_label = "Position error norm (km)"
        velocity_label = "Velocity error norm (mm/s)"
        filename = "wls_state_error_norms.svg"
    else:
        position_key = "position_error_rms_km"
        velocity_key = "velocity_error_rms_mm_s"
        title = "WLS State Error Component RMS"
        position_label = "Position RMS component error (km)"
        velocity_label = "Velocity RMS component error (mm/s)"
        filename = "wls_state_error_rms.svg"

    width = 980
    height = 720
    left = 110
    plot_width = 760
    panel_height = 240
    panel_tops = [90, 415]
    iterations = [row["iteration"] for row in summary_rows]
    x_low, x_high = padded_domain([float(value) for value in iterations])
    colors = ["#1f77b4", "#ff7f0e"]

    parts = svg_header(width, height, title)
    for panel_index, (top, key, ylabel, color) in enumerate(
        zip(panel_tops, [position_key, velocity_key], [position_label, velocity_label], colors)
    ):
        y_values = [row[key] for row in summary_rows]
        y_low, y_high = padded_domain(y_values)
        parts.append(
            f'<rect x="{left}" y="{top}" width="{plot_width}" height="{panel_height}" '
            'fill="#ffffff" stroke="#333333" />'
        )
        for tick in axis_ticks(y_low, y_high, 5):
            y = scale(tick, y_low, y_high, top + panel_height, top)
            parts.append(f'<line class="grid" x1="{left}" y1="{y:.2f}" x2="{left + plot_width}" y2="{y:.2f}" />')
            parts.append(
                f'<text x="{left - 10}" y="{y + 4:.2f}" font-size="12" '
                f'text-anchor="end">{tick:.3g}</text>'
            )
        for tick in iterations:
            x = scale(float(tick), x_low, x_high, left, left + plot_width)
            parts.append(f'<line class="grid" x1="{x:.2f}" y1="{top}" x2="{x:.2f}" y2="{top + panel_height}" />')
            if panel_index == 1:
                parts.append(
                    f'<text x="{x:.2f}" y="{top + panel_height + 22}" font-size="12" '
                    f'text-anchor="middle">{tick}</text>'
                )

        points = []
        for row in summary_rows:
            x = scale(float(row["iteration"]), x_low, x_high, left, left + plot_width)
            y = scale(row[key], y_low, y_high, top + panel_height, top)
            points.append(f"{x:.2f},{y:.2f}")
            parts.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="4" fill="{color}" />')
        parts.append(
            f'<polyline fill="none" stroke="{color}" stroke-width="2.2" points="{" ".join(points)}" />'
        )
        parts.append(
            f'<text x="28" y="{top + panel_height / 2:.2f}" font-size="14" text-anchor="middle" '
            f'transform="rotate(-90 28,{top + panel_height / 2:.2f})">{escape(ylabel)}</text>'
        )

    parts.append(
        f'<text x="{left + plot_width / 2:.2f}" y="{height - 24}" '
        'font-size="14" text-anchor="middle">WLS iteration</text>'
    )
    parts.append("</svg>")

    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / filename
    output_path.write_text("\n".join(parts) + "\n")
    print(f"Wrote {output_path}")


def write_svg_fallback(
    residuals: list[dict],
    summary: list[dict],
    output_dir: Optional[Path],
    normalized: bool,
    state_error_norm: bool,
) -> None:
    fallback_dir = output_dir or (SCRIPT_DIR / "wls_iteration_plots")
    print("Matplotlib is not installed; writing SVG plots instead of opening interactive windows.")
    for measurement_type in ("range", "range_rate", "vlbi"):
        write_svg_residual_plot(residuals, measurement_type, fallback_dir, normalized)
    write_svg_state_error(summary, fallback_dir, state_error_norm)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot Voyager WLS residuals for range, range-rate, VLBI, and state-error history."
    )
    parser.add_argument(
        "--residuals",
        type=Path,
        default=DEFAULT_RESIDUALS,
        help=f"Per-iteration residual CSV, default: {DEFAULT_RESIDUALS}",
    )
    parser.add_argument(
        "--summary",
        type=Path,
        default=DEFAULT_SUMMARY,
        help=f"Per-iteration summary CSV, default: {DEFAULT_SUMMARY}",
    )
    parser.add_argument(
        "--save-dir",
        type=Path,
        default=None,
        help="Optional directory for PNG outputs.",
    )
    parser.add_argument(
        "--no-show",
        action="store_true",
        help="Do not open Matplotlib windows. Useful with --save-dir.",
    )
    parser.add_argument(
        "--normalized",
        action="store_true",
        help="Plot measurement residuals divided by their raw observation sigma.",
    )
    parser.add_argument(
        "--state-error-norm",
        action="store_true",
        help="Plot 3D position/velocity error norms instead of component RMS errors.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    residuals = read_residuals(args.residuals)
    summary = read_summary(args.summary)

    if plt is None:
        write_svg_fallback(
            residuals,
            summary,
            args.save_dir,
            args.normalized,
            args.state_error_norm,
        )
        return

    for measurement_type in ("range", "range_rate", "vlbi"):
        plot_measurement_residuals(
            residuals,
            measurement_type,
            args.save_dir,
            args.normalized,
        )
    plot_state_error(summary, args.save_dir, args.state_error_norm)

    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    main()
