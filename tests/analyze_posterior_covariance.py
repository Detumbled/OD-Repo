#!/usr/bin/env python3
"""Analyze and plot the 3-sigma position covariance ellipsoid from WLS CSV output."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_COVARIANCE_PATH = SCRIPT_DIR / "posterior_covariance.csv"


def read_posterior_covariance(path: Path) -> np.ndarray:
    """Read and validate the header-free 6x6 posterior covariance CSV."""
    try:
        covariance = np.loadtxt(path, delimiter=",", dtype=np.float64)
    except OSError as error:
        raise RuntimeError(f"Cannot read posterior covariance CSV: {path}") from error
    except ValueError as error:
        raise RuntimeError(f"Posterior covariance CSV is not a numeric matrix: {path}") from error

    if covariance.shape != (6, 6):
        raise ValueError(
            f"Expected a 6x6 posterior covariance matrix, received shape {covariance.shape}."
        )
    if not np.all(np.isfinite(covariance)):
        raise ValueError("Posterior covariance matrix contains non-finite values.")

    return covariance


def position_eigensystem(covariance: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return ascending position variances and their unit principal directions."""
    position_covariance = covariance[:3, :3]
    if not np.allclose(position_covariance, position_covariance.T, rtol=1.0e-10, atol=1.0e-12):
        raise ValueError("Position covariance submatrix must be symmetric.")

    eigenvalues, eigenvectors = np.linalg.eigh(position_covariance)
    scale = max(float(np.max(np.abs(eigenvalues))), 1.0)
    negative_tolerance = 100.0 * np.finfo(np.float64).eps * scale
    if np.any(eigenvalues < -negative_tolerance):
        raise ValueError(
            "Position covariance has materially negative eigenvalues; it is not positive semidefinite."
        )

    return np.maximum(eigenvalues, 0.0), eigenvectors


def ellipsoid_surface(
    semi_axes_km: np.ndarray,
    principal_directions: np.ndarray,
    latitude_samples: int,
    longitude_samples: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Vectorize the rotated ellipsoid surface construction without per-point loops."""
    azimuth, elevation = np.meshgrid(
        np.linspace(0.0, 2.0 * np.pi, longitude_samples, endpoint=True),
        np.linspace(-0.5 * np.pi, 0.5 * np.pi, latitude_samples, endpoint=True),
    )
    unit_sphere = np.stack(
        (
            np.cos(elevation) * np.cos(azimuth),
            np.cos(elevation) * np.sin(azimuth),
            np.sin(elevation),
        ),
        axis=0,
    ).reshape(3, -1)
    surface = principal_directions @ (semi_axes_km[:, np.newaxis] * unit_sphere)
    return tuple(component.reshape(elevation.shape) for component in surface)


def create_figure(
    semi_axes_3sigma_km: np.ndarray,
    principal_directions: np.ndarray,
    latitude_samples: int,
    longitude_samples: int,
):
    """Create a rotatable Matplotlib 3D view in position-error coordinates."""
    try:
        import matplotlib.pyplot as plt
    except ModuleNotFoundError as error:
        raise RuntimeError(
            "Matplotlib is required for the interactive 3D plot. "
            "Install it with 'python -m pip install matplotlib'."
        ) from error

    x, y, z = ellipsoid_surface(
        semi_axes_3sigma_km,
        principal_directions,
        latitude_samples,
        longitude_samples,
    )
    figure = plt.figure(figsize=(9, 8))
    axis = figure.add_subplot(projection="3d")
    axis.plot_surface(x, y, z, cmap="Blues", alpha=0.65, linewidth=0, antialiased=True)

    axis_colors = ("#d62728", "#2ca02c", "#1f77b4")
    for index, (semi_axis, color) in enumerate(zip(semi_axes_3sigma_km, axis_colors, strict=True)):
        endpoint = semi_axis * principal_directions[:, index]
        axis.plot(
            (-endpoint[0], endpoint[0]),
            (-endpoint[1], endpoint[1]),
            (-endpoint[2], endpoint[2]),
            color=color,
            linewidth=2.5,
            label=f"Principal axis {index + 1}: {semi_axis:.3g} km",
        )

    maximum_extent = max(float(np.max(semi_axes_3sigma_km)), 1.0)
    axis.set_xlim(-maximum_extent, maximum_extent)
    axis.set_ylim(-maximum_extent, maximum_extent)
    axis.set_zlim(-maximum_extent, maximum_extent)
    axis.set_box_aspect((1.0, 1.0, 1.0))
    axis.set_xlabel("Position error X (km)")
    axis.set_ylabel("Position error Y (km)")
    axis.set_zlabel("Position error Z (km)")
    axis.set_title("3-sigma Position-Error Covariance Ellipsoid")
    axis.legend(loc="upper left")
    figure.tight_layout()
    return figure, plt


def print_results(
    eigenvalues_km2: np.ndarray,
    semi_axes_1sigma_km: np.ndarray,
    semi_axes_3sigma_km: np.ndarray,
    principal_directions: np.ndarray,
) -> None:
    print("Position covariance eigensystem (ascending eigenvalue):")
    for index in range(3):
        direction = principal_directions[:, index]
        print(
            f"  Axis {index + 1}: variance = {eigenvalues_km2[index]:.12e} km^2, "
            f"1-sigma = {semi_axes_1sigma_km[index]:.12e} km, "
            f"3-sigma = {semi_axes_3sigma_km[index]:.12e} km"
        )
        print(f"           direction = [{direction[0]: .12e}, {direction[1]: .12e}, {direction[2]: .12e}]")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Analyze a WLS 6x6 posterior covariance CSV and display its 3-sigma position ellipsoid."
    )
    parser.add_argument(
        "covariance_csv",
        nargs="?",
        type=Path,
        default=DEFAULT_COVARIANCE_PATH,
        help=f"Header-free 6x6 covariance CSV (default: {DEFAULT_COVARIANCE_PATH})",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="Optional path for a static image of the interactive plot.",
    )
    parser.add_argument(
        "--no-show",
        action="store_true",
        help="Build the plot without opening it; useful together with --output.",
    )
    parser.add_argument(
        "--latitude-samples",
        type=int,
        default=48,
        help="Latitude samples for the rendered surface (default: 48).",
    )
    parser.add_argument(
        "--longitude-samples",
        type=int,
        default=96,
        help="Longitude samples for the rendered surface (default: 96).",
    )
    arguments = parser.parse_args()
    if arguments.latitude_samples < 3 or arguments.longitude_samples < 3:
        parser.error("Surface sample counts must both be at least 3.")
    if arguments.no_show and arguments.output is None:
        parser.error("--no-show requires --output so the plot is still generated.")
    return arguments


def main() -> None:
    arguments = parse_args()
    covariance = read_posterior_covariance(arguments.covariance_csv)
    eigenvalues_km2, principal_directions = position_eigensystem(covariance)
    semi_axes_1sigma_km = np.sqrt(eigenvalues_km2)
    semi_axes_3sigma_km = 3.0 * semi_axes_1sigma_km
    print_results(
        eigenvalues_km2,
        semi_axes_1sigma_km,
        semi_axes_3sigma_km,
        principal_directions,
    )

    figure, pyplot = create_figure(
        semi_axes_3sigma_km,
        principal_directions,
        arguments.latitude_samples,
        arguments.longitude_samples,
    )
    if arguments.output is not None:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(arguments.output, dpi=180, bbox_inches="tight")
        print(f"Wrote ellipsoid plot image: {arguments.output}")
    if not arguments.no_show:
        pyplot.show()


if __name__ == "__main__":
    main()
