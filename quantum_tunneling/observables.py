"""Diagnostics and observables for eigenstates and dynamics."""
from __future__ import annotations

from typing import Tuple
import numpy as np

Array = np.ndarray


def localization_metrics(x: Array, psi: Array, dx: float) -> Tuple[float, float, float]:
    rho = np.abs(psi) ** 2
    mean_x = float(np.sum(x * rho) * dx)
    mean_x2 = float(np.sum((x * x) * rho) * dx)
    sigma = float(np.sqrt(max(mean_x2 - mean_x * mean_x, 0.0)))
    ipr = float(np.sum(rho * rho) * dx)
    return mean_x, sigma, ipr


def forbidden_probability(Vx: Array, En: float, psi: Array, dx: float) -> float:
    forb = Vx > En
    return float(np.sum(np.abs(psi[forb]) ** 2) * dx)


# probability current or flux
def probability_current(psi: Array, dx: float, hbar: float = 1.0, m: float = 1.0) -> Array:
    """Return local probability current j(x) = (hbar/m) Im(psi* dpsi/dx)."""
    dpsi_dx = np.gradient(psi, dx)
    return (hbar / m) * np.imag(np.conj(psi) * dpsi_dx)


def compute_probability(psi: Array, dx: float, mask: Array | None = None) -> float:
    if mask is None:
        return float(np.sum(np.abs(psi) ** 2) * dx)
    return float(np.sum(np.abs(psi[mask]) ** 2) * dx)
