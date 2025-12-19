"""Turning points and WKB action integrals."""
from __future__ import annotations

from typing import Optional
import numpy as np

Array = np.ndarray


# Helper function to find turning points
def find_turning_points(
    x: Array,
    Vx: Array,
    E: float,
    x_min: Optional[float] = None,
    x_max: Optional[float] = None,
    atol: float = 1e-12,
    merge_tol: Optional[float] = None,
) -> Array:
    x = np.asarray(x, dtype=float)
    Vx = np.asarray(Vx, dtype=float)
    if x.ndim != 1 or Vx.ndim != 1 or x.size != Vx.size:
        raise ValueError("x and Vx must be 1D arrays of the same length.")

    if x[0] > x[-1]:
        x = x[::-1]
        Vx = Vx[::-1]

    mask = np.ones_like(x, dtype=bool)
    if x_min is not None:
        mask &= x >= x_min
    if x_max is not None:
        mask &= x <= x_max

    idx = np.where(mask)[0]
    if idx.size < 2:
        return np.array([], dtype=float)

    xw = x[idx]
    fw = Vx[idx] - E

    if merge_tol is None:
        dx = np.median(np.diff(xw))
        merge_tol = 1.5 * dx

    roots = []

    near = np.where(np.abs(fw) <= atol)[0]
    if near.size > 0:
        start = near[0]
        prev = near[0]
        for k in near[1:]:
            if k == prev + 1:
                prev = k
            else:
                roots.append(0.5 * (xw[start] + xw[prev]))
                start = prev = k
        roots.append(0.5 * (xw[start] + xw[prev]))

    good = (np.abs(fw[:-1]) > atol) & (np.abs(fw[1:]) > atol)
    sign_change = good & (fw[:-1] * fw[1:] < 0)

    for i in np.where(sign_change)[0]:
        x0, x1 = xw[i], xw[i + 1]
        f0, f1 = fw[i], fw[i + 1]
        xr = x0 - f0 * (x1 - x0) / (f1 - f0)
        roots.append(xr)

    if not roots:
        return np.array([], dtype=float)

    roots = np.array(sorted(roots), dtype=float)

    merged = [roots[0]]
    for r in roots[1:]:
        if abs(r - merged[-1]) > merge_tol:
            merged.append(r)
        else:
            merged[-1] = 0.5 * (merged[-1] + r)
    return np.array(merged, dtype=float)

# Helper function to calculate action integral
def action_integral(x: Array, Vx: Array, E: float, m: float = 1.0) -> float:
    mask = Vx > E
    xb = x[mask]
    barrier = Vx[mask] - E
    integrand = np.sqrt(np.maximum(0.0, 2.0 * m * barrier))
    return float(np.trapz(integrand, xb))


def wkb_transmission(S: float, hbar: float = 1.0) -> float:
    return float(np.exp(-2.0 * S / hbar))
