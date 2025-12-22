"""External field helpers and barrier geometry."""
from __future__ import annotations

from typing import Tuple, Optional
import numpy as np

Array = np.ndarray


def apply_field(Vx: Array, x: Array, F: float) -> Array:
    return Vx - F * x


def barrier_top(x: Array, Vx: Array, xmin: None | float = None, xmax: None | float = None) -> Tuple[float, float]:
    if xmin is not None:
        mask = x >= xmin
        x = x[mask]
        Vx = Vx[mask]
    if xmax is not None:
        mask = x <= xmax
        x = x[mask]
        Vx = Vx[mask]
    idx = np.argmax(Vx)
    return float(x[idx]), float(Vx[idx])


def classify_barrier(E: float, Vx: Array,) -> str:
    vmax = float(np.max(Vx))
    if E >= vmax:
        return "over_the_barrier"
    return "not_over_the_barrier"
