"""External field helpers and barrier geometry."""
from __future__ import annotations

from typing import Tuple, Optional
import numpy as np

Array = np.ndarray


def apply_field(Vx: Array, x: Array, F: float) -> Array:
    return Vx - F * x


def barrier_top(x: Array, Vx: Array) -> Tuple[float, float]:
    idx = np.argmax(Vx)
    return float(x[idx]), float(Vx[idx])


def classify_barrier(E: float, Vx: Array,) -> str:
    vmax = float(np.max(Vx))
    if E >= vmax:
        return "over_the_barrier"
    return "not_over_the_barrier"
