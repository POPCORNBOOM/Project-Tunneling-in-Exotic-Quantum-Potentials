"""External field helpers and barrier geometry."""
from __future__ import annotations

from typing import Tuple, Optional
import numpy as np

Array = np.ndarray


def apply_field(Vx: Array, x: Array, F: float) -> Array:
    return Vx - F * x


def barrier_top(x: Array, Vx: Array, x_min: float = float('-inf'), x_max: float = float('inf')) -> int:
    mask = (x >= x_min) & (x <= x_max)
    x = x[mask]
    Vx = Vx[mask]
    return np.argmax(Vx)
