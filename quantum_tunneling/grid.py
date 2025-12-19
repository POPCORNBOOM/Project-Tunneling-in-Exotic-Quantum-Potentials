"""Grid and domain helpers."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Tuple
import numpy as np

Array = np.ndarray


@dataclass(frozen=True)
class GridSpec:
    L: float
    N: int

    @property
    def dx(self) -> float:
        return 2.0 * self.L / (self.N + 1)

    def build_grid(self) -> Tuple[Array, float]:
        dx = self.dx
        x = np.linspace(-self.L + dx, self.L - dx, self.N)
        return x, dx


def finite_difference_laplacian(N: int, dx: float) -> Tuple[Array, Array]:
    pref = 1.0 / (dx * dx)
    main = 2.0 * pref * np.ones(N)
    off = -pref * np.ones(N - 1)
    return main, off
