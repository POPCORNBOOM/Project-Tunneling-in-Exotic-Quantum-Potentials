"""Potential definitions and factories for exotic 1D models."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Dict, Tuple
import numpy as np


Array = np.ndarray


@dataclass(frozen=True)
class PotentialSpec:
    kind: str
    params: Dict[str, float]

    def build(self) -> Callable[[Array], Array]:
        fn = POTENTIAL_REGISTRY.get(self.kind)
        if fn is None:
            raise ValueError(f"Unknown potential kind: {self.kind}")
        return lambda x: fn(x, **self.params)


def V_cusp(x: Array, V0: float = 1.0, alpha: float = 0.5) -> Array:
    return V0 * np.power(np.abs(x), alpha)


def V_exponential(x: Array, V0: float = 1.0, a: float = 1.0) -> Array:
    return -V0 * np.exp(-np.abs(x) / a)


def V_soft_barrier(x: Array, V0: float = 1.0) -> Array:
    return V0 / (1.0 + x * x)


def V_piecewise(x: Array, pieces: Tuple[Tuple[float, float, float], ...]) -> Array:
    """Example rough/fractal-like piecewise segments: list of (x0, x1, slope)."""
    out = np.zeros_like(x, dtype=float)
    for x0, x1, s in pieces:
        mask = (x >= x0) & (x < x1)
        out[mask] = s * (x[mask] - x0)
    return out


def V_rough(x: Array, V0: float = 1.0, k0: float = 1.0, levels: int = 4, decay: float = 0.5) -> Array:
    """Deterministic rough/fractal-like potential via multi-scale sines.

    Parameters
    ----------
    V0 : float
        Overall amplitude.
    k0 : float
        Base wave number; higher levels use powers of two.
    levels : int
        Number of harmonic levels to include.
    decay : float
        Amplitude decay per level (geometric).
    """
    x = np.asarray(x, dtype=float)
    V = np.zeros_like(x)
    amp = V0
    k = k0
    for _ in range(max(1, int(levels))):
        V += amp * np.sin(k * x)
        amp *= decay
        k *= 2.0
    return V

def V_Square_Well(x: Array, offset: float = 0.0, width: float = 1.0, depth: float = 1.0) -> Array:
    """Square well potential: V(x) = -depth for |x - offset| < width/2, else 0."""
    out = np.zeros_like(x, dtype=float)
    mask = np.abs(x - offset) < (width / 2)
    out[mask] = -depth
    return out

POTENTIAL_REGISTRY: Dict[str, Callable[..., Array]] = {
    "cusp": V_cusp,
    "exp_well": V_exponential,
    "soft_barrier": V_soft_barrier,
    "piecewise": V_piecewise,
    "rough": V_rough,
    "square": V_Square_Well,
}
