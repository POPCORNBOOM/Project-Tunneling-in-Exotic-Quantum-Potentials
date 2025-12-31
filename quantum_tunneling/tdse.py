"""Minimal TDSE propagator with a single, intuitive interface.

Usage:
    frames = run_tdse_frames(psi0, V_func, x, dx, duration, dt, record_interval=10)
    # frames: list of {"t": float, "psi": array}

We provide only a split-operator FFT scheme (fast, Kosloff-style).
CAP (complex absorbing potential) is optional via the `cap` array.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, List, Dict, Any
import numpy as np

Array = np.ndarray


@dataclass(frozen=True)
class PropSpec:
    dt: float
    hbar: float = 1.0
    m: float = 1.0


def build_cap(x: Array, x_start: float, x_end: float, strength: float = 1.0, power: int = 2) -> Array:
    """Symmetric CAP: W(x)=0 inside |x|<x_start, rises to strength at |x|>=x_end."""
    w = np.zeros_like(x, dtype=float)
    mask_left = x <= -x_start
    mask_right = x >= x_start
    if np.any(mask_left):
        w[mask_left] = strength * ((np.abs(x[mask_left]) - x_start) / (x_end - x_start)) ** power
    if np.any(mask_right):
        w[mask_right] = strength * ((np.abs(x[mask_right]) - x_start) / (x_end - x_start)) ** power
    return w


def _k_half_phase(N: int, dx: float, dt: float, hbar: float, m: float) -> Array:
    k = 2.0 * np.pi * np.fft.fftfreq(N, d=dx)
    return np.exp(-1j * (hbar * k * k) / (2.0 * m) * (dt / 2.0))


def _split_step(psi: Array, Vx: Array, k_half: Array, dt: float, hbar: float, cap: Array | None) -> Array:
    V_eff = Vx if cap is None else (Vx - 1j * cap)
    psi_k = np.fft.fft(psi)
    psi_k *= k_half
    psi = np.fft.ifft(psi_k)
    psi *= np.exp(-1j * V_eff * (dt / hbar))
    psi_k = np.fft.fft(psi)
    psi_k *= k_half
    psi = np.fft.ifft(psi_k)
    return psi


def run_tdse_frames(
    psi0: Array,
    Vx: Array,
    dx: float,
    duration: float,
    dt: float,
    record_interval: int = 1,
    cap: Array | None = None,
    hbar: float = 1.0,
    m: float = 1.0,
) -> List[Dict[str, Any]]:
    """Propagate psi in time and return frames.

    Parameters
    ----------
    psi0 : array
        Initial wavefunction on grid x.
    V_func : callable
        V(x) -> array on the same grid.
    x : array
        Spatial grid.
    dx : float
        Grid spacing.
    duration : float
        Total evolution time.
    dt : float
        Time step.
    record_interval : int
        Save every N steps (including t=0).
    cap : array or None
        Optional CAP profile W(x); used as -i W in V_eff.
    hbar, m : float
        Physical constants.

    Returns
    -------
    frames : list of dict
        Each frame has keys: "t" (float), "psi" (ndarray).
    """
    psi = psi0.copy()

    N = psi.size
    n_steps = int(np.ceil(duration / dt))
    k_half = _k_half_phase(N, dx, dt, hbar, m)

    frames: List[Dict[str, Any]] = [{"t": 0.0, "psi": psi.copy()}]

    for step in range(1, n_steps + 1):
        psi = _split_step(psi, Vx, k_half, dt, hbar, cap)
        if step % record_interval == 0:
            frames.append({"t": step * dt, "psi": psi.copy()})

    return frames
