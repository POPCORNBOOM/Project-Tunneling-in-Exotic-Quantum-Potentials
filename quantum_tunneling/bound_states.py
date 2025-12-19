"""Finite-difference eigensolver for 1D bound states."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Tuple
import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import eigsh

from .grid import GridSpec

Array = np.ndarray


@dataclass(frozen=True)
class SolverSpec:
    k: int = 6
    hbar: float = 1.0
    m: float = 1.0


def solve_bound_states(
    V_func: Callable[[Array], Array],
    grid: GridSpec,
    solver: SolverSpec = SolverSpec(),
) -> Tuple[Array, Array, Array, Array, float]:
    """Return x, V(x), energies, eigenvectors, dx."""
    L, N = grid.L, grid.N
    k, hbar, m = solver.k, solver.hbar, solver.m

    dx = grid.dx
    x = np.linspace(-L + dx, L - dx, N)
    Vx = V_func(x)

    pref = (hbar * hbar) / (2.0 * m * dx * dx)
    main_diag = 2.0 * pref + Vx
    off_diag = -pref * np.ones(N - 1)
    H = diags([off_diag, main_diag, off_diag], offsets=[-1, 0, 1], format="csc")

    E, psi = eigsh(H, k=k, which="SA")
    idx = np.argsort(E)
    E = E[idx]
    psi = psi[:, idx]

    for i in range(k):
        norm = np.sqrt(np.sum(np.abs(psi[:, i]) ** 2) * dx)
        psi[:, i] /= norm

    return x, Vx, E, psi, dx
