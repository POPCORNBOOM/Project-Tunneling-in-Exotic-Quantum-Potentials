"""High-level task pipelines for the tunneling project."""
from __future__ import annotations

from typing import Dict, Any, Callable
import numpy as np

from .potentials import PotentialSpec
from .grid import GridSpec
from .bound_states import solve_bound_states, SolverSpec
from .fields import apply_field, barrier_top, classify_barrier
from .wkb import find_turning_points, action_integral, wkb_transmission
from .observables import localization_metrics, forbidden_probability, compute_probability, probability_current
from .tdse import run_tdse_frames, build_cap

Array = np.ndarray


def run_bound_states(cfg: Dict[str, Any]) -> Dict[str, Any]:
    pot_spec = PotentialSpec(**cfg["potential"])
    grid = GridSpec(**cfg["grid"])
    solver = SolverSpec(**cfg.get("solver", {}))

    V_func = pot_spec.build()
    x, Vx, E, psi, dx = solve_bound_states(V_func, grid, solver)

    metrics = []
    forb = []
    for n in range(solver.k):
        mean_x, sigma, ipr = localization_metrics(x, psi[:, n], dx)
        metrics.append({"mean_x": mean_x, "sigma": sigma, "ipr": ipr})
        forb.append(forbidden_probability(Vx, E[n], psi[:, n], dx))

    return {
        "x": x,
        "Vx": Vx,
        "E": E,
        "psi": psi,
        "dx": dx,
        "metrics": metrics,
        "forbidden": forb,
    }


def run_wkb_slice(cfg: Dict[str, Any], res: Dict[str, Any], state_index: int = 0, F: float = 0.1) -> Dict[str, Any]:
    x = res["x"]
    Vx = res["Vx"]
    E = res["E"][state_index]
    Vtilt = apply_field(Vx, x, F)
    tps = find_turning_points(x, Vtilt, E, x_min=0.0)
    if tps.size < 2:
        return {"barrier": False, "message": "No closed barrier (likely over-the-barrier)."}
    x1, x2 = tps[0], tps[1]
    S = action_integral(x[(x >= x1) & (x <= x2)], Vtilt[(x >= x1) & (x <= x2)], E)
    T = wkb_transmission(S)
    return {
        "barrier": True,
        "turning_points": (x1, x2),
        "S": S,
        "T_wkb": T,
    }


def run_field_scan(cfg: Dict[str, Any], res: Dict[str, Any], state_index: int = 0, F_grid: Array | None = None) -> Dict[str, Any]:
    x = res["x"]
    Vx = res["Vx"]
    E = float(res["E"][state_index])
    if F_grid is None:
        F_grid = np.linspace(0.01, cfg.get("F_max", 0.5), cfg.get("F_steps", 20))
    records = []
    for F in F_grid:
        Vtilt = apply_field(Vx, x, F)
        tps = find_turning_points(x, Vtilt, E, x_min=0.0)
        status = classify_barrier(E, Vtilt)
        if tps.size >= 2 and status == "not_over_the_barrier":
            x1, x2 = tps[0], tps[1]
            mask = (x >= x1) & (x <= x2)
            S = action_integral(x[mask], Vtilt[mask], E)
            T = wkb_transmission(S)
            xtop, vtop = barrier_top(x, Vtilt)
            records.append({"F": float(F), "status": "not_over_the_barrier", "turning_points": (x1, x2), "S": S, "T": T, "barrier_top": (xtop, vtop)})
        else:
            xtop, vtop = barrier_top(x, Vtilt)
            records.append({"F": float(F), "status": f"over_the_barrier or turning_point_out_of_range", "turning_points": (), "barrier_top": (xtop, vtop)})
    return {"records": records, "F_grid": np.array(F_grid, dtype=float)}



def run_tdse(cfg: Dict[str, Any], res: Dict[str, Any], state_index: int = 0) -> Dict[str, Any]:
    tdse_cfg = cfg.get("tdse", {})
    pot_spec = PotentialSpec(**cfg["potential"])
    base_V = pot_spec.build()
    F = float(tdse_cfg.get("F", 0.0))
    V_func = (lambda x: apply_field(base_V(x), x, F)) if F != 0.0 else base_V

    x = res["x"]
    dx = res["dx"]
    psi0 = res["psi"][:, state_index]
    
    scan = run_field_scan(cfg, res, state_index=state_index, F_grid=[cfg["tdse"]["F"]])
    record = scan['records'][0]

    cap_cfg = tdse_cfg.get("cap")
    cap = None
    cap_markers = {}
    if cap_cfg:
        cap = build_cap(x, **cap_cfg)
        cap_markers = {"left_cap": -abs(cap_cfg.get("x_start", 0.0)), "right_cap": abs(cap_cfg.get("x_start", 0.0))}

    frames = run_tdse_frames(
        psi0,
        V_func,
        x,
        dx,
        duration=float(tdse_cfg["duration"]),
        dt=float(tdse_cfg["dt"]),
        record_interval=int(tdse_cfg.get("record_interval", 10)),
        cap=cap,
        hbar=float(tdse_cfg.get("hbar", 1.0)),
        m=float(tdse_cfg.get("m", 1.0)),
    )

    # Define mask_in_well by calculating the indices of x within the well region, that is res["Vx"] < E[state_index]
    mask_in_well = res["Vx"] < float(res["E"][state_index])
    surv = np.array([compute_probability(f["psi"], dx, mask=mask_in_well) for f in frames], dtype=float)
    #flux = np.array([boundary_flux(f["psi"], dx, hbar=float(tdse_cfg.get("hbar", 1.0)), m=float(tdse_cfg.get("m", 1.0))) for f in frames])
    # calculate the flux at turning points x1 and x2
    x1, x2 = record.get("turning_points", (None, None))
    jxt = np.array([probability_current(f["psi"], dx, hbar=float(tdse_cfg.get("hbar", 1.0)), m=float(tdse_cfg.get("m", 1.0))) for f in frames])
    j1t = np.array([jx[int(np.searchsorted(x, x1))] if x1 is not None else 0.0 for jx in jxt], dtype=float)
    j2t = np.array([jx[int(np.searchsorted(x, x2))] if x2 is not None else 0.0 for jx in jxt], dtype=float)
    total = np.array([compute_probability(f["psi"], dx) for f in frames], dtype=float)
    return {
        "frames": frames,
        "survival": surv,
        "x1": x1,
        "x2": x2,
        "j1": j1t,
        "j2": j2t,
        "total": total,
        "cap": cap,
        "cap_markers": cap_markers,
    }
