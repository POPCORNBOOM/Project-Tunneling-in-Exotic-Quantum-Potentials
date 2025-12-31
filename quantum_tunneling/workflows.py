"""更干净的 workflow 重写版本。

功能覆盖原 workflows 中的主流程：本征态求解、WKB 单点评估、外场扫描、TDSE 演化。
"""
from __future__ import annotations

from typing import Dict, Any, Iterable, Tuple
from dataclasses import dataclass
import numpy as np

from .potentials import PotentialSpec
from .grid import GridSpec
from .bound_states import solve_bound_states, SolverSpec
from .fields import apply_field, barrier_top
from .wkb import barrier_check
from .observables import localization_metrics, forbidden_probability, compute_probability
from .tdse import run_tdse_frames, build_cap

Array = np.ndarray


def _require(cfg: Dict[str, Any], key: str) -> Any:
    if key not in cfg:
        raise KeyError(f"Missing required config field: {key}")
    return cfg[key]


def _as_float_array(vals: Iterable[float]) -> Array:
    arr = np.asarray(list(vals), dtype=float).ravel()
    if arr.size == 0:
        raise ValueError("F_grid is empty; provide at least one field value.")
    return arr


def run_bound_states(cfg: Dict[str, Any]) -> Dict[str, Any]:
    """
    参数
    ---
    cfg: Dict[str, Any]
        - `potential`: 传入 `PotentialSpec(**kwargs)` 的字典。
        - `grid`: 传入 `GridSpec(L, N)` 的字典，网格步长 `dx=2L/(N+1)`。
        - `solver` (可选): 传入 `SolverSpec(k, hbar, m)` 的字典。

    返回
    ---
    Dict[str, Any]
        - `x` (ndarray): 网格点。
        - `Vx` (ndarray): 势能取值。
        - `E` (ndarray): 升序本征能。
        - `psi` (ndarray): 列为本征态，已归一化且中心点相位为正。
        - `dx` (float): 网格间距。
        - `metrics` (List[dict]): 每个本征态的 `mean_x`, `sigma`, `ipr`。
        - `forbidden` (List[float]): 禁阻区概率质量。
    """
    pot_cfg = _require(cfg, "potential")
    grid_cfg = _require(cfg, "grid")

    pot_spec = PotentialSpec(**pot_cfg)
    grid = GridSpec(**grid_cfg)
    solver = SolverSpec(**cfg.get("solver", {}))

    V_func = pot_spec.build()
    x, Vx, E, psi, dx = solve_bound_states(V_func, grid, solver)

    metrics = []
    forbidden = []
    for idx in range(solver.k):
        mean_x, sigma, ipr = localization_metrics(x, psi[:, idx], dx)
        metrics.append({"mean_x": mean_x, "sigma": sigma, "ipr": ipr})
        forbidden.append(forbidden_probability(Vx, E[idx], psi[:, idx], dx))

    return {
        "x": x,
        "Vx": Vx,
        "E": E,
        "psi": psi,
        "dx": dx,
        "metrics": metrics,
        "forbidden": forbidden,
    }




def run_field_scan(
    cfg: Dict[str, Any],
    res: Dict[str, Any],
    state_index: int = 0,
    F_grid: Iterable[float] | None = None,
    x_min: float = 0.0,
) -> Dict[str, Any]:
    """
    参数
    ---
    cfg: Dict[str, Any]
        用于读取默认外场范围：`F_max`(默认0.5)、`F_steps`(默认20)。
    res: Dict[str, Any]
        `run_bound_states` 输出，需含 `x`, `Vx`, `E`。
    state_index: int
        使用的本征态索引。
    F_grid: Iterable[float] | None
        外场网格；未提供时用 cfg 构造线性网格（含端点）。
    x_min: float
        转折点搜索的下限。

    返回
    ---
    Dict[str, Any]
        - `records` (List[dict]): 每个外场的 WKB 结果，附加键 `F`。
        - `F_grid` (ndarray): 使用的外场序列。
    """
    if F_grid is None:
        F_max = float(cfg.get("F_max", 0.5))
        F_steps = int(cfg.get("F_steps", 20))
        if F_steps < 1:
            raise ValueError("F_steps must be >= 1")
        F_grid_arr = np.linspace(0.01, F_max, F_steps)
    else:
        F_grid_arr = _as_float_array(F_grid)

    records = []
    for F in F_grid_arr:
        x = res["x"]
        Vx = res["Vx"]
        Vtilt = apply_field(Vx, x, F)
        E_n = float(res["E"][state_index])
        rec = barrier_check(x, Vtilt, E_n, x_min=x_min)
        rec["F"] = F
        records.append(rec)

    return {"records": records, "F_grid": F_grid_arr}


def run_tdse(
    cfg: Dict[str, Any],
    res: Dict[str, Any],
    state_index: int = 0,
) -> Dict[str, Any]:
    """
    参数
    ---
    cfg: Dict[str, Any]
        必需字段：
        - `potential`: 传入 `PotentialSpec(**kwargs)` 的字典。
        - `tdse`: 包含 `F`, `duration`, `dt`；可选 `record_interval`(默认10)、`hbar`, `m`, `cap`。
    res: Dict[str, Any]
        `run_bound_states` 输出，需含 `x`, `dx`, `psi`, `E`, `Vx`。
    state_index: int
        作为初态的本征态索引。

    返回
    ---
    Dict[str, Any]
        - `frames`: TDSE 演化帧列表，每帧含 `t`, `psi`。
        - `survival`: 生存概率序列（势阱区域内）。
        - `wkb`: 与所选外场对应的 WKB 记录（含转折点等）。
        - `cap`: CAP 轮廓或 None。
        - `cap_markers`: 左右 CAP 起点标记，便于绘图。
    """
    tdse_cfg = _require(cfg, "tdse")
    pot_cfg = _require(cfg, "potential")

    pot_spec = PotentialSpec(**pot_cfg)
    base_V = pot_spec.build()

    F = float(tdse_cfg.get("F", 0.0))
    V_func = (lambda x: apply_field(base_V(x), x, F)) if F != 0.0 else base_V

    x = res["x"]
    dx = float(res["dx"])
    psi0 = res["psi"][:, state_index]

    cap_cfg = tdse_cfg.get("cap")
    cap = build_cap(x, **cap_cfg) if cap_cfg else None
    cap_markers = {}
    if cap_cfg:
        x_start = abs(cap_cfg.get("x_start", 0.0))
        cap_markers = {"left_cap": -x_start, "right_cap": x_start}

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

    mask_in_well = res["Vx"] < float(res["E"][state_index])
    survival = np.array([compute_probability(f["psi"], dx, mask=mask_in_well) for f in frames], dtype=float)

    barrier_info = barrier_check(x, apply_field(res["Vx"], x, F), float(res["E"][state_index]), x_min=0.0)

    x1, x2 = barrier_info.get("turning_points", (None, None)) if barrier_info.get("barrier") else (None, None)

    return {
        "frames": frames,
        "survival": survival,
        "barrier_info": barrier_info,
        "x1": x1,
        "x2": x2,
        "cap": cap,
        "cap_markers": cap_markers,
    }
