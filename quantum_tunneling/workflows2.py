"""High-level task pipelines for the tunneling project."""
from __future__ import annotations

from typing import Dict, Any, Callable
import numpy as np

from .potentials import PotentialSpec
from .grid import GridSpec
from .bound_states import solve_bound_states, SolverSpec
from .fields import apply_field, barrier_top, classify_barrier
from .wkb import find_roots, action_integral, wkb_transmission
from .observables import localization_metrics, forbidden_probability, compute_probability, probability_current
from .tdse import run_tdse_frames, build_cap

Array = np.ndarray


def run_bound_states(cfg: Dict[str, Any]) -> Dict[str, Any]:
    """
    参数
    ---
    cfg: Dict[str, Any]
        高层配置：
        - `potential`: 传递给 `PotentialSpec(**kwargs)` 的字典（势函数类型及其参数）。
        - `grid`: 传递给 `GridSpec(L, N)` 的字典，`L` 为半区间长度，`N` 为内部点数，步长 `dx=2L/(N+1)`。
        - `solver` (可选): 传递给 `SolverSpec(k=6, hbar=1.0, m=1.0)` 的字典，`k` 为要求的本征态数量。

    返回
    ---
    Dict[str, Any]
        - `x` (ndarray, shape (N,)): 均匀网格点。
        - `Vx` (ndarray, shape (N,)): 势能取值。
        - `E` (ndarray, shape (k,)): 由小到大排序的本征能量。
        - `psi` (ndarray, shape (N, k)): 归一化本征态，中心点相位已统一为正。
        - `dx` (float): 网格间距。
        - `metrics` (List[dict]): 长度为 k，每个字典含 `mean_x`, `sigma`, `ipr`。
        - `forbidden` (List[float]): 每个本征态在禁阻区的概率质量。
    """
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
    """
    参数
    ---
    cfg: Dict[str, Any]
        高层配置，需包含 `potential` 与 `grid` 定义（用于势和网格范围参考）。
    res: Dict[str, Any]
        `run_bound_states` 的输出，必须含 `x`, `Vx`, `E`。
    state_index: int
        选择使用的本征态索引（0 <= index < len(E)）。
    F: float
        均匀外场强度，势按 `Vx - F*x` 线性倾斜。

    返回
    ---
    Dict[str, Any]
        - `barrier` (bool): 是否检测到闭合势垒（至少两个转折点）。
        - `turning_points` (tuple): `(x1, x2)`，当 `barrier=True` 时为升序浮点；否则为空。
        - `S` (float): 作用量积分，仅在 `barrier=True` 时给出。
        - `T_wkb` (float): WKB 透射率 `exp(-2S)`，仅在 `barrier=True` 时给出。
        - `message` (str): 当 `barrier=False` 时的说明（如过势垒）。
    """
    x = res["x"]
    Vx = res["Vx"]
    E = res["E"][state_index]
    Vtilt = apply_field(Vx, x, F)
    tps = find_roots(x, Vtilt, E, x_min=0.0)
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
    """
    参数
    ---
    cfg: Dict[str, Any]
        高层配置；若未提供 `F_grid`，将从中读取 `F_max`, `F_steps` 生成线性外场网格（默认 0.01 到 0.5）。
    res: Dict[str, Any]
        `run_bound_states` 输出，需含 `x`, `Vx`, `E`。
    state_index: int
        选择扫描的本征态索引（用于参考能量）。
    F_grid: Array | None
        可选外场网格；若为 `None` 自动生成等距网格。

    返回
    ---
    Dict[str, Any]
        - `records` (List[dict]): 长度与 `F_grid` 相同。每条记录包含：
          - `F` (float): 外场值。
          - `status` (str): 势垒状态，`not_over_the_barrier` 或 `over_the_barrier or turning_point_out_of_range`。
          - `turning_points` (tuple): 当存在两转折点且未过势垒时给出 `(x1, x2)`；否则为空 tuple。
          - `S` (float): 作用量积分，仅在未过势垒时提供。
          - `T` (float): WKB 透射率，仅在未过势垒时提供。
          - `barrier_top` (tuple): `(x_top, V_top)`，总是提供。
        - `F_grid` (ndarray, shape (nF,)): 实际使用的外场网格。
    """
    x = res["x"]
    Vx = res["Vx"]
    E = float(res["E"][state_index])
    if F_grid is None:
        F_grid = np.linspace(0.01, cfg.get("F_max", 0.5), cfg.get("F_steps", 20))
    records = []
    for F in F_grid:
        Vtilt = apply_field(Vx, x, F)
        tps = find_roots(x, Vtilt, E, x_min=0.0)
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
    """
    参数
    ---
    cfg: Dict[str, Any]
        高层配置，需包含 `potential` 与 `tdse`：
        - `potential`: 传递给 `PotentialSpec(**kwargs)` 的字典。
        - `tdse`: 必需字段 `F`(外场)、`duration`、`dt`；可选 `record_interval`(默认10)、`hbar`、`m`、`cap`。
          其中 `cap` 将直接传给 `build_cap(x, x_start, x_end, strength, power)` 构造 CAP。
    res: Dict[str, Any]
        `run_bound_states` 输出，用作初态与网格参考（需含 `x`, `dx`, `psi`, `E`, `Vx`）。
    state_index: int
        选择作为初态的本征态索引。

    返回
    ---
    Dict[str, Any]
        - `frames` (List[dict]): 由 `run_tdse_frames` 给出的时间序列，每帧含 `t`, `psi`。
        - `survival` (ndarray): 在势阱区域 (V < E_state) 的生存概率随时间的序列。
        - `x1`, `x2` (float | None): 对应外场下的转折点（若存在）；来自单点 `run_field_scan` 结果。
        - `cap` (ndarray | None): 复吸收势轮廓；未启用时为 None。
        - `cap_markers` (dict): 便于可视化的左右 CAP 起点标记，键 `left_cap`/`right_cap`。
    """
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
    return {
        "frames": frames,
        "survival": surv,
        "x1": x1,
        "x2": x2,
        "cap": cap,
        "cap_markers": cap_markers,
    }
