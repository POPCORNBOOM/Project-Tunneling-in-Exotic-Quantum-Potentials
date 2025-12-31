"""Turning points and WKB action integrals."""
from __future__ import annotations

from .util import find_roots
from .fields import barrier_top

from typing import Optional, Dict, Any, Tuple
import numpy as np

Array = np.ndarray
# 这个函数可以用来依据V(x)和E_n检查势垒情况。如果存在闭合势垒，则计算WKB透射率。
def barrier_check(
    x: Array,
    Vx: Array,
    E_n: float,
    x_min: float = 0.0,
    x_max: float = np.inf,
) -> Dict[str, Any]:
    """
    参数
    ---
    x: ndarray
        网格点。
    Vx: ndarray
        势能取值。
    E_n: float
        本征能。
    x_min: float
        搜索转折点的 x 起点（默认只取 x>=0）。
    x_max: float
        搜索转折点的 x 终点（默认无穷大）。

    返回
    ---
    Dict[str, Any]
        - `barrier` (bool): 是否存在闭合势垒（至少两转折点且未过顶）。
        - `status` (str): 势垒状态
            - `over_the_barrier`: 能量超过势垒顶。
            - `no_closed_barrier`: 未检测到闭合势垒。
            - `can_tunnel`: 存在闭合势垒，可计算 WKB 透射率。
        * 当 `barrier` 为 True 时，附加键：
            - `turning_points` (Tuple[float, float]): 两转折点位置。
            - `S` (float): WKB 作用量积分。
            - `T_wkb` (float): WKB 透射率 `exp(-2S)`。
            
    """
    tps = find_roots(x, Vx, E_n, x_min=x_min, x_max=x_max)
    
    vtop_idx = barrier_top(x, Vx, x_min=x_min, x_max=x_max)
    xtop = float(x[vtop_idx])
    vtop = float(Vx[vtop_idx])
    
    result = {
            "status": f"over_the_barrier({E_n:.4f} > {vtop:.4f})",
            "barrier": False,
            "barrier_top": (xtop, vtop),
        }
    
    if E_n > vtop: # 越过势垒顶
        return result
        
    if len(tps) < 2: # 不存在闭合势垒，单调V或另一个转折点在考虑区域之外，需要增大L或减小x_min
        result["status"] = "no_closed_barrier"
        return result

    x1, x2 = float(tps[0]["x"]), float(tps[-1]["x"])
    mask = (x >= x1) & (x <= x2)
    S = float(action_integral(x[mask], Vx[mask], E_n))
    T = float(wkb_transmission(S))
    result["status"] = f"can_tunnel({E_n:.4f} < {vtop:.4f}, Tunneling Points: {x1:.4f} to {x2:.4f})"
    result["barrier"] = True
    result["turning_points"] = (x1, x2)
    result["S"] = S
    result["T_wkb"] = T
    return result

# Helper function to calculate action integral
def action_integral(x: Array, Vx: Array, E: float, m: float = 1.0) -> float:
    mask = Vx > E
    xb = x[mask]
    barrier = Vx[mask] - E
    integrand = np.sqrt(np.maximum(0.0, 2.0 * m * barrier))
    return float(np.trapz(integrand, xb))


def wkb_transmission(S: float, hbar: float = 1.0) -> float:
    return float(np.exp(-2.0 * S / hbar))
