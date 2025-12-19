import numpy as np
from typing import Optional
Array = np.ndarray
# 找根函数

def find_root_points(
    x: Array,
    Fx: Array,
    Y: float,
    x_min: Optional[float] = None,
    x_max: Optional[float] = None,
) -> list:
    x = np.asarray(x, dtype=float)
    Fx = np.asarray(Fx, dtype=float)
    if x.ndim != 1 or Fx.ndim != 1 or x.size != Fx.size:
        raise ValueError("x and Fx must be 1D arrays of the same length.")

    if x[0] > x[-1]:
        x = x[::-1]
        Fx = Fx[::-1]

    mask = np.ones_like(x, dtype=bool)
    if x_min is not None:
        mask &= x >= x_min
    if x_max is not None:
        mask &= x <= x_max

    idx = np.where(mask)[0]
    if idx.size < 2:
        return []

    xw = x[idx]
    fw = Fx[idx] - Y

    roots = []

    near = np.where(np.abs(fw) <= 1e-12)[0]
    if near.size > 0:
        start = near[0]
        prev = near[0]
        for k in near[1:]:
            if k == prev + 1:
                prev = k
            else:
                id_frac = idx[0] + 0.5 * (start + prev)
                nearest = int(round(id_frac))
                x_val = 0.5 * (xw[start] + xw[prev])
                roots.append({"id": float(id_frac), "nearest": nearest, "x": float(x_val)})
                start = prev = k
        id_frac = idx[0] + 0.5 * (start + prev)
        nearest = int(round(id_frac))
        x_val = 0.5 * (xw[start] + xw[prev])
        roots.append({"id": float(id_frac), "nearest": nearest, "x": float(x_val)})

    good = (np.abs(fw[:-1]) > 1e-12) & (np.abs(fw[1:]) > 1e-12)
    sign_change = good & (fw[:-1] * fw[1:] < 0)
    for k in np.where(sign_change)[0]:
        xL, xR = xw[k], xw[k + 1]
        fL, fR = fw[k], fw[k + 1]
        x_root = xL - fL * (xR - xL) / (fR - fL)
        id_frac = idx[0] + k + (x_root - xL) / (xR - xL)
        nearest = int(round(id_frac))
        roots.append({"id": float(id_frac), "nearest": nearest, "x": float(x_root)})

    return roots


if __name__ == "__main__":
    x = np.linspace(-2*np.pi, 2*np.pi, 1000)
    Fx = np.sin(x)
    Y = 0.99999
    roots = find_root_points(x, Fx, Y)
    for r in roots:
        print(r)