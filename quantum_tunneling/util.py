"""Utility helpers."""
import numpy as np
from typing import Optional

Array = np.ndarray


def find_roots(
	x: Array,
	Yx: Array,
	y_target: float,
	x_min: Optional[float] = None,
	x_max: Optional[float] = None,
	atol: float = 1e-12,
	merge_tol: Optional[float] = None,
):
	"""
	通用一维找根：在 `x` 上对 `Yx - y_target` 的零点进行线性插值求解。

	参数
	---
	x, Yx: ndarray
		样本点与对应函数值，需等长一维。
	y_target: float
		目标值，求解 Yx == y_target 的交点。
	x_min, x_max: float | None
		可选区间裁剪。
	atol: float
		视为“已命中零”的绝对公差。
	merge_tol: float | None
		合并近邻根的阈值；缺省时按网格步长中值 * 1.5。

	返回
	---
	List[dict]
		每个根对应一个字典：
		- `x`: 根位置（插值或邻段中点）。
		- `idx`: 最接近该根的原始索引（按裁剪后的 xw）。
	"""
	x = np.asarray(x, dtype=float)
	Yx = np.asarray(Yx, dtype=float)
	if x.ndim != 1 or Yx.ndim != 1 or x.size != Yx.size:
		raise ValueError("x and Yx must be 1D arrays of the same length.")

	# ensure ascending for consistent slicing
	if x[0] > x[-1]:
		x = x[::-1]
		Yx = Yx[::-1]

	mask = np.ones_like(x, dtype=bool)
	if x_min is not None:
		mask &= x >= x_min
	if x_max is not None:
		mask &= x <= x_max

	idx = np.where(mask)[0]
	if idx.size < 2:
		return np.array([], dtype=float)

	xw = x[idx]
	fw = Yx[idx] - y_target

	if merge_tol is None:
		dx = np.median(np.diff(xw))
		merge_tol = 1.5 * dx if dx > 0 else atol

	roots = []  # store (x_root, idx_nearest)

	# exact/near hits
	near = np.where(np.abs(fw) <= atol)[0]
	if near.size > 0:
		start = near[0]
		prev = near[0]
		for k in near[1:]:
			if k == prev + 1:
				prev = k
			else:
				roots.append(0.5 * (xw[start] + xw[prev]))
				start = prev = k
		xr = 0.5 * (xw[start] + xw[prev])
		i_mid = int(0.5 * (start + prev))
		roots.append((xr, int(idx[i_mid])))

	# sign-change crossings (linear interpolation)
	good = (np.abs(fw[:-1]) > atol) & (np.abs(fw[1:]) > atol)
	sign_change = good & (fw[:-1] * fw[1:] < 0)
	for i in np.where(sign_change)[0]:
		x0, x1 = xw[i], xw[i + 1]
		f0, f1 = fw[i], fw[i + 1]
		xr = x0 - f0 * (x1 - x0) / (f1 - f0)
		# choose nearer original index among the bracket
		nearest = idx[i] if abs(f0) <= abs(f1) else idx[i + 1]
		roots.append((xr, int(nearest)))

	if not roots:
		return np.array([], dtype=float)

	# sort by x
	roots = sorted(roots, key=lambda t: t[0])

	merged = [roots[0]]
	for r in roots[1:]:
		if abs(r[0] - merged[-1][0]) > merge_tol:
			merged.append(r)
		else:
			# merge position; keep nearer index to merged position
			x_new = 0.5 * (merged[-1][0] + r[0])
			i_new = merged[-1][1] if abs(x_new - x[merged[-1][1]]) <= abs(x_new - x[r[1]]) else r[1]
			merged[-1] = (x_new, i_new)

	return [{"x": float(xr), "idx": int(i)} for xr, i in merged]
