"""Plotting helpers (minimal scaffolding)."""
from __future__ import annotations

import warnings

import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.animation import FuncAnimation
from matplotlib.colors import Normalize
from typing import List, Dict, Any, Optional
from IPython.display import HTML
from tqdm.auto import tqdm

import numpy as np
from .fields import apply_field

Array = np.ndarray


def plot_potential_and_states(x: Array, Vx: Array, E: Array, psi: Array, scale: float = 0.5, F = 0):
    plt.figure()
    
    plt.plot(x, Vx, label="V(x)")
    
    for n, En in enumerate(E):
        if n < 8:
            plt.plot(x, scale * psi[:, n] + En, label=f"psi_{n} + E_{n}")
        else:
            plt.plot(x, scale * psi[:, n] + En, color="gray", alpha=0.5)
            
    if F != 0:
        Vtilt = apply_field(Vx, x, F)
        plt.plot(x, Vtilt, label=f"V(x) - F*x (F={F})")
        
    plt.xlabel("x")
    plt.ylabel("Energy")
    plt.legend()
    plt.tight_layout()


def plot_turning_points(x: Array, Vx: Array, E: float, tps: Array):
    plt.figure()
    plt.plot(x, Vx, label="V(x)")
    plt.plot(x, np.full_like(x, E), linestyle="--", label="E")
    for tp in tps:
        plt.axvline(tp, color="r", linestyle=":")
    plt.legend()
    plt.tight_layout()
def create_wavefunction_animation(
    frames: List[Dict],
    config: Dict[str, Any],
    interval: int = 50,
    save_path: Optional[str] = None,
    x: Optional[Array] = None,
    Vx: Optional[Array] = None,
    psi_scale: float = 1.0,
    psi_offset: float = 0.0,
    markers: Optional[Dict[str, float]] = None,
    figsize: tuple = (10, 6),
    notebook: bool = True,
    progress: bool = True,
    phase_stride: int = 1,
    generate_progress: bool = True,

):
    """简化版波函数/概率密度动画。

    必需的 config:
    - mode: "probability density" 或 "wave function"
    - 若 mode=="wave function"，可指定 real/im/magnitude 颜色 (None 表示不画)，phase(True/False) 是否用相位着色填充。
    - phase_stride: 相位填充的抽样步长，>1 可减少多边形数量以提升速度。
    - Vx: 可选势能曲线（长度需与 psi 一致）；不再缩放 Vx，而是通过 psi_scale/psi_offset 调整波函数以对齐势能。
    - psi_scale / psi_offset: 对波函数与密度的整体缩放与平移（默认不变）。
    - markers: 可选 dict，支持键 left_cap/right_cap/left_tp/right_tp 对应四条竖线位置（x 值）。
    """

    valid_modes = ["probability density", "wave function"]
    if config.get("mode") not in valid_modes:
        raise ValueError(f"mode 必须是 {valid_modes} 之一")
    if config["mode"] == "wave function":
        for key in ["real", "im", "magnitude", "phase"]:
            if key not in config:
                raise ValueError(f"wave function 模式需要 '{key}' 参数")

    psi0 = np.asarray(frames[0]["psi"])
    n_points = len(psi0)
    if x is None:
        x = np.arange(n_points)
    else:
        x = np.asarray(x)
        if x.shape[0] != n_points:
            raise ValueError("x 的长度必须与 psi 匹配")

    Vx_plot = None
    if Vx is not None:
        Vx_arr = np.asarray(Vx)
        if Vx_arr.shape[0] != n_points:
            raise ValueError("Vx 的长度必须与 psi 匹配")
        Vx_plot = Vx_arr.astype(float)

    if notebook and not save_path:
        try:
            get_ipython().run_line_magic("matplotlib", "inline")
        except Exception:
            pass

    fig, ax = plt.subplots(figsize=figsize)

    lines: List[Any] = []
    phase_fills: List[Any] = []
    pot_line = None
    marker_lines: List[Any] = []

    # y 轴范围基于 psi（先不含 Vx，后续合并），并对 psi 应用 scale/offset
    if config["mode"] == "probability density":
        probs_all = [(np.abs(np.asarray(f["psi"])) ** 2) * psi_scale + psi_offset for f in frames]
        y_min = 0.0
        y_max = float(np.max(probs_all)) if len(probs_all) else 1.0
        y_max = 1.0 if y_max == 0 else y_max * 1.1
    else:
        real_vals = []
        imag_vals = []
        abs_vals = []
        for f in frames:
            psi_arr = np.asarray(f["psi"])
            real_vals.extend([(psi_arr.real * psi_scale + psi_offset).min(), (psi_arr.real * psi_scale + psi_offset).max()])
            imag_vals.extend([(psi_arr.imag * psi_scale + psi_offset).min(), (psi_arr.imag * psi_scale + psi_offset).max()])
            abs_vals.extend([(np.abs(psi_arr) * psi_scale + psi_offset).min(), (np.abs(psi_arr) * psi_scale + psi_offset).max()])
        y_min = min(real_vals + imag_vals + abs_vals)
        y_max = max(real_vals + imag_vals + abs_vals)
        padding = (y_max - y_min) * 0.1 if y_max > y_min else 1.0
        y_min, y_max = y_min - padding, y_max + padding

    if config["mode"] == "probability density":
        line_prob, = ax.plot([], [], "b-", linewidth=2, label="|ψ|^2")
        lines.append(line_prob)
        ax.set_ylabel("Probability Density", fontsize=12)
        ax.set_title("Probability Density Evolution", fontsize=14, fontweight="bold")
    else:
        if config.get("real") is not None:
            line_r, = ax.plot([], [], color=config["real"], linewidth=2, label="Re(ψ)")
            lines.append(line_r)
        if config.get("im") is not None:
            line_i, = ax.plot([], [], color=config["im"], linewidth=2, label="Im(ψ)")
            lines.append(line_i)
        if config.get("magnitude") is not None:
            line_mag, = ax.plot([], [], color=config["magnitude"], linewidth=2, label="|ψ|")
            lines.append(line_mag)
        if lines:
            ax.legend(loc="upper right")
        ax.set_ylabel("Amplitude", fontsize=12)
        ax.set_title("Wave Function Evolution", fontsize=14, fontweight="bold")

        phase_norm = Normalize(vmin=0, vmax=2 * np.pi) if config.get("phase") else None
        cmap = cm.hsv if config.get("phase") else None

    # 若提供 Vx，则直接绘制（不缩放势能），波函数通过 psi_scale/psi_offset 自行对齐
    if Vx_plot is not None:
        y_min = min(y_min, Vx_plot.min())
        y_max = max(y_max, Vx_plot.max())
        pot_line, = ax.plot(x, Vx_plot, color="gray", alpha=0.6, linestyle="--", label="V(x)")
        # 保证图例包含势能线
        if config["mode"] == "probability density" or lines:
            ax.legend(loc="upper right")

    ax.set_ylim(y_min, y_max)
    ax.set_xlim(x.min(), x.max())
    ax.set_xlabel("Position", fontsize=12)
    ax.grid(True, alpha=0.3)

    # 竖线标记
    marker_styles = {
        "left_cap": {"color": "#2ca02c", "linestyle": ":", "label": "CAP L"},
        "right_cap": {"color": "#2ca02c", "linestyle": ":", "label": "CAP R"},
        "left_tp": {"color": "#d62728", "linestyle": "--", "label": "TP L"},
        "right_tp": {"color": "#d62728", "linestyle": "--", "label": "TP R"},
    }
    if markers:
        for key, style in marker_styles.items():
            if key in markers:
                line = ax.axvline(markers[key], color=style["color"], linestyle=style["linestyle"], alpha=0.7, label=style["label"])
                marker_lines.append(line)
        if marker_lines:
            ax.legend(loc="upper right")
    time_text = ax.text(
        0.02,
        0.95,
        "",
        transform=ax.transAxes,
        fontsize=12,
        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.8),
    )

    def init():
        for ln in lines:
            ln.set_data([], [])
        for fill in phase_fills:
            try:
                fill.remove()
            except Exception:
                pass
        phase_fills.clear()
        time_text.set_text("")
        return lines + [time_text] + phase_fills

    def update(frame_idx):
        frame = frames[frame_idx]
        psi = np.asarray(frame["psi"])
        t = frame["t"]

        if config["mode"] == "probability density":
            prob = np.abs(psi) ** 2
            prob = prob * psi_scale + psi_offset
            lines[0].set_data(x, prob)
        else:
            line_idx = 0
            if config.get("real") is not None:
                lines[line_idx].set_data(x, psi.real * psi_scale + psi_offset)
                line_idx += 1
            if config.get("im") is not None:
                lines[line_idx].set_data(x, psi.imag * psi_scale + psi_offset)
                line_idx += 1
            if config.get("magnitude") is not None:
                mag = np.abs(psi) * psi_scale + psi_offset
                lines[line_idx].set_data(x, mag)
                line_idx += 1

                if config.get("phase"):
                    phase = np.angle(psi) % (2 * np.pi)
                    for fill in phase_fills:
                        try:
                            fill.remove()
                        except Exception:
                            pass
                    phase_fills.clear()
                    step = max(1, int(phase_stride))
                    for i in range(0, len(x) - 1, step):
                        avg_phase = 0.5 * (phase[i] + phase[min(i + step, len(x) - 1)])
                        color = cmap(phase_norm(avg_phase))
                        fill = ax.fill_between(
                            x[i : i + step + 1], 0, mag[i : i + step + 1], color=color, alpha=0.3, edgecolor="none"
                        )
                        phase_fills.append(fill)

        time_text.set_text(f"Time: {t:.3f}")
        return lines + [time_text] + phase_fills

    frame_indices = range(len(frames))
    if generate_progress:
        frame_indices = tqdm(frame_indices, desc="Generating animation", leave=False)

    ani = FuncAnimation(fig, update, frames=frame_indices, init_func=init, interval=interval, blit=False, repeat=True)

    if save_path:
        progress_cb = None
        pbar = None
        if progress:
            pbar = tqdm(total=len(frames), desc="Saving animation", leave=False)

            def _progress(i, n):
                if pbar.total != n:
                    pbar.total = n
                pbar.n = i
                pbar.refresh()
                if i >= n:
                    pbar.close()

            progress_cb = _progress

        try:
            if save_path.endswith(".gif"):
                ani.save(save_path, writer="pillow", fps=1000 / interval, progress_callback=progress_cb)
            elif save_path.endswith(".mp4"):
                ani.save(save_path, writer="ffmpeg", fps=1000 / interval, progress_callback=progress_cb)
            else:
                ani.save(save_path, fps=1000 / interval, progress_callback=progress_cb)
            print(f"动画已保存到: {save_path}")
        finally:
            if pbar is not None:
                pbar.close()

        plt.close(fig)
        return ani

    if notebook:
        plt.close(fig)
        try:
            return HTML(ani.to_jshtml())
        except Exception:
            update(0)
            plt.show()
            return None

    plt.show()
    return ani

