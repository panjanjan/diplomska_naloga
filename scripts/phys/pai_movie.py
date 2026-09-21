#!/usr/bin/env python3
"""Zgenerira MP4 z dinamiko glavnih osi vztrajnosti domen.

Vhod:  TEST/d{n}_p{k}.xvg (6 datotek, vsaka: čas [ps], vx, vy, vz — enotski vektor)
Izhod: outputs/pai_movie.mp4

Frame-i se rišejo vzporedno po več procesih (--workers); posnetek se nato
sestavi z ffmpeg.
"""
# ------------------------------------------------------------------------------
import argparse
import glob
import io
import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
from mpl_toolkits.mplot3d.art3d import Line3DCollection

_G = None


def read_axes(pattern):
    """Prebere d{n}_p{k}.xvg in vrne dict: (domain, axis) -> matrika len(t)×3."""
    files = sorted(glob.glob(pattern))
    if len(files) != 6:
        raise SystemExit(f"Pričakovano 6 datotek, najdenih {len(files)}: {files}")

    axes = {}
    for f in files:
        dat = np.loadtxt(f)
        time = dat[:, 0]
        vec = dat[:, 1:4]
        key = os.path.basename(f).removesuffix(".xvg")
        axes[key] = (time, vec)
    return axes


def unwrap_sign(vecs):
    """Odstrani poljubnost predznaka: vrne kontinuiran potek smeri."""
    cur = vecs[0]
    out = np.empty_like(vecs)
    out[0] = cur
    for i in range(1, len(vecs)):
        nxt = vecs[i]
        if np.dot(cur, nxt) < 0:
            nxt = -nxt
        out[i] = nxt
        cur = nxt
    return out


def arrow_segments(v, head_ratio=0.15, head_width=0.07):
    """Puščica od izhodišča do v: gred + glava kot 3 segmenta."""
    v = np.asarray(v, float)
    tip = v
    shaft = np.array([[0.0, 0.0, 0.0], tip])
    norm = float(np.linalg.norm(v))
    if norm == 0:
        return [shaft]
    axis = v / norm
    ref = np.array([1.0, 0.0, 0.0]) if abs(axis[0]) < 0.9 \
        else np.array([0.0, 1.0, 0.0])
    perp = np.cross(axis, ref)
    perp = perp / np.linalg.norm(perp)
    base = tip - axis * head_ratio * norm
    left = base + perp * head_width
    right = base - perp * head_width
    return [shaft, np.array([tip, left]), np.array([tip, right])]


def _init_worker(state):
    global _G
    _G = state


def _draw_background(ax, keys, trails, colors):
    ax.set_box_aspect((1, 1, 1))
    for a in (ax.xaxis, ax.yaxis, ax.zaxis):
        a.set_pane_color((1, 1, 1))
        a.set_ticklabels([])
    for lab in (ax.set_xlabel, ax.set_ylabel, ax.set_zlabel):
        lab("")
    ax.set_xlim(-1, 1); ax.set_ylim(-1, 1); ax.set_zlim(-1, 1)
    for k in keys:
        tpath = trails[k]
        ax.plot(tpath[:, 0], tpath[:, 1], tpath[:, 2],
                color=colors[k], alpha=0.25, lw=0.7, zorder=1)


def _render_chunk(chunk_idxs):
    """Poriše del frame-ov in vrne seznam PNG-bajtov (v vrstnem redu)."""
    axes, keys, trails, colors, labels, time_all, dpi, rotate = _G
    base_azim, base_elev = -60, 30

    fig = plt.figure(figsize=(6, 6), dpi=dpi)
    ax = fig.add_subplot(111, projection="3d")
    _draw_background(ax, keys, trails, colors)

    proxies = [Line2D([0], [0], color=colors[k], lw=3, label=labels[k])
               for k in keys]
    ax.legend(handles=proxies, loc="upper left", fontsize=8, frameon=True)
    time_txt = ax.text2D(0.02, 0.95, "", transform=ax.transAxes, fontsize=12)

    arrows = {}
    for k in keys:
        coll = Line3DCollection([np.zeros((2, 3))], colors=colors[k],
                                linewidths=3, zorder=5)
        ax.add_collection3d(coll)
        arrows[k] = coll

    pngs = []
    for j, i in chunk_idxs:
        for k in keys:
            arrows[k].set_segments(arrow_segments(axes[k][1][i]))
        ax.view_init(elev=base_elev, azim=base_azim + j * rotate)
        time_txt.set_text(f"t = {time_all[i] / 1000:.1f} ns"
                          f"  ({i} / {len(time_all)})")
        buf = io.BytesIO()
        fig.savefig(buf, format="png")
        pngs.append(buf.getvalue())
    plt.close(fig)
    return pngs


def _build_pipeline(out, codec, fps):
    cmd = ["ffmpeg", "-y", "-loglevel", "error",
           "-f", "image2pipe", "-vcodec", "png", "-framerate", str(fps),
           "-i", "pipe:",
           "-c:v", codec, "-b:v", "2000k", "-pix_fmt", "yuv420p", out]
    return subprocess.Popen(cmd, stdin=subprocess.PIPE)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in-pattern", default="TEST/d*_p*.xvg",
                    help="globska maska vhodnih datotek")
    ap.add_argument("--out", default="outputs/pai_movie.mp4")
    ap.add_argument("--step", type=int, default=2,
                    help="vsaka n-ta frame vnese v posnetek")
    ap.add_argument("--fps", type=int, default=30)
    ap.add_argument("--workers", type=int, default=0,
                    help="število procesov (0 = avtomatsko)")
    ap.add_argument("--chunk", type=int, default=50,
                    help="koliko frame-ov poriše en proces hkrati")
    ap.add_argument("--codec", default="libopenh264",
                    help="ffmpeg video kodek (brez libx264 v ffmpeg-free)")
    ap.add_argument("--dpi", type=int, default=120)
    ap.add_argument("--rotate", type=float, default=0.0,
                    help="vrtenje kamere [°/frame], 0 = fiksna")
    ap.add_argument("--no-unwrap", action="store_true",
                    help="ne popravljaj predznaka osi")
    args = ap.parse_args()

    axes = read_axes(args.in_pattern)
    keys = sorted(axes)
    time_all, _ = axes[keys[0]]
    idxs = np.arange(0, len(time_all), args.step)

    colors = {
        "d1_p1": "#1f77b4", "d1_p2": "#6fb1de", "d1_p3": "#c6dcef",
        "d2_p1": "#d62728", "d2_p2": "#e89b9b", "d2_p3": "#f2cccc",
    }
    labels = {
        "d1_p1": "Domena 1 · os 1", "d1_p2": "Domena 1 · os 2",
        "d1_p3": "Domena 1 · os 3", "d2_p1": "Domena 2 · os 1",
        "d2_p2": "Domena 2 · os 2", "d2_p3": "Domena 2 · os 3",
    }

    if not args.no_unwrap:
        axes = {k: (t, unwrap_sign(v)) for k, (t, v) in axes.items()}

    trails = {k: v[idxs] for k, (_, v) in axes.items()}

    os.makedirs(os.path.dirname(args.out), exist_ok=True)

    # razdeli frame-e na večje kupčke za vzporedno risanje
    frame_tasks = [(j, i) for j, i in enumerate(idxs)]
    chunks = [frame_tasks[s:s + args.chunk]
              for s in range(0, len(frame_tasks), args.chunk)]

    state = (axes, keys, trails, colors, labels,
             time_all, args.dpi, args.rotate)
    n_workers = args.workers or min(os.cpu_count() or 1, len(chunks))

    proc = _build_pipeline(args.out, args.codec, args.fps)
    done = 0
    n_frames = len(frame_tasks)
    print(f"{n_frames} frame-ov, {len(chunks)} kupčkov, {n_workers} procesov",
          flush=True)

    if n_workers > 1:
        with ProcessPoolExecutor(max_workers=n_workers, initializer=_init_worker,
                                 initargs=(state,)) as ex:
            for pngs in ex.map(_render_chunk, chunks):
                for png in pngs:
                    proc.stdin.write(png)
                done += len(pngs)
                print(f"{done}/{n_frames}", flush=True)
    else:
        _init_worker(state)
        for chunk in chunks:
            pngs = _render_chunk(chunk)
            for png in pngs:
                proc.stdin.write(png)
            done += len(pngs)
            print(f"{done}/{n_frames}", flush=True)

    proc.stdin.close()
    rc = proc.wait()
    if rc:
        sys.exit(f"ffmpeg je končal z napako: rc={rc}")
    print(f"zapisano: {args.out}")


if __name__ == "__main__":
    main()