#!/usr/bin/env python3
"""MP4 film dinamike proteina + prve glavne osi vztrajnosti domen.

Vhod:
  * PDB      topologija proteina                  (--pdb)
  * XTC      trajektorija (že centrirana)         (--xtc)
  * d1_p1.xvg / d2_p1.xvg smeri prve osi [ps, vx, vy, vz]   (--in-pattern)

Prva glavna os (enotski vektor iz xvg) se nariše skozi masni center domene,
protein kot točke težkih atomov. Frame-i se rišejo vzporedno (--workers).
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

import MDAnalysis as mda

_G = None

DOMAIN_STYLE = {
    "d1": {"axis": "#164e7d", "points": "#7da6c9"},
    "d2": {"axis": "#c0392b", "points": "#d9a4a4"},
}


def read_axes(pattern):
    """Prebere d{n}_p1.xvg; vrne matrike vektorjev per prefix in čas."""
    files = sorted(glob.glob(pattern))
    if not files:
        raise SystemExit(f"Najdenih ni bilo datotek po maski: {pattern}")
    axes = {}
    time = None
    for f in files:
        dat = np.loadtxt(f)
        key = os.path.basename(f).split("_")[0]
        axes[key] = dat[:, 1:4]
        if time is None:
            time = dat[:, 0]
    return axes, time


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


def axis_segments(uvec, center, length):
    """Segmenta osi skozi center: gred + glavi na obeh koncih."""
    axis = np.asarray(uvec, float)
    norm = np.linalg.norm(axis)
    if norm == 0:
        return [np.array([center, center])]
    axis = axis / norm
    ref = np.array([1.0, 0.0, 0.0]) if abs(axis[0]) < 0.9 \
        else np.array([0.0, 1.0, 0.0])
    perp = np.cross(axis, ref)
    perp = perp / np.linalg.norm(perp)
    half = length / 2.0
    head = 0.18 * length
    width = 0.10 * half
    tip_p = center + axis * half
    tip_m = center - axis * half
    base_p = tip_p - axis * head
    base_m = tip_m + axis * head
    return [
        np.array([center - axis * half, center + axis * half]),
        np.array([tip_p, base_p + perp * width]),
        np.array([tip_p, base_p - perp * width]),
        np.array([tip_m, base_m + perp * width]),
        np.array([tip_m, base_m - perp * width]),
    ]


def _init_worker(state):
    global _G
    _G = state


def _render_chunk(chunk_idxs):
    """Poriše del frame-ov in vrne seznam PNG-bajtov (v vrstnem redu)."""
    st = _G
    u = mda.Universe(st["pdb"], st["xtc"])
    leaf = u.atoms

    fig = plt.figure(figsize=(6, 6), dpi=st["dpi"])
    ax = fig.add_subplot(111, projection="3d")
    bx = st["box"]
    ax.set_box_aspect((bx[0], bx[1], bx[2]))
    for a in (ax.xaxis, ax.yaxis, ax.zaxis):
        a.set_pane_color((1, 1, 1))
        a.set_ticklabels([])
    for lab in (ax.set_xlabel, ax.set_ylabel, ax.set_zlabel):
        lab("")
    ax.set_xlim(*st["ranges"][0]); ax.set_ylim(*st["ranges"][1])
    ax.set_zlim(*st["ranges"][2])

    heavy = {}
    for prefix, d in st["domains"].items():
        ag = u.select_atoms(f"resid {d['r0']}:{d['r1']}")
        heavy[prefix] = ag.select_atoms("not element H")

    # statične poti krajišč osi (tip = com + v * length)
    for prefix in st["axis_keys"]:
        d = st["domains"][prefix]
        tips = d["coms"] + d["vecs"] * d["length"]
        ax.plot(tips[:, 0], tips[:, 1], tips[:, 2],
                color=st["colors"][prefix], alpha=0.3, lw=0.8, zorder=1)

    # legenda: en vnos na domeno z osjo
    proxies = []
    for prefix in st["axis_keys"]:
        c = st["colors"][prefix]
        proxies.append(Line2D([0], [0], color=c, lw=2.5, marker="o", ms=7,
                              mfc=c, mec=c,
                              label=f"Domena {prefix[1:]}"))
    if proxies:
        ax.legend(handles=proxies, loc="upper left", fontsize=8, frameon=True)

    time_txt = ax.text2D(0.02, 0.95, "", transform=ax.transAxes, fontsize=12)

    # točke proteina (težki atomi) — posodabljamo s _offsets3d
    scatters = {}
    for prefix, d in st["domains"].items():
        s = ax.scatter([0], [0], [0], s=6, lw=0, color=st["colors"][prefix],
                       alpha=0.55, depthshade=False, zorder=3)
        scatters[prefix] = s

    # puščice osi — Line3DCollection, posodabljamo s set_segments
    arrows = {}
    for prefix in st["axis_keys"]:
        coll = Line3DCollection([np.zeros((2, 3))],
                                colors=st["colors"][prefix],
                                linewidths=3, zorder=5)
        ax.add_collection3d(coll)
        arrows[prefix] = coll

    pngs = []
    time_all = st["time"]
    n_total = len(time_all)
    for k, i in chunk_idxs:
        i = int(i)
        u.trajectory[i]
        prot_com = leaf.center_of_mass()
        for prefix in st["domains"]:
            coords = heavy[prefix].positions - prot_com
            scatters[prefix]._offsets3d = (coords[:, 0], coords[:, 1],
                                           coords[:, 2])
        for prefix in st["axis_keys"]:
            d = st["domains"][prefix]
            arrows[prefix].set_segments(
                axis_segments(d["vecs"][i], d["coms"][i], d["length"]))
        ax.view_init(elev=30, azim=-60 + k * st["rotate"])
        time_txt.set_text(f"t = {time_all[i] / 1000:.1f} ns"
                          f"  ({i} / {n_total})")
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
    ap.add_argument("--pdb", default="atlas_db/PDB_chained/1dd3_A.pdb")
    ap.add_argument("--xtc", default="TEST/1dd3_A_R1.xtc")
    ap.add_argument("--domains", default="d1:1-49,d2:50-128",
                    help="domene: prefix:resid-za-delek, ločene z vejico")
    ap.add_argument("--in-pattern", default="TEST/d*_p1.xvg",
                    help="smeri glavnih osi (samo prva os, p1)")
    ap.add_argument("--out", default="outputs/pai_movie.mp4")
    ap.add_argument("--step", type=int, default=2)
    ap.add_argument("--fps", type=int, default=30)
    ap.add_argument("--workers", type=int, default=0)
    ap.add_argument("--chunk", type=int, default=50)
    ap.add_argument("--codec", default="libopenh264",
                    help="ffmpeg video kodek (brez libx264 v ffmpeg-free)")
    ap.add_argument("--dpi", type=int, default=120)
    ap.add_argument("--rotate", type=float, default=0.0,
                    help="vrtenje kamere [°/frame], 0 = fiksna")
    ap.add_argument("--no-unwrap", action="store_true",
                    help="ne popravljaj predznaka osi")
    args = ap.parse_args()

    axes, time_all = read_axes(args.in_pattern)

    domains = {}
    for item in args.domains.split(","):
        prefix, rr = item.split(":")
        r0s, r1s = rr.split("-")
        domains[prefix] = {"r0": int(r0s), "r1": int(r1s)}
    for p in sorted(domains):
        if p in axes:
            vecs = axes[p] if args.no_unwrap else unwrap_sign(axes[p])
            domains[p]["vecs"] = vecs
    axis_keys = [p for p in sorted(domains) if "vecs" in domains[p]]

    colors = {}
    for p in sorted(domains):
        if p in DOMAIN_STYLE:
            colors[p] = DOMAIN_STYLE[p]["axis"]
        elif p in axis_keys:
            colors[p] = "#333333"
        else:
            colors[p] = "#888888"

    # en prehod čez trajektorijo: obseg, COM domen in dolžina osi
    u = mda.Universe(args.pdb, args.xtc)
    if u.trajectory.n_frames != len(time_all):
        raise SystemExit(f"trajektorija {u.trajectory.n_frames} frame-ov, "
                         f"xvg {len(time_all)}; ne ujemata se")

    mins = np.full(3, np.inf)
    maxs = np.full(3, -np.inf)
    agroups = {}
    heavies = {}
    for p, d in domains.items():
        agroups[p] = u.select_atoms(f"resid {d['r0']}:{d['r1']}")
        heavies[p] = agroups[p].select_atoms("not element H")
        if "vecs" in d:
            d["coms"] = np.empty((len(time_all), 3))
        d["length"] = 0.8 * float(
            np.abs(heavies[p].positions - agroups[p].center_of_mass()).max())
    for i, ts in enumerate(u.trajectory):
        prot_com = u.atoms.center_of_mass()
        c = u.atoms.positions - prot_com
        mins = np.minimum(mins, c.min(0))
        maxs = np.maximum(maxs, c.max(0))
        for p, d in domains.items():
            if "coms" in d:
                d["coms"][i] = agroups[p].center_of_mass() - prot_com

    margin = 0.05 * (maxs - mins).max()
    ranges = [(mins[k] - margin, maxs[k] + margin) for k in range(3)]
    dims = [rng[1] - rng[0] for rng in ranges]
    box = [dim / max(dims) for dim in dims]

    os.makedirs(os.path.dirname(args.out), exist_ok=True)

    frame_tasks = [(k, i) for k, i
                   in enumerate(np.arange(0, u.trajectory.n_frames, args.step))]
    chunks = [frame_tasks[s:s + args.chunk]
              for s in range(0, len(frame_tasks), args.chunk)]
    n_frames = len(frame_tasks)

    state = {
        "pdb": args.pdb, "xtc": args.xtc, "time": time_all,
        "domains": domains, "axis_keys": axis_keys, "colors": colors,
        "ranges": ranges, "box": box, "dpi": args.dpi, "rotate": args.rotate,
    }
    n_workers = args.workers or min(os.cpu_count() or 1, len(chunks))
    proc = _build_pipeline(args.out, args.codec, args.fps)

    print(f"{n_frames} frame-ov, {len(chunks)} kupčkov, {n_workers} procesov",
          flush=True)
    done = 0
    if n_workers > 1:
        with ProcessPoolExecutor(max_workers=n_workers,
                                 initializer=_init_worker,
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