#!/usr/bin/env python3
"""Interaktivni HTML z dinamiko proteina in prve glavne osi domen.

Enako kot pai_movie.py (MP4), a izhod je samostojen HTML:
  * slider + play/pause
  * prizor lahko z miško rotiraš in povečaš

Frame-i se vzorčijo po --stride (privzeto vsak 10. -> 1 ns koraki).
"""
# ------------------------------------------------------------------------------
import argparse
import glob
import os

import MDAnalysis as mda
import numpy as np
import plotly.graph_objects as go

DOMAIN_STYLE = {
    "d1": {"axis": "#164e7d", "points": "#7da6c9"},
    "d2": {"axis": "#c0392b", "points": "#d9a4a4"},
}


def read_axes(pattern):
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", default="atlas_db/PDB_chained/1dd3_A.pdb")
    ap.add_argument("--xtc", default="TEST/1dd3_A_R1.xtc")
    ap.add_argument("--domains", default="d1:1-49,d2:50-128",
                    help="domene: prefix:resid-za-delek, ločene z vejico")
    ap.add_argument("--in-pattern", default="TEST/d*_p1.xvg",
                    help="smeri glavnih osi (samo p1)")
    ap.add_argument("--out", default="outputs/pai_movie.html")
    ap.add_argument("--stride", type=int, default=10,
                    help="vsak n-ti frame gre v animacijo")
    ap.add_argument("--frame-duration", type=int, default=150,
                    help="ms na frame pri predvajanju")
    ap.add_argument("--heavy", action="store_true",
                    help="namesto Cα vse težke atome (večji HTML)")
    args = ap.parse_args()

    axes, time_all = read_axes(args.in_pattern)

    d = {}
    for item in args.domains.split(","):
        prefix, rr = item.split(":")
        r0s, r1s = rr.split("-")
        d[prefix] = {"r0": int(r0s), "r1": int(r1s),
                     "vec": unwrap_sign(axes[prefix])}
    keys = sorted(d)

    u = mda.Universe(args.pdb, args.xtc)
    if u.trajectory.n_frames != len(time_all):
        raise SystemExit(f"trajektorija {u.trajectory.n_frames} frame-ov, "
                         f"xvg {len(time_all)}; ne ujemata se")

    ca = {}
    for p in keys:
        ag = u.select_atoms(f"resid {d[p]['r0']}:{d[p]['r1']}")
        if args.heavy:
            ca[p] = ag.select_atoms("not element H")
        else:
            ca[p] = ag.select_atoms("name CA")
        # dolžina osi = 0.8 × obseg domene okoli COM
        ref = ag.select_atoms("not element H")
        d[p]["length"] = 0.8 * float(
            np.abs(ref.positions - ag.center_of_mass()).max())

    idxs = np.arange(0, u.trajectory.n_frames, args.stride)
    n = len(idxs)
    print(f"{u.trajectory.n_frames} frame-ov -> {n} v animaciji", flush=True)

    pts = {p: np.empty((n, ca[p].n_atoms, 3)) for p in keys}
    coms = {p: np.empty((n, 3)) for p in keys}
    for j, i in enumerate(idxs):
        u.trajectory[int(i)]
        prot_com = u.atoms.center_of_mass()
        for p in keys:
            pts[p][j] = ca[p].positions - prot_com
            coms[p][j] = ca[p].center_of_mass() - prot_com

    allpts = np.vstack([pts[p].reshape(-1, 3) for p in keys])
    D = float(np.abs(allpts).max()) * 1.15

    # --- trace-i ------------------------------------------------------------
    fig = go.Figure()

    def trail_trace(p):
        c = DOMAIN_STYLE[p]["axis"]
        tips = coms[p] + d[p]["vec"][idxs] * d[p]["length"] / 2
        return go.Scatter3d(
            x=tips[:, 0], y=tips[:, 1], z=tips[:, 2], mode="lines",
            line=dict(color=c, width=1.5), opacity=0.35,
            hoverinfo="skip", showlegend=False)

    def scatter(p, j):
        c = DOMAIN_STYLE[p]["points"]
        return go.Scatter3d(
            x=pts[p][j, :, 0], y=pts[p][j, :, 1], z=pts[p][j, :, 2],
            mode="markers", marker=dict(size=3, color=c, opacity=0.6),
            name=f"Domena {p[1:]}", hoverinfo="skip")

    def axis_line(p, j):
        c = DOMAIN_STYLE[p]["axis"]
        v = d[p]["vec"][idxs[j]]
        L = d[p]["length"] / 2
        com = coms[p][j]
        return go.Scatter3d(
            x=[com[0] - v[0] * L, com[0] + v[0] * L],
            y=[com[1] - v[1] * L, com[1] + v[1] * L],
            z=[com[2] - v[2] * L, com[2] + v[2] * L],
            mode="lines", line=dict(color=c, width=6),
            name=f"Prva os domena {p[1:]}", hoverinfo="skip")

    # vsi frame-i in začetno stanje imajo ENAKO število trace-ov, da se
    # animacija ujema po indeksih: [trail d1, trail d2, scatter d1, scatter d2,
    #                                axis d1, axis d2]
    def frame_data(j):
        return ([trail_trace(p) for p in keys]
                + [scatter(p, j) for p in keys]
                + [axis_line(p, j) for p in keys])

    for tr in frame_data(0):
        fig.add_trace(tr)

    def frame_layout(j):
        return go.Layout(annotations=[dict(
            text=f"t = {time_all[idxs[j]] / 1000:.1f} ns",
            x=0.02, y=0.98, xref="paper", yref="paper",
            showarrow=False, font=dict(size=14))])

    frames = [go.Frame(
        data=frame_data(j), name=str(j), layout=frame_layout(j))
        for j in range(n)]
    fig.frames = frames

    # --- slider + play/pause ------------------------------------------------
    slider_steps = [dict(
        method="animate",
        args=[[str(j)], dict(mode="immediate",
                             frame=dict(duration=0, redraw=True),
                             transition=dict(duration=0))],
        label=f"{(j * args.stride) / 10:.0f}") for j in range(n)]
    sliders = [dict(active=0,
                    currentvalue=dict(prefix="t = ", font=dict(size=12),
                                      xanchor="left"),
                    steps=slider_steps)]

    buttons = [
        dict(label="▶ Predvajaj", method="animate",
             args=[None, dict(frame=dict(duration=args.frame_duration,
                                         redraw=True),
                              fromcurrent=True,
                              transition=dict(duration=0))]),
        dict(label="⏸", method="animate",
             args=[[None], dict(mode="immediate",
                                frame=dict(duration=0, redraw=False))]),
    ]

    fig.layout.update(
        title="Dinamika proteina in prva os vztrajnosti domen",
        scene=dict(
            xaxis=dict(range=[-D, D], title="x [Å]"),
            yaxis=dict(range=[-D, D], title="y [Å]"),
            zaxis=dict(range=[-D, D], title="z [Å]"),
            aspectmode="cube"),
        updatemenus=[dict(type="buttons", buttons=buttons,
                          x=0.0, y=1.08, xanchor="left", yanchor="top")],
        sliders=sliders)

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    fig.write_html(args.out, include_plotlyjs="inline")
    print(f"zapisano: {args.out} "
          f"({os.path.getsize(args.out) / 1e6:.1f} MB)")


if __name__ == "__main__":
    main()