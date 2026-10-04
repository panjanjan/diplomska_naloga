"""Angles between domain inertia axes (MDAnalysis version).

Python counterpart of principal_axes_of_inertia.r: writes
outputs/PAI_2/{protein}_angles.csv with columns
frame,R1_p1..R3_p3, where Rn_pm is the angle between the first
principal axis of domain 1 and the m-th principal axis of
domain 2 in replicate n.

Domain segments come from outputs/two_domains.csv (a domain may
span several disjoint residue ranges). Proteins run in parallel,
one worker per protein, like R mclapply().
"""

import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from csv import DictReader, writer
from pathlib import Path

import MDAnalysis as mda
import numpy as np

ROOT = os.getenv("ROOT")
if not isinstance(ROOT, str):
    raise TypeError("ROOT environment variable not set")

PDB_DIR = Path(ROOT, "atlas_db", "PDB_chained")
TRAJ_DIR = Path(ROOT, "atlas_db", "trajectories")
DOMAINS_CSV = Path(ROOT, "outputs", "two_domains.csv")
OUT_DIR = Path(ROOT, "outputs", "PAI_2")

# like R: min(detectCores() - 1, 10)
N_WORKERS = max(min((os.cpu_count() or 2) - 1, 10), 1)


def load_domains(path):
    """Map {protein: {1: [(start, end)], 2: [...]}} from CSV rows."""
    domains = {}
    with open(path, newline="") as fh:
        for row in DictReader(fh):
            segs = domains.setdefault(row["protein"], {}).setdefault(
                int(row["domain"]), []
            )
            segs.append((int(row["start"]), int(row["end"])))
    for segmap in domains.values():
        for segs in segmap.values():
            segs.sort()
    return domains


def domain_sel(segments):
    """CA selection covering every (start, end) segment of a domain."""
    inner = " or ".join(f"resid {s}:{e}" for s, e in segments)
    return f"protein and name CA and ({inner})"


def principal_axes(inertia):
    """Principal axes as columns, ordered like R eigen() (decreasing)."""
    vals, vecs = np.linalg.eig(inertia)
    order = np.argsort(vals.real)[::-1]
    return vecs[:, order].real


def angle(d1p, d2p):
    """Angle in [0, pi/2]; abs() fixes eigenvector sign ambiguity."""
    cos = abs(float(np.dot(d1p, d2p)))
    cos /= float(np.linalg.norm(d1p) * np.linalg.norm(d2p))
    return float(np.arccos(np.clip(cos, -1.0, 1.0)))


def run_replicate(pdb_file, traj_file, d1_sel, d2_sel):
    """d1p1-d2p{1,2,3} angles for every frame of one trajectory."""
    print(traj_file, flush=True)
    u = mda.Universe(pdb_file, traj_file)
    d1_atoms = u.select_atoms(d1_sel)
    d2_atoms = u.select_atoms(d2_sel)
    if d1_atoms.n_atoms == 0 or d2_atoms.n_atoms == 0:
        raise ValueError(f"empty domain selection in {traj_file}")

    angles = np.empty((u.trajectory.n_frames, 3))
    for i, _ts in enumerate(u.trajectory):
        d1_axes = principal_axes(d1_atoms.moment_of_inertia())
        d2_axes = principal_axes(d2_atoms.moment_of_inertia())
        for j in range(3):
            angles[i, j] = angle(d1_axes[:, 0], d2_axes[:, j])
    return angles


def pad_replicates(replicates):
    """Repeat the last row of shorter replicates, like the R script."""
    n_frames = max(map(len, replicates))
    padded = []
    for rep in replicates:
        if len(rep) < n_frames:
            tail = np.repeat(rep[-1:], n_frames - len(rep), axis=0)
            rep = np.vstack([rep, tail])
        padded.append(rep)
    return padded


def write_csv(out, replicates):
    """Angles to {protein}_angles.csv; 1-based frame column."""
    header = ["frame"]
    for n in (1, 2, 3):
        header += [f"R{n}_p{j}" for j in (1, 2, 3)]
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as fh:
        csv = writer(fh)
        csv.writerow(header)
        for i in range(len(replicates[0])):
            csv.writerow(
                [i + 1] + [rep[i, j] for rep in replicates for j in range(3)]
            )


def run_protein(protein, segmap):
    """Full pipeline for one protein; returns the output path."""
    pdb_file = PDB_DIR / f"{protein}.pdb"
    if not pdb_file.is_file():
        raise FileNotFoundError(f"no PDB for {protein}")
    traj_files = sorted(
        f for f in TRAJ_DIR.glob(f"{protein}*.xtc")
        if f.name.startswith(protein + "_")
    )
    if len(traj_files) != 3:
        raise FileNotFoundError(
            f"expected 3 trajectories for {protein}, "
            f"found {len(traj_files)}"
        )
    if sorted(segmap) != [1, 2]:
        raise ValueError(f"expected domains 1 and 2 for {protein}")
    d1_sel = domain_sel(segmap[1])
    d2_sel = domain_sel(segmap[2])

    replicates = pad_replicates(
        [
            run_replicate(str(pdb_file), str(traj), d1_sel, d2_sel)
            for traj in traj_files
        ]
    )
    out = OUT_DIR / f"{protein}_angles.csv"
    write_csv(out, replicates)
    return str(out)


def main():
    """One worker per protein; a failed protein logs and continues."""
    domains = load_domains(DOMAINS_CSV)
    proteins = sorted(domains)
    print(f"{len(proteins)} proteins, {N_WORKERS} workers", flush=True)
    with ProcessPoolExecutor(max_workers=N_WORKERS) as ex:
        pending = {
            ex.submit(run_protein, p, domains[p]): p for p in proteins
        }
        for fut in as_completed(pending):
            protein = pending[fut]
            try:
                out = fut.result()
            except Exception as err:
                print(f"{protein}: FAILED: {err}", flush=True)
            else:
                print(f"{protein}: written: {out}", flush=True)


if __name__ == "__main__":
    main()
