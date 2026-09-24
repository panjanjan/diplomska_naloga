#!/usr/bin/env python3
"""Generira GROMACS indeksne datoteke (ndx) za domene iz two_domains.csv.

Vsak protein iz PDB_chained/{protein}.pdb dobi grupe:
  [ System ]  vsi atomi
  [ D1 ]      atomi domen z oznako domain=1
  [ D2 ]      atomi domen z oznako domain=2

Domena je lahko sestavljena iz več nizov (start..end). Poženi pred
com_distances_gmx.fish: prvi zgradi indekse, drugi jih samo uporabi.
"""
import argparse
import csv
import os
import sys


def read_pdb(pdb):
    """Vrni (seriali, resi) za ATOM vrstice: s[6:11], r[22:26]. Ujema awk."""
    serials = []
    resids = []
    with open(pdb) as fh:
        for line in fh:
            if line[:4] == "ATOM":
                serials.append(int(line[6:11]))
                resids.append(int(line[22:26]))
    return serials, resids


def fmt15(nums):
    """Zapiši številke po 15 na vrstico (enako kot fish fmt15)."""
    lines = []
    for i in range(0, len(nums), 15):
        lines.append(" ".join(str(n) for n in nums[i:i + 15]))
    return "\n".join(lines) + "\n"


def build_index(protein, pdb, ranges, out_dir, force):
    out = os.path.join(out_dir, f"{protein}.ndx")
    if not force and os.path.exists(out):
        print(f"{protein}: skip (exists)")
        return True

    serials, resids = read_pdb(pdb)
    domains = {1: set(), 2: set()}
    for domain, start, end in ranges:
        if domain not in domains:
            print(f"{protein}: unexpected domain {domain}", file=sys.stderr)
            return False
        for s, r in zip(serials, resids):
            if start <= r <= end:
                domains[domain].add(s)

    if not domains[1] or not domains[2]:
        print(f"{protein}: empty domain index", file=sys.stderr)
        return False

    all_atoms = sorted(set(serials))
    with open(out, "w") as fh:
        fh.write("[ System ]\n")
        fh.write(fmt15(all_atoms))
        fh.write("[ D1 ]\n")
        fh.write(fmt15(sorted(domains[1])))
        fh.write("[ D2 ]\n")
        fh.write(fmt15(sorted(domains[2])))
    print(f"{protein}: ok -> {out}")
    return True


def main():
    root = os.environ.get("ROOT", os.getcwd())
    ap = argparse.ArgumentParser(
        description="Indeksne datoteke domen za com_distances_gmx.fish")
    ap.add_argument("--domains", default=os.path.join(root, "outputs",
                                                      "two_domains.csv"),
                    help="CSV z domenami (protein,domain,start,end)")
    ap.add_argument("--pdb-dir", default=os.path.join(root, "atlas_db",
                                                      "PDB_chained"))
    ap.add_argument("--out-dir", default=os.path.join(root, "outputs",
                                                      "COM_gmx", "index"))
    ap.add_argument("--proteins", nargs="*",
                    help="samo ti proteini; privzeto vsi iz CSV")
    ap.add_argument("--force", action="store_true",
                    help="prepiši obstoječe ndx")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    by_protein = {}
    with open(args.domains) as fh:
        for row in csv.DictReader(fh):
            by_protein.setdefault(row["protein"], []).append(
                (int(row["domain"]), int(row["start"]), int(row["end"])))

    proteins = args.proteins or sorted(by_protein)
    n_ok, n_fail = 0, 0
    for protein in proteins:
        if protein not in by_protein:
            print(f"{protein}: not in {args.domains}", file=sys.stderr)
            n_fail += 1
            continue
        pdb = os.path.join(args.pdb_dir, f"{protein}.pdb")
        if not os.path.exists(pdb):
            print(f"{protein}: missing {pdb}", file=sys.stderr)
            n_fail += 1
            continue
        if build_index(protein, pdb, by_protein[protein],
                       args.out_dir, args.force):
            n_ok += 1
        else:
            n_fail += 1

    print(f"done: {n_ok} ok, {n_fail} failed")
    sys.exit(1 if n_fail else 0)


if __name__ == "__main__":
    main()