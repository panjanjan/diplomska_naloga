#!/bin/python3
import csv
import glob
import json
import os
from pathlib import Path
from typing import TypedDict, cast

"""
Format podatkov ---------------------------------------------------------------------------------------
zgradil bo CSV iz JSON datotek.

`sword_results.csv` vsebuje podatke iz JSON datotek, ki jih ustvari SWORD2.
Vrstice se združujejo po proteinih in proteini po svojih particijah.

protein:   PDB koda ter veriga, ki je bila uporabljena za SWORD2
aindex:    ambiguity index proteina
partition: indeks particije. Optimalna ima 0, alternativne 1 ali več
quality:   ocena particije
domain:    indeks domene. Prva domena 1, druga domena 2, ...
AUL:       AUL vrednost domene, ponovljena na vsaki vrstici PU-ja
start:     prva aminokislina proteinske enote (PU)
end:       zadnja aminokislina proteinske enote (PU)

Ena vrstica na proteinsko enoto (PU), ne na domeno. Domeno lahko sestavlja
več PU-jev, ki si ne sledijo nujno zvezno vzdolž zaporedja.

Primer za 1a62_A:

protein aindex partition quality domain AUL start end
1a62_A  1      0         0       1      81  1     130 <---| opt.
1a62_A  1      1         1       1      70  1     47 <----| alt. 1
1a62_A  1      1         1       2      0   48    94      |
1a62_A  1      1         1       3      46  95    130     |
1a62_A  1      2         3       1      72  1     47 <----| alt. 2
1a62_A  1      2         3       2      8   48    130     |
1a62_A  1      3         1       1      76  1     130 <---| alt. 3
1a62_A  1      3         1       2      0   48    94      |
"""

ROOT = os.getenv("ROOT")
if ROOT is None:
    raise RuntimeError("ROOT environment variable is not set")

field_names: list[str] = [
    "protein",
    "aindex",
    "partition",
    "quality",
    "domain",
    "AUL",
    "start",
    "end",
]

csv_name: Path = Path(ROOT, "outputs", "sword_results.csv")


class PU(TypedDict):
    AUL: int


class Domain(TypedDict):
    AUL: int
    PUs: dict[str, PU]


class Partition(TypedDict):
    Quality: str
    Domains: dict[str, Domain]


def process_partition(part: Partition) -> list[list[int]]:
    quality = len(part["Quality"])
    dom_list: list[list[int]] = []

    for i, domain in enumerate(part["Domains"].values()):
        aul: int = domain["AUL"]
        # vsak PU dobi svojo vrstico, saj domene niso nujno zvezne
        for pu in domain["PUs"]:
            pu_s, pu_e = pu.split("-")
            dom_list.append([quality, i + 1, aul, int(pu_s), int(pu_e)])

    return dom_list


def process_report(fname: str) -> list[list[str | int]]:
    # globalni podatki
    pname: str = fname.removeprefix("outputs/sword_output/").removesuffix("/SWORD2_summary.json")[:6]

    with open(fname) as fp:
        parsed: object = json.load(fp)

    if not isinstance(parsed, dict):
        raise TypeError(f"expected a JSON object in {fname}")
    data: dict[str, object] = parsed

    aidx = data["Ambiguity index"]
    if not isinstance(aidx, str):
        raise ValueError(f"Ambiguity index is not a string in {fname}")

    # pojdi čez vsako particijo in izvleci podatke
    part_reports: list[list[str | int]] = []
    pid: int = 0
    for key in data:
        if "partition" in key:
            raw = data[key]
            if not isinstance(raw, dict):
                raise TypeError(f"partition {key} is not an object in {fname}")
            part: Partition = cast(Partition, raw)
            # vrne N seznamov, ki predstavljajo podatke o domenah
            dom_list = process_partition(part)
            for items in dom_list:
                part_reports.append([pname, len(aidx), pid, *items])
            pid += 1

    return part_reports


def main() -> None:
    with open(csv_name, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(field_names)
        for fname in glob.glob("outputs/sword_output/*/*/*.json"):
            writer.writerows(process_report(fname))
    print("written to", csv_name)


if __name__ == "__main__":
    main()
