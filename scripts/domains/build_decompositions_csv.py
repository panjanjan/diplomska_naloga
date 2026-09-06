#!/bin/python3
import csv
import glob
import json
import os
from pathlib import Path

# zgradil bo CSV iz JSON datotek. vrstice se grupirajo po proteinih in
# particijah optimalna particija ima vrednost 0. Domain predstavlja
# število zaporedne domene. 1 = prva domena. range start-end domene je glede
# na prvo aminokislino prve PU in zadnjo aminokislino zadnje PU. AUL je od
# domene. A-index in Qualtiy sta 0, če je prazen string, drugače je število
# zvezdic '*'.
#
# primer za 1a62_A:
#
#   protein aindex partition quality domain AUL start end
#   1a62_A  1      0         0       1      81  1     130<---| optimal
#   1a62_A  1      1         1       1      70  1     47<----| alt. 1
#   1a62_A  1      1         1       2      0   48    94     |
#   1a62_A  1      1         1       3      46  95    130    |
#   1a62_A  1      2         3       1      72  1     47<----| alt. 2
#   1a62_A  1      2         3       2      8   48    130    |
#   1a62_A  1      3         1       1      76  1     130<---| alt. 3
#   1a62_A  1      3         1       2      0   48    94     |
#
ROOT = os.getenv("ROOT")

field_names = [
    "protein",
    "aindex",
    "partition",
    "quality",
    "domain",
    "AUL",
    "start",
    "end",
]

csv_name = Path(ROOT, "outputs", "sword_results.csv")


def domain_bounds(domain) -> tuple[int, int]:
    begin = 0
    end = 0

    for i, pu in enumerate(domain["PUs"]):
        pu_s, pu_e = pu.split("-")
        pu_s, pu_e = int(pu_s), int(pu_e)
        # če je prva, vzami start index
        if i == 0:
            begin = pu_s
        # najdi zadnjo aminokislino
        if pu_e > end:
            end = pu_e

    return begin, end


def process_partition(part) -> list:
    quality = len(part["Quality"])
    domains = part["Domains"]
    dom_list = []

    for i, domain in enumerate(domains.values()):
        aul = domain["AUL"]
        begin, end = domain_bounds(domain)
        dom_list.append([quality, i + 1, aul, begin, end])

    return dom_list


def process_report(fname: str) -> list:
    # globalni podatki
    pname = fname.removeprefix("outputs/sword_output/").removesuffix("/SWORD2_summary.json")[:6]
    fp = open(fname, "r")
    obj = json.load(fp)
    aidx = obj["Ambiguity index"]

    # pojdi čez vsako particijo in izvleci podatke
    part_reports = []
    pid = 0
    for key in obj:
        if key.find("partition") != -1:
            # vrne N seznamov, ki predstavljajo podatke o domenah
            dom_list = process_partition(obj[key])
            for items in dom_list:
                part_reports.append([pname, len(aidx), pid, *items])
            pid += 1

    fp.close()
    return part_reports


if __name__ == "__main__":
    with open(csv_name, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(field_names)
        for fname in glob.glob("outputs/sword_output/*/*/*.json"):
            report = process_report(fname)
            for line in report:
                writer.writerow(line)
    print("written to", csv_name)
