# \<naslov diplomske naloge>

Repository vsebuje vso kodo, analize, izračune, s katerimi sem raziskoval proteine iz
[ATLAS MD](https://www.dsimb.inserm.fr/ATLAS/about.html) za diplomsko nalogo.

# Struktura

- [`atlas_db`](./atlas_db/): podatki iz baze.
- [`outputs`](./outputs/): rezultati analiz.
- [`scripts`](./scripts/): vsebuje svoj [README](./scripts/README.md)
- [`.envrc`](./.envrc): globalne spremenljivke potrebne za vse skripte. Priporočam direnv, da se vedno sourca-jo.

# Setup

Testirano na Fedora Linux 44 (Workstation Edition).

- Python 3.14.6: <https://www.python.org/downloads/release/python-3146/>
    - SWORD2: <https://github.com/DSIMB/SWORD2>.
    - mdconvert: <https://mdtraj.org/1.9.4/mdconvert.html>
- R 4.6.1: <https://www.r-project.org/>
    - bio3d 2.4-5 <https://thegrantlab.org/bio3d/>
    - dplyr 1.2.1: <https://dplyr.tidyverse.org//>
    - ggplot2 4.0.3: <https://ggplot2.tidyverse.org/>
    - patchwork 1.3.2: <https://patchwork.data-imaginist.com/index.html> - NI POTREBNO
    - plotly 4.12.0: <https://plotly.com/r/> - NI POTREBNO
- GROMACS 2026.3: <https://manual.gromacs.org/current/install-guide/index.html> - NI POTREBNO

# Pipeline

Glej [scripts/README](./scripts/README.md).
