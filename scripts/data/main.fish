#!/bin/bash
mkdir -p "$ROOT/atlas_db"

# ./calculate_atlas_size.fish
"$ROOT/scripts/data/download_atlas.fish"

# TODO: preveri če dela
"$ROOT/scripts/data/extract.fish" --query "*_RMSF.tsv" --target "RMSF"
"$ROOT/scripts/data/extract.fish" --query "*.pdb"      --target "PDB"
"$ROOT/scripts/data/extract.fish" --query "*.tpr"      --target "trajectories"
"$ROOT/scripts/data/extract.fish" --query "*.xtc"      --target "trajectories"

"$ROOT/scripts/data/mdconvert_xtc.fish"
