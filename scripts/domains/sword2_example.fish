#!/bin/fish
# FIX!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# FIX!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# FIX!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# FIX!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
set protein "1a62.pdb"

# dobi ID od trenutnega proteina
set protein_id (basename "$protein" .pdb)
echo "$protein_id"
#> 1a62

# najdi vnos proteina v seznamu baze
set db_entry (grep "$protein_id" "$ROOT/"atlas_db/2024_11_18_ATLAS_pdb.txt)
echo "$db_entry"
#> 1a62_A

# dobi uporabljeno verigo
set chain (string split "_" "$db_entry" -f 2)
echo "$chain"
#> A

conda activate "$CONDA_ENV_NAME"
"$SWORD_PATH" -i "$DB/$protein" -o "$SWO" -c "$chain"
conda deactivate