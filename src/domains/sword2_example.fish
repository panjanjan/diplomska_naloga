#!/bin/fish
set test_dir "TEST_SWORD"
test -d "$test_dir" || mkdir "$test_dir"
echo "output: $test_dir"

set protein "$ROOT/atlas_db/PDB_chained/1a62_A.pdb"
echo "protein: $protein"

# dobi uporabljeno verigo
set chain (path basename --no-extension $protein | cut -d '_' -f2)
echo "chain: $chain"
echo ""

conda activate $SWORD_CONDA_ENV

$SWORD_PATH -i $protein -o $test_dir -c $chain

conda deactivate
