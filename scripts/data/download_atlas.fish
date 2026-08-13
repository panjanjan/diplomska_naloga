#!/usr/bin/fish
# prenese datoteke iz ATLAS baze

pushd "$ROOT/atlas_db"

set n 5
set prot_list "2024_11_18_ATLAS_pdb.txt"
set url "https://www.dsimb.inserm.fr/ATLAS/api/ATLAS/analysis"
set header "accept: application/octet-stream"
set data_dir "analysis"

test -d $data_dir || mkdir -p $data_dir

# -P: število istočasnih procesov
# -I: placeholder za input
cat "$prot_list" | xargs -P "$n" -I {} curl -X "GET" "$url/{}" -H "$header" -o "$data_dir/{}.zip"

popd
