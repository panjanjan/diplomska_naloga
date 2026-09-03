#!/usr/bin/fish
# prenese datoteke iz ATLAS baze

pushd "$ROOT/atlas_db"

argparse t/test -- $argv
or exit

# število paralelnih procesov za xargs
set n 5

# seznam proteinov
set proteins "2024_11_18_ATLAS_pdb.txt"

set data_dir "analysis"
set test_dir "TEST"

set url "https://www.dsimb.inserm.fr/ATLAS/api/ATLAS/analysis"
set header "accept: application/octet-stream"

function run -a target_dir prot_list
  echo "storing into $(pwd)/$target_dir"
  echo "$prot_list"

  test -d "$target_dir" || mkdir -p "$target_dir"

  echo "$prot_list" |\
    xargs -P "$n" -I {} \
    curl \
      -X "GET" "$url/{}" \
      -H "$header" \
      -o "$target_dir/{}.zip"
end

if set -ql _flag_test
  # uporabi par proteinov za test
  run "$test_dir" "$(head -n 5 "$proteins")"
else
  run "$data_dir" "$(cat "$prot_list")"
end

popd
