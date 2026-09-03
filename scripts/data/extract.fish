#!/bin/fish
# iz vseh zip datotek proteinov skopira neke željene datoteke

pushd "$ROOT/atlas_db"

argparse query target -- $argv
or exit

if ! set -ql _flag_query || ! set -ql _flag_target
  echo "usage: extract.fish <query> <target>"
  echo ""
  echo "query: pattern to match files using 'find'"
  echo "target: output directory for extracted files"
  exit
else
  set query "$argv[1]"
  set target "$argv[2]"
end

test -d "$target" || mkdir -p "$target"

# sem unzipa datoteke, skopira ven željene in izbriše nepotrebne
rm -rf tmp/*

function process_zip -a name
  set base (path basename --no-extension "$name")

  # preveri če datoteka že obstaja
  ls "$target" | grep -q "$base" && return

  unzip -qd "tmp/$base" "$name"

  # absoluten path do datoteke za mv
  set fname (find "tmp/$base" -name "$query")
  mv "$fname" "$target"

  # počisti
  rm -r "tmp/$base"
end

# izogibam sem paralelizacije tega, ker se lahko hitro zafila prostor
for zipf in (find analysis -name "*.zip")
  process_zip "$zipf"
end

popd
