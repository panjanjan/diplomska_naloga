#!/bin/fish
# iz vseh zip datotek proteinov skopira neke željene datoteke

pushd "$ROOT/atlas_db"

argparse q/query t/target T/test -- $argv;
or exit

if ! set -ql _flag_query; or ! set -ql _flag_target
  echo "usage: extract.fish -q <str> -t <str> [ -T/--test ]"
  exit
else
  set query "$argv[1]"
  set target "$argv[2]"
end

test -d "$target"; or mkdir -p "$target"

# sem unzipa datoteke, skopira ven željene in izbriše nepotrebne
rm -rf "tmp/*"
mkdir -p tmp

set files (ls analysis/*.zip)

# vzami subset za testiranje
if set -ql _flag_test
  set files (string split " " $files | head -n 5)
end

set n (count $files)
set i 0

# izogibam sem paralelizacije tega, ker se lahko hitro zafila prostor
for zipf in $files
  set i (math $i + 1)
  echo -n "[$i/$n] $zipf ... "

  set base (path basename --no-extension "$zipf")

  # preskoči tiste, ki že obstajajo
  if count $target/$base* > /dev/null
    echo "(skip)"
    continue
  end

  unzip -qd "tmp/$base" "$zipf"

  # absoluten path do željene datoteke za mv
  set query_file (find "tmp/$base" -name "$query")

  mv "$query_file" "$target"

  # počisti
  rm -r "tmp/$base"
  echo "done"
end

popd
