#!/bin/fish
# iz vseh zip datotek proteinov skopira neke željene datoteke

pushd "$ROOT/atlas_db"

# argument parsing
argparse q/query t/target T/test -- $argv;
or exit

if ! set -ql _flag_query; or ! set -ql _flag_target
  echo "usage: extract.fish -q <str> -t <str> [ -T/--test ]"
  exit
else
  set query "$argv[1]"
  set target "$argv[2]"
end

echo "using query: $query, target: $(find $ROOT -name $target)"
while read --nchars 1 -l response --prompt-str="Proceed? (y/n): "; or return 1
  switch $response
    case "y" "Y"
      break
    case "n" "N"
      exit
    case '*'
      echo "invalid input"
      continue
  end
end

# ustvari target directory če še ne obstaja
test -d "$target"; or mkdir "$target"

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

  # samo ime proteina z verigo, npr. "1dd3_A"
  set base (path basename --no-extension "$zipf")

  # preskoči tiste, ki že obstajajo
  if count $target/$base* > /dev/null
    echo "(skip)"
    continue
  end

  # extractaj v začasni directory
  # q : quiet
  # d : directory
  # n : no overwriting
  unzip -qnd "tmp/$base" "$zipf"

  # absoluten path do željene datoteke za mv
  set query_files (find "tmp/$base" -name "$query")

  # prestavi pomembne datoteke, zbriši ostanek
  # unqoatano ker je lahko več datotek skupaj
  mv $query_files "$target"
  rm -r "tmp/$base"
  echo "done"
end

popd
