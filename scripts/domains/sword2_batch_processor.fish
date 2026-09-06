#!/bin/fish
pushd "$ROOT"

# uporabi 3 datoteke namesto vseh
argparse t/test -- $argv

# error messages proteinov, ki so failali pri analizi
set log "$SWO/failed_proteins.log"
set tmp_log "$SWO/tmp.log"

# seznam proteinov, ki so uspešno prešli analizo. Na začetku je prazen.
# SWORD2 lahko fail-a po tem, ko že ustvari directory, ampak preden naredi karkoli
set processed "$SWO/processed_proteins.log"

echo "using the following parameters:

* inputs directory      $DB
* outputs directory     $SWO
* SWORD2 path           $SWORD_PATH
* conda environment     $SWORD_CONDA_ENV

using the following log files:

* failed proteins       $log
* successfull proteins  $processed
"

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

echo
echo "==> initializing log files"
# ustvari samo ko še ne obstajajo
# resetiraj tmp log v vsakem primeru
touch "$log"
touch "$tmp_log"
touch "$processed"
echo -n > "$tmp_log"

echo
echo "==> checking past analysis integrity ($processed)"
# preveri, da ima res vsak protein iz seznama svoj output directory
# ls fail-a če zvezdica ne najde directorija - take odstrani s seznama
for protein in (cat "$processed")
    echo -ne "\r$protein: "

    test -e "$SWO/$protein"*

    if test $status -ne 0
        echo "no output directory. removing from list"
    else
        echo "$protein" >> "$tmp_log"
        echo -n "pass"
    end
end

# NOTE: tmp se izbriše tukaj
mv "$tmp_log" "$processed"

echo
echo "==> obtaining list of PDB files"
# find vrne space-separated list of values
# vzame subset datotek če testiraš
set files (find "$DB" -name "*.pdb")

if set -ql _flag_test
    echo "TEST MODE: using 3 files"
    set files (string split " " $files | head -n 3)
end

set n (count $files)
set i 0
set failed 0

echo
echo "==> activating conda environment"
conda activate "$SWORD_CONDA_ENV"

echo
echo "==> running analyses"

for protein in $files
    set i          (math $i + 1)
    set prot_fname (path basename --no-extension "$protein")
    set chain      (string split "_" "$prot_fname" -f 2)

    echo -n "[$i/$n] $protein: "

    # če je protein že obdelan, pojdi na naslednjega
    grep -q "$prot_fname" "$processed" &&\
        echo "skip" && continue ||\
        echo -n "... "

    # shrani stderr v tmp datoteko
    set current_tmp "$SWO/tmp_$i.log"
    touch "$current_tmp"

    "$SWORD_PATH" \
        -i "$protein" \
        -c "$chain" \
        -o "$SWO/$prot_fname" > /dev/null 2> "$current_tmp"

    if test $status -ne 0
        set failed (math $failed + 1)
        echo "error. status: $status"

        # v primeru, da gre nekaj narobe, shrani v log
        echo -e "--- $(date "+%Y-%m-%d %H:%M") $protein ---\n" >> "$log"
        cat "$current_tmp"                                     >> "$log"
        echo ""                                                >> "$log"
    else
        # v primeru, da gre vse prav, dodaj protein na seznam
        echo "done"
        echo "$prot_fname" >> "$processed"
    end

    rm "$current_tmp"
end

echo
echo "$failed/$n proteins failed"
echo "logs written to $log"
echo "results in $SWO"

# počisti za sabo
conda deactivate
popd
