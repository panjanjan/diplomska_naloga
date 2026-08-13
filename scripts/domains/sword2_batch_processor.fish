#!/bin/fish
# FIX!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# FIX!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# FIX!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# FIX!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

pushd "$ROOT"
set out "$ROOT/outputs"

# error messages proteinov, ki so failali pri analizi
set log "$out/failed_proteins.log"
set tmp_log "$out/tmp.log"

# seznam proteinov, ki so uspešno prešli analizo
# SWORD2 lahko fail-a po tem, ko že ustvari directory,
# zato je potreben bolj robusten sistem kot le "ls $SWO | grep protein"
set processed "$out/processed_proteins.log"

echo "using the following parameters:
* data directory        $DB
* output directory      $SWO
* SWORD2 path           $SWORD_PATH
* conda environment     $CONDA_ENV_NAME

using the following log files:
* failed proteins       $log, $tmp_log
* successfull proteins  $processed
"

echo "==> initializing log files"
for log_file in "$log" "$tmp_log" "$processed"
    test -w "$log_file" || touch "$log_file"
end

# resetiraj tmp log v vsakem primeru
echo "" > "$tmp_log"

# preveri, da ima res vsak protein iz seznama svoj output directory
echo "==> checking past analysis integrity ($processed)"
for protein in (cat "$processed")
    echo -n "$protein: "
    ls "$SWO" | grep -q "$protein"
    if test $status -ne 0
        echo "no output directory. removing from list"
        grep -v "$protein" "$processed" > "$processed" 2> /dev/null
    else
        echo "pass"
    end
end

echo "==> obtaining list of PDB files"
set files "$DB/*.pdb"
set n (count "$files")
set i 0
set failed 0

echo "==> activating conda environment"
conda activate "$CONDA_ENV_NAME"

echo "==> running analyses"

for protein in "$files"
    set i (math $i + 1)
    set prot_fname (path basename "$protein" | string replace ".pdb" "")
    set chain (string split "_" "$prot_fname" -f 2)
    echo -n "[$i/$n] $protein: "

    # če je protein že obdelan, pojdi na naslednjega
    grep -q "$prot_fname" "$processed" &&\
        echo "skip" && continue ||\
        echo -n "... "

    # shrani stderr v tmp datoteko
    "$SWORD_PATH" \
        -i "$protein" \
        -c "$chain" \
        -o "$SWO" \
    > /dev/null 2> "$tmp_log"

    if test $status -ne 0
        # v primeru, da gre nekaj narobe, shrani log
        echo "error. status: $status"
        set failed (math $failed + 1)
        echo -e "--- $(date "+%Y-%m-%d %H:%M") $protein ---\n" >> "$log"
        cat "$tmp_log" >> "$log"
        echo "" >> "$log"
    else
        # v primeru, da gre vse prav, dodaj protein na seznam
        echo "done"
        echo "$prot_fname" >> "$processed"
    end
end

echo "$failed/$n proteins failed"
echo "logs written to $log"
echo "results in $SWO"

# clean afer yourself
rm "$tmp_log"
conda deactivate
popd
