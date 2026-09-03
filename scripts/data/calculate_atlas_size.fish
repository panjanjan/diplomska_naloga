#!/bin/fish
# seznam proteinov v bazi. prenesi če še ne obstaja
if test ! -e "$ROOT/atlas_db/2024_11_18_ATLAS_pdb.txt"
    curl \
        "https://www.dsimb.inserm.fr/ATLAS/data/download/distributions/2024_11_18_ATLAS_pdb.txt"\
        -o "$ROOT/atlas_db/2024_11_18_ATLAS_pdb.txt"
end

set files (cat "$ROOT/atlas_db/2024_11_18_ATLAS_pdb.txt")
set total 0
set n (count $files)
set i 0
set failed 0

# izračuna velikost kompleta analysis podatkov za vsak protein
for protein in $files
    set i (math "$i + 1")
    echo -ne "\r[$i/$n] Checking: $protein. Total: $total_mb MB"

    # https://www.dsimb.inserm.fr/ATLAS/api/docs#/Downloads/download_atlas_analysis_ATLAS_analysis__pdb_chain__get
    set url "https://www.dsimb.inserm.fr/ATLAS/api/ATLAS/analysis/$protein"
    set size (
        curl -s -I -X 'GET' "$url" -H 'accept: application/octet-stream' |\
        grep content-length |\
        awk '{print $2}' |\
        string trim
    )

    if test -z "$size"
        set failed (math "$failed + 1")
        set size 0
        cat "$protein" >> "failed.log"
    end

    set total (math "$total + $size")
    set total_mb (math "$total / (1024^2)")
end

echo ""
echo "Total: $total_mb MB"
echo "Failed: $failed. Check failed.log"

echo "Total: $total_mb MB" > "atlas_db_size.txt"
