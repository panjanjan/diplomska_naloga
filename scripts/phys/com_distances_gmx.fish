#!/bin/fish
# Izračuna razdalje med masnima centroma obeh domen za vsak protein iz
# two_domains.csv skozi trajektorije v atlas_db/trajectories z GROMACS.
#
# Trajektorije morajo biti XTC (glej extract.fish) in PDB datoteke (PDB_chained)
# morajo imeti enako število atomov kot trajektorije. Vsak frame predhodno
# superimpoziramo na domeno 1.
#
# Pred prvim zagonom poženi build_ndx.py, ki zgradi indeksne datoteke
# (outputs/COM_gmx/index/{protein}.ndx).
#
# Izhodi:
#   outputs/COM_gmx/{protein}_dist.csv  frame,R1,R2,R3 v nm
#   outputs/COM_gmx/xvg/{protein}_R{i}.xvg  če je -k
#   outputs/COM_gmx/failed.log             neuspeli proteini
#
# Uporaba:
#   com_distances_gmx.fish [-t/--test n] [-j/--jobs n] [-k/--keep-xvg]
# ---------------------------------------------------------------------
argparse -S \
    't/test='   \
    'j/jobs='   \
    'k/keep-xvg' \
    -- $argv
or exit 1

pushd "$ROOT"

set gmx /usr/local/gromacs/bin/gmx
set traj_dir "$ROOT/atlas_db/trajectories"
set pdb_dir  "$ROOT/atlas_db/PDB_chained"
set domains_csv "$ROOT/outputs/two_domains.csv"
set out_dir   "$ROOT/outputs/COM_gmx"
set index_dir "$out_dir/index"
set tmp_dir   "$out_dir/tmp"
set log       "$out_dir/failed.log"

test -d "$out_dir"   || mkdir -p "$out_dir"
test -d "$index_dir" || mkdir -p "$index_dir"
test -d "$tmp_dir"   || mkdir -p "$tmp_dir"
touch "$log"

set jobs 4
if set -ql _flag_jobs
    set jobs (math $_flag_jobs)
end

# argparse lokalne spremenljivke niso vidne v ozadnih opravilih,
# zato jih prenesemo v globalne pred zagonom
set -g _keep_xvg 0
if set -ql _flag_keep_xvg
    set -g _keep_xvg 1
end

# ------------------------------------------------ razširi izbiro proteinov
set proteins (awk -F, 'NR>1 {print $1}' "$domains_csv" | sort -u)
if set -ql _flag_test
    set n (math $_flag_test)
    set proteins $proteins[1..$n]
end
set n (count $proteins)
if test $n -eq 0
    echo "no proteins in $domains_csv" >&2
    exit 1
end

echo "processing $n proteins (jobs: $jobs)"
echo "outputs -> $out_dir"

# ------------------------------------------------ en protein, 3 replicati
function run_protein -a protein
    set pdb "$pdb_dir/$protein.pdb"
    set idx "$index_dir/$protein.ndx"
    set csv "$out_dir/$protein"_dist.csv
    set pdir "$tmp_dir/$protein"

    if test -f "$csv"
        echo "$protein: skip (exists)"
        return 0
    end

    if not test -f "$idx"
        echo "--- $protein: missing index ($idx) ---" >> "$log"
        echo "  poženi build_ndx.py pred com_distances_gmx.fish" >> "$log"
        return 1
    end

    mkdir -p "$pdir"

    set dists
    set ok 1
    for i in 1 2 3
        set xtc  "$traj_dir/$protein"_R$i.xtc
        if not test -f "$xtc"
            echo "$protein: missing R$i ($xtc)" >> "$log"
            set ok 0
            break
        end

set whole "$pdir/whole_R$i.xtc"
        set fit   "$pdir/fit_R$i.xtc"
        if test $_keep_xvg -eq 1
            mkdir -p "$out_dir/xvg"
            set xvg "$out_dir/xvg/$protein"_R$i.xvg
        else
            set xvg "$pdir/dist_R$i.xvg"
        end

        printf 'System\n' | "$gmx" trjconv -s "$pdb" -f "$xtc" -o "$whole" -pbc whole \
            > /dev/null 2> "$pdir/err_R$i.log"
        or begin
            echo "--- $protein R$i trjconv whole ---" >> "$log"
            cat "$pdir/err_R$i.log" >> "$log"
            set ok 0
            break
        end

        printf 'D1\nSystem\n' | "$gmx" trjconv -s "$pdb" -f "$whole" -o "$fit" \
            -fit rot+trans -n "$idx" > /dev/null 2>> "$pdir/err_R$i.log"
        or begin
            echo "--- $protein R$i trjconv fit ---" >> "$log"
            cat "$pdir/err_R$i.log" >> "$log"
            set ok 0
            break
        end

        "$gmx" distance -s "$pdb" -f "$fit" -n "$idx" \
            -select 'com of group "D1" plus com of group "D2"' \
            -oall "$xvg" -xvg none > /dev/null 2>> "$pdir/err_R$i.log"
        or begin
            echo "--- $protein R$i distance ---" >> "$log"
            cat "$pdir/err_R$i.log" >> "$log"
            set ok 0
            break
        end

        # stolpec razdalj za združevanje
        awk '!/^[#@]/ {print $2}' "$xvg" > "$pdir/d_$i.txt"
        set -a dists "$pdir/d_$i.txt"
    end

    if test $ok -eq 0
        rm -rf "$pdir"
        return 1
    end

    # združi R1-R3 v CSV (ime frame/R1/R2/R3, razdalje v nm)
    if test (count $dists) -eq 3
        echo "frame,R1,R2,R3" > "$csv"
        paste -d ' ' $dists > "$pdir/dists.txt"
        awk '{printf "%d,%s,%s,%s\n", NR, $1, $2, $3}' "$pdir/dists.txt" >> "$csv"
        echo "$protein: ok -> $csv"
    else
        echo "$protein: incomplete replicates" >> "$log"
    end

    rm -rf "$pdir"
end

# ------------------------------------------------ vzporedno izvajanje
set running 0
for protein in $proteins
    run_protein $protein &
    set running (math $running + 1)
    if test $running -ge $jobs
        wait
        set running 0
    end
end
wait

# ------------------------------------------------ poročilo
set missing 0
for protein in $proteins
    if not test -f "$out_dir/$protein"_dist.csv
        echo "missing: $protein"
        set missing (math $missing + 1)
    end
end
echo
echo "done: $n proteins, $missing missing (see $log)"

popd
