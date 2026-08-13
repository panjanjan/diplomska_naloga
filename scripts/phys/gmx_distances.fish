#!/bin/fish
# z gromacs izračuna razdalje med masnima centroma domen skozi trejektorijo

pushd "$ROOT"

set traj "$ROOT/atlas_db/TRAJ"
set target "$ROOT/atlas_db/COM"
set target_tmp "$target/tmp"

test -d $target || mkdir -p $target
test -d $target_tmp || mkdir -p $target_tmp
rm -rf $target_tmp/*

function run -a protein
    # 1. dobi trajektorijo (xtc) in topologijo (tpr)
    unzip -q "atlas_db/analysis/$protein" -d $target_tmp
    find $target_tmp -name "*.xtc" -exec cp {} $target \;
    find $target_tmp -name "*.tpr" -exec cp {} $target \;

    # 2. dobi indekse od domen
    $ROOT/scripts/domain_inds.r $protein | read -L s1 s2 e1 e2
    set sel1 "res_com of resnr $s1 to $e1"
    set sel2 "res_com of resnr $s2 to $e2"

    # 3. dobi razdalje
    for i in 1 2 3
        set traj_f {$traj}/{$protein}_R{$i}.xtc
        set top_f {$traj}/{$protein}_R{$i}.tpr
        set out_f {$traj}/{$protein}_dist_R{$i}

        echo traj_f top_f out_f
        continue

        gmx pairdist \
            -f "$traj_f" \
            -s "$top_f" \
            -ref "$sel1" \
            -sel "$sel2" \
            -xvg none \
            -o "$out_f" 2> /dev/null
    end

    rm -rf $target_tmp/*
end

set protein_list (cut -d ',' -f1 two_domains.csv)
set i 1
set n (count $protein_list)

for protein in $protein_list
    echo -e "\r\x1b[K[$i/$n] $protein ..."
    ls $target | grep -q "^$protein.*xvg\$" || run $protein
    set i (math $i + 1)
end

rm -rf "$target_tmp"

popd
