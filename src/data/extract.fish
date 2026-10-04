#!/bin/fish
# prejšnji extract je produkt mojega igranja ustvarjanja cli toola
# ta naredi kar potrebujem. to je vse.

# point je, da se ta script pokliče z xargs, zato sprejme ime datoteke
# preko -i/--input, npr: find ... | xargs -P 4 -I '{}' script -i '{}'
function _validate_input
    if not test -f "$_flag_value"
        echo "input nonexisting: $_flag_value" >&2
        return 1
    end
end

argparse -S 'i/input=!_validate_input' -- $argv; or exit 1

if not set -q _flag_input
    echo "./extract.fish -i/--input file" >&2
    exit 1
end

set zip_file "$_flag_input"
set pdb_dir "$ROOT/atlas_db/PDB"
set traj_dir "$ROOT/atlas_db/trajectories"
set tmp_dir "$ROOT/atlas_db/tmp"

test -d "$pdb_dir"; or mkdir -p "$pdb_dir"
test -d "$traj_dir"; or mkdir -p "$traj_dir"
test -d "$tmp_dir"; or mkdir -p "$tmp_dir"

function process_zip -a zipf
    set -l base (path basename --no-extension "$zipf")

    # - 1 PDB datoteka: npr. 1dd3_A.pdb
    # - 3 XTC datoteke: npr. 1dd3_A_R{1,2,3}.xtc
    set -l pdb_file "$pdb_dir/$base.pdb"
    set -l xtc_files "$traj_dir/$base"_R1.xtc "$traj_dir/$base"_R2.xtc "$traj_dir/$base"_R3.xtc

    # ne unzippaj če vse že obstaja
    set -l all_exist 1
    for f in $pdb_file $xtc_files
        if not test -f "$f"
            set all_exist 0
            break
        end
    end

    if test $all_exist -eq 1
        echo "$base: skip"
        return 0
    end

    # extractaj v začasni directory
    # q : quiet
    # d : directory
    set -l tmp_base "$tmp_dir/$base"
    mkdir -p "$tmp_base"

    set -l wanted_in_zip "$base.pdb" "$base"_R1.xtc "$base"_R2.xtc "$base"_R3.xtc
    if not unzip -q -o -d "$tmp_base" "$zipf" $wanted_in_zip
        echo "$base: unzip failed" >&2
        rm -rf "$tmp_base"
        return 1
    end

    # premakni vsako datoteko, če ne obstaja v target
    if not test -f "$pdb_file"; and test -f "$tmp_base/$base.pdb"
        mv "$tmp_base/$base.pdb" "$pdb_file"
    end
    for i in 1 2 3
        set -l target "$traj_dir/$base"_R$i.xtc
        set -l src "$tmp_base/$base"_R$i.xtc
        if not test -f "$target"; and test -f "$src"
            mv "$src" "$target"
        end
    end

    rm -rf "$tmp_base"
    echo "$base: done"
end

process_zip "$zip_file"
