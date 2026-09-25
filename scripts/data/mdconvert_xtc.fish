#!/bin/fish
# bio3d ima možnost branja dcd, za xtc nima
# xtc trajektorije pretvori v dcd format
#
# python dependency: https://mdtraj.org/1.9.4/mdconvert.html
#
# sprejme ime xtc datoteke skozi stdin, da lahko uporabim xargs
# s skripto.
# --------------------------------------------------------------
mdconvert -h &> /dev/null || echo "mdconvert missing. run 'uv tool install
mdconvert'"

# --------------------------------------------------------------
function _validate_input
    if not test -f "$_flag_value"
        echo "input nonexisting: $_flag_value" >&2
        return 1
    end
end

argparse -S 'i/input=!_validate_input' -- $argv; or exit 1

if not set -q _flag_input
    echo "./mdconvert_xtc.fish -i/--input file" >&2
    exit 1
end

# --------------------------------------------------------------
set xtcfile "$_flag_input"
set pdb_dir "$ROOT/atlas_db/PDB_chained"
set traj_dir "$ROOT/atlas_db/trajectories"

# --------------------------------------------------------------
function process_xtc -a name
    set -l base_name (path basename "$name" --no-extension)
    set -l protein (string replace -r "_R.*" "" "$base_name")

    # preveri za obstoj PDB datoteke
    set -l pdbfile "$pdb_dir/$protein".pdb
    if not test -f "$pdbfile"
        echo "$base_name: PDB nonexistant" >&2
        return 1
    end

    # preveri za obstoj DCD trajektorije, da ne požene
    # mdconvert ponovno
    set -l dcdfile "$traj_dir/$base_name".dcd
    if test -f "$dcdfile"
        echo "$base_name: DCD exists, skip" >&2
        return 0
    end

    if not mdconvert -o "$dcdfile" -t "$pdbfile" "$name" &> /dev/null
        echo "$base_name: mdconvert failed" >&2
        return 1
    end

    echo "$base_name: done"
end

# --------------------------------------------------------------
process_xtc "$xtcfile"
