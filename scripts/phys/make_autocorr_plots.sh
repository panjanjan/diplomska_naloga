#!/bin/sh
# Regenerate autocorr_plots.html from outputs/ACF/*.png.
# Each protein has 3 PNGs (<protein>_bb|ca|noh_acf.png).
# Layout: one row per protein, 3 side-by-side columns.
# Run from repo root: sh scripts/phys/make_autocorr_plots.sh
set -eu

ACF_DIR="$ROOT/outputs/ACF"
OUT="$ROOT/outputs/autocorr_plots.html"

echo '<!DOCTYPE html>' > "$OUT"
echo '<html lang="en">' >> "$OUT"
echo '<head>' >> "$OUT"
echo '    <meta charset="UTF-8">' >> "$OUT"
echo '    <title>Autocorrelation Plots</title>' >> "$OUT"
echo '    <style>' >> "$OUT"
echo '        body { font-family: sans-serif; margin: 20px; }' >> "$OUT"
echo '        .row { display: flex; align-items: center;' >> "$OUT"
echo '               margin-bottom: 30px;' >> "$OUT"
echo '               border-bottom: 1px solid #ccc;' >> "$OUT"
echo '               padding-bottom: 10px; }' >> "$OUT"
echo '        .protein-label { width: 100px; font-weight: bold;' >> "$OUT"
echo '                         flex-shrink: 0; }' >> "$OUT"
echo '        .img-container { display: flex; flex-grow: 1;' >> "$OUT"
echo '                         justify-content: space-around; }' >> "$OUT"
echo '        .img-container div { text-align: center; width: 32%; }' >> "$OUT"
echo '        img { width: 100%; height: auto; }' >> "$OUT"
echo '    </style>' >> "$OUT"
echo '</head>' >> "$OUT"
echo '<body>' >> "$OUT"
echo '    <h1>Autocorrelation Plots</h1>' >> "$OUT"

# Loop over one PNG per protein (*_bb_acf.png) to get the
# protein id, then emit all 3 columns (bb | ca | noh) in a
# single row.
for bb in "$ACF_DIR"/*_bb_acf.png; do
  [ -e "$bb" ] || continue
  protein=$(basename -- "$bb" _bb_acf.png)
  ca="$ACF_DIR/${protein}_ca_acf.png"
  noh="$ACF_DIR/${protein}_noh_acf.png"
  if [ ! -f "$ca" ] || [ ! -f "$noh" ]; then
    echo "skip $protein: missing ca/noh png" >&2
    continue
  fi
  echo '    <div class="row">' >> "$OUT"
  echo "        <div class=\"protein-label\">$protein</div>" >> "$OUT"
  echo '        <div class="img-container">' >> "$OUT"
  echo "            <div><p>Backbone</p><img src=\"$ACF_DIR/${protein}_bb_acf.png\" alt=\"$protein BB\"></div>" >> "$OUT"
  echo "            <div><p>C-alpha</p><img src=\"$ACF_DIR/${protein}_ca_acf.png\" alt=\"$protein CA\"></div>" >> "$OUT"
  echo "            <div><p>Non-H</p><img src=\"$ACF_DIR/${protein}_noh_acf.png\" alt=\"$protein NOH\"></div>" >> "$OUT"
  echo '        </div>' >> "$OUT"
  echo '    </div>' >> "$OUT"
done

echo '</body>' >> "$OUT"
echo '</html>' >> "$OUT"
echo "wrote $OUT"
