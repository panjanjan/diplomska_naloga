> Posodobi glede na kapitana in core.

Skripte so organizirane po fazi analize. Vsak direktorij vsebuje `main.sh`, ki lahko služi kot referenca na vrstni red skript.

[`00_data`](./00_data/): pridobivanje podatkov

- [`calculate_atlas_size.fish`](./00_data/calculate_atlas_size.fish): izračuna velikost baze. V primeru da na disku ni dovolj prostora, bo problem naslednji korak.
- [`download_atlas.fish`](./00_data/download_atlas.fish): prenese `analysis` del ATLAS MD baze
- [`extract.fish`](./00_data/extract.fish): prenešene so zip datoteke. Ker ne potrebujemo vseh vključenih datotek, ta skripta izvleče samo določene glede na parameter. Trenutno so pomembne RMSF (`*_RMSF.tsv`), PDB (`*.pdb`), trajektorije (`*.xtc`) in topologije (`*.tpr`)
- [`mdconvert_xtc.fish`](./00_data/mdconvert_xtc.fish): trajektorije pretvori v format, ki ga je bio3d zmožen prebrati.

[`01_domains`](./01_domains/): določanje in izbiranje proteinov s kvalitetnimi dekompozicijami (google SWORD2)

- [`introduce_chain.r`](./00_data/introduce_chain.r): v PDB datoteke, ki jih sprejme spodnji script, doda podatek o verigi. Ta korak je pomemben, saj SWORD2 drugače ne bo deloval.
- [`sword2_example.fish`](./01_domains/sword2_example.fish): da lahko preveriš če ti deluje SWORD2.
- [`sword2_batch_processor.fish`](./01_domains/sword2_batch_processor.fish): požene SWORD2 nad proteini, da določi domene in ostale metrike.
- [`build_decompositions_csv.py`](./01_domains/build_decompositions_csv.py): iz JSON datotek, ki jih vrne SWORD2, naredi pregleden csv.
- [`filtering.rmd`](./01_domains/filtering.rmd): prečisti rezultate, da ohrani proteine s kvalitetnimi in primernimi dekompozicijami za nadaljne analize.
- [`two_domains.r`](./01_domains/two_domains.r): naredi pregleden csv, kjer so shranjeni podatki o mejah domen. Pomemben input za nadaljne skripte.
- [`domain_inds.r`](./01_domains/domain_inds.r): tebe lahko izbrišem*?

[`02_rmsf`](./02_rmsf/): analize RMSF vrednosti

- [`rmsf_analysis.r`](./02_rmsf/rmsf_analysis.r): zastarelo. Gitignore?
- [`rmsf_analysis.rmd`](./02_rmsf/rmsf_analysis.rmd): utility. Gitignore?
- [`rmsf_test.r`](./02_rmsf/rmsf_test.r): REFACTOR, ODSTRANI V2. izvede statistični test nad RMSF-ji domen
- [`rmsf_ratios.r`](./02_rmsf/rmsf_ratios.r): preveri razmerja med notranjo in zunanjo fleksibilnostjo domen.

[`03_phys`](./03_phys/): analize, ki vključujejo več..fizike

- [`com_distances_bio3d.r`](./03_phys/com_distances_bio3d.r): določi masne centre domen in izračuna njune razdalje skozi trajektorijo.
- [`com_distances_gmx.r`](./03_phys/com_distances_gmx.r): DEPRECATED
- [`gmx_distances.fish`](./03_phys/gmx_distances.fish): DEPRECATED
- [`pai_analysis.r`](./03_phys/pai_analysis.r): čaka na nadaljno implementacijo
- [`principal_axes_of_inertia.r`](./03_phys/principal_axes_of_inertia.r): določi vztrajnostne osi domen in izračuna kote med njimi skozi trajektorijo.
- [`sde.r`](./03_phys/sde.r): spektralna analiza razdalj, kotov.
