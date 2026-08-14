> ZASTARELO

[`data`](./data/): pridobivanje podatkov

- [`calculate_atlas_size.fish`](./data/calculate_atlas_size.fish): izračuna velikost baze. V primeru da na disku ni dovolj prostora, bo problem naslednji korak.
- [`download_atlas.fish`](./data/download_atlas.fish): prenese `analysis` del ATLAS MD baze
- [`extract.fish`](./data/extract.fish): prenešene so zip datoteke. Ker ne potrebujemo vseh vključenih datotek, ta skripta izvleče samo določene glede na parameter. Trenutno so pomembne RMSF (`*_RMSF.tsv`), PDB (`*.pdb`), trajektorije (`*.xtc`) in topologije (`*.tpr`)
- [`mdconvert_xtc.fish`](./data/mdconvert_xtc.fish): trajektorije pretvori v format, ki ga je bio3d zmožen prebrati.

[`domains`](./domains/): določanje in izbiranje proteinov s kvalitetnimi dekompozicijami (google SWORD2)

- [`introduce_chain.r`](./data/introduce_chain.r): v PDB datoteke, ki jih sprejme spodnji script, doda podatek o verigi. Ta korak je pomemben, saj SWORD2 drugače ne bo deloval.
- [`sword2_example.fish`](./domains/sword2_example.fish): da lahko preveriš če ti deluje SWORD2.
- [`sword2_batch_processor.fish`](./domains/sword2_batch_processor.fish): požene SWORD2 nad proteini, da določi domene in ostale metrike.
- [`build_decompositions_csv.py`](./domains/build_decompositions_csv.py): iz JSON datotek, ki jih vrne SWORD2, naredi pregleden csv.
- [`filtering.rmd`](./domains/filtering.rmd): prečisti rezultate, da ohrani proteine s kvalitetnimi in primernimi dekompozicijami za nadaljne analize.
- [`two_domains.r`](./domains/two_domains.r): naredi pregleden csv, kjer so shranjeni podatki o mejah domen. Pomemben input za nadaljne skripte.
- [`domain_inds.r`](./domains/domain_inds.r): tebe lahko izbrišem*?

[`rmsf`](./rmsf/): analize RMSF vrednosti

- [`rmsf_analysis.r`](./rmsf/rmsf_analysis.r): zastarelo. Gitignore?
- [`rmsf_analysis.rmd`](./rmsf/rmsf_analysis.rmd): utility. Gitignore?
- [`rmsf_test.r`](./rmsf/rmsf_test.r): REFACTOR, ODSTRANI V2. izvede statistični test nad RMSF-ji domen
- [`rmsf_ratios.r`](./rmsf/rmsf_ratios.r): preveri razmerja med notranjo in zunanjo fleksibilnostjo domen.

[`phys`](./phys/): analize, ki vključujejo več..fizike

- [`com_distances_bio3d.r`](./phys/com_distances_bio3d.r): določi masne centre domen in izračuna njune razdalje skozi trajektorijo.
- [`com_distances_gmx.r`](./phys/com_distances_gmx.r): DEPRECATED
- [`gmx_distances.fish`](./phys/gmx_distances.fish): DEPRECATED
- [`pai_analysis.r`](./phys/pai_analysis.r): čaka na nadaljno implementacijo
- [`principal_axes_of_inertia.r`](./phys/principal_axes_of_inertia.r): določi vztrajnostne osi domen in izračuna kote med njimi skozi trajektorijo.
- [`sde.r`](./phys/sde.r): spektralna analiza razdalj, kotov.
