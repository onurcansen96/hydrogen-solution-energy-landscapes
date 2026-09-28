# Hydrogen Solution-Energy Landscapes in Ferritic Fe

Open research data for **1,648 alternative single-H calculations** across
25 Fe host structures: a rotated single crystal and four grain-boundary
families, each with pure Fe, Cr, Cu, Ti, or a vacancy environment.
The calculations are from **DFT/VASP**, not LAMMPS simulations. LAMMPS-style
dump files are used only as a portable structure/site-property format.

**Paper status: unpublished.** The paper citation and DOI are pending and will
be added later. This repository is a data release, not a published-paper claim.

## Start Here

1. Download this repository with **Code > Download ZIP**, or clone it.
2. Read the existing JSON in `generated_metadata/`; no conversion is required.
3. For volumetric data, also download the CHGCAR archives from the
   [v1.0.0 release](https://github.com/onurcansen96/hydrogen-solution-energy-landscapes/releases/tag/v1.0.0).
   Git clone and GitHub's source-code ZIP do **not** include these archives.
4. Use [DATA_GUIDE.md](DATA_GUIDE.md) for definitions, units, crystallography,
   calculation protocols, source authority, and scientific limitations.

**Do not load all H rows as one physical configuration.** Each type-2 row in
an `*_allH` dump is an alternative H site. A physical starting structure is
all non-H rows plus exactly **one** selected H row. Frozen (`F`) means only H
relaxes; non-frozen (`NF`) means all ions relax, at fixed cell in both cases.

## Coverage

| Raw filename prefix | Paper label | Environment | H sites | CHGCAR files |
|---|---|---|---:|---:|
| `SC` | `SCR` | Rotated single crystal | 31 | 5 |
| `S3` | Sigma 3 | Symmetric-tilt GB | 113 | 5 |
| `S5` | Sigma 5 | Symmetric-tilt GB | 425 | **0** |
| `S27` | Sigma 27 | Symmetric-tilt GB | 471 | 5 |
| `S5m` | Sigma 5m | Mixed/near-CSL GB | 608 | 5 |
| Total | | 25 structure/chemistry combinations | **1,648** | **20** |

Every family has `PureFe`, `FeCr`, `FeCu`, `FeTi`, and `FeVac` site maps.
`SC_*` is the paper's **120-atom rotated SCR** reference, not its separate
128-atom bulk SC reference. Atom type 1 is Fe, type 2 is H, and type 3 is the
named solute where present. A vacancy is an absent Fe atom.

## What Is Included

| File or directory | Purpose |
|---|---|
| 25 `*_allH` files | Raw aggregate candidate-site maps, host cells/positions, and site properties |
| `Solution_Energy_Descriptors.txt` | Comma-separated source table, including corrected energies and Voronoi topology; the extension is `.txt`, not `.csv` |
| `generated_metadata/manifest_index.json` | Entry point: structure list, site-file paths, summaries, and output checksums |
| `generated_metadata/dataset.conceptual.json` | All 25 H-free hosts in conceptual-dictionary form |
| `generated_metadata/*/*.conceptual.json` | Per-host conceptual descriptions |
| `generated_metadata/*/*.hydrogen_sites.json` | One record per alternative H site, with separate F/NF properties and provenance |
| `generated_metadata/computational_settings.json` | Reported DFT and structure-preparation settings |
| `generated_metadata/field_dictionary.json` | Field meanings, units, source mappings, and normalization policy |
| `generated_metadata/schemas/` | JSON schemas for the site collections and manifest |
| `generated_metadata/source_manifest.json` | Scientific source artifacts, hashes, and publication-packaging notes |
| `generated_metadata/validation_report.json` | Supplied scientific validation report; its original timestamp is retained |
| `charge_density_manifest.json` | Inventory, download URLs, byte sizes, SHA-256 hashes, and CHGCAR header information |
| `examples/read_dataset.py` | Runnable reader and one-H starting-configuration example |
| `convert_h_solution_landscape.py` | Optional converter; NumPy required |
| `tools/verify_release.py` | Read-only integrity checks; Python standard library only |
| `provenance/` | Original author validation/source manifests and publication audit; see its README |

There are **no generated CSV outputs, standalone CONTCAR files, or notebook
files** in this release. The manuscript, representative INCAR, and external
descriptor-check CSV mentioned in historical provenance are **not included**.
Third-party source checkouts, caches, and local scratch files are not published.

## Read the Data

Python 3.10 or later is recommended. Reading the JSON does not require NumPy,
Jupyter, or the conceptual-dictionary package.

```bash
python examples/read_dataset.py
python tools/verify_release.py
```

In Python or a Jupyter notebook, with the repository root as working directory:

```python
import json
from pathlib import Path

root = Path("generated_metadata")
manifest = json.loads((root / "manifest_index.json").read_text(encoding="utf-8"))
entry = next(s for s in manifest["structures"]
             if s["paper_structure_id"] == "SCR_PureFe_allH")
collection = json.loads((root / entry["site_file"]).read_text(encoding="utf-8"))
site = collection["sites"][0]
frozen = site["calculations"]["frozen"]["properties"]
print(site["source_key"], frozen["solution_energy"]["value"])
print(frozen["solution_energy_zpe_corrected"]["value"])
```

`SE`/`SE_NF`, `ZPE`/`ZPE_NF`, and `SE_ZPEC_F`/`SE_ZPEC_NF` are in eV.
`BC`/`BC_NF` are reported Bader charge-transfer values in units of e, with the
source sign convention preserved. `CD` is a pre-insertion H-free-host
charge-density descriptor in e/angstrom^3. Consult the guide before analysis.

## Charge-Density Downloads

The release has four ZIP archives: `CHGCAR_SC.zip`, `CHGCAR_S3.zip`,
`CHGCAR_S27.zip`, and `CHGCAR_S5m.zip`. Each contains the five chemistry files
under their original `CHGCAR_<family>_<chemistry>` names, plus licensing and
attribution information. `SC` in these filenames means SCR in the paper.

There are **no `CHGCAR_S5_*` files** in the supplied data. The inventory maps
filenames to host identifiers and records header species/cells/grid dimensions;
it does not assert that every stored `CD` scalar has been independently
reproduced from these grids. The scalar converter does not read CHGCAR.
Charge integration and paper slice-plot scripts are not supplied.

After extracting the four archives into a directory, check the raw grids with:

```bash
python tools/verify_release.py --chgcar-dir /path/to/extracted/files
```

## Optional Regeneration

The distributed metadata is already usable. To regenerate in a **separate
copy** of the repository:

```bash
python -m pip install -r requirements.txt
python convert_h_solution_landscape.py --overwrite
```

`--overwrite` explicitly permits replacement of `generated_metadata/` only.
The scientific defaults are encoded in the supplied converter; it does not
infer them from a new manuscript. Without optional external sources, their
historical cross-checks cannot be repeated. See the guide for these distinctions.

## Quality and Reproducibility

Missing values are JSON `null`, never imputed. Six `SC_FeCr` sites have a
documented dump/table NF-ZPE conflict; both the reported corrected energies
and the conflict information are retained. The report also documents missing
Bader/ZPE values, unwrapped coordinates, and ignored properties on host rows.

The release supports analysis and new plots from the reported scalar data,
initial one-H structure reconstruction, and inspection of the supplied grids.
It does **not** claim that every manuscript figure or DFT result is exactly
reproducible: plotting/model-training scripts, per-site relaxed structures,
full calculation outputs, and some calculation details are absent.

## License, Attribution, and Citation

Data, JSON metadata, CHGCAR release assets, and documentation are licensed
under **CC BY 4.0** ([LICENSE-DATA](LICENSE-DATA)). Python code is licensed
under **MIT** ([LICENSE-CODE](LICENSE-CODE)). See [LICENSE.md](LICENSE.md).

Credit the creators listed in [CITATION.cff](CITATION.cff), identify the
repository and release version, and indicate modifications. The CFF describes
this **dataset**, not the unpublished paper. The paper citation and DOI remain
pending. Please report data questions through the repository's Issues page,
including the structure name and source atom ID.
