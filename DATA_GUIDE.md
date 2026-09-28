# Scientific Data Guide

This guide retains the scientific definitions supplied with the author's
metadata export. Packaging statements have been updated for the public
release. Start with [README.md](README.md) for downloads and runnable examples.
Scientific definitions and settings below are author-supplied descriptions,
not an independent re-analysis of the unpublished manuscript.

This release contains H-free atomistic host structures, alternative candidate
hydrogen positions, density-functional-theory-derived hydrogen properties, and
local geometric/electronic descriptors for a rotated single crystal and four
grain-boundary families in pure Fe, Fe-Cr, Fe-Cu, Fe-Ti, and vacancy-containing
Fe.

**Publication status:** unpublished. A formal citation and DOI have not yet
been assigned. The generated metadata therefore contains:

- `publication_status: "unpublished"`
- `citation: null`
- `doi: null`

The current manuscript is used only as an internal methodology source. It is
not presented as a published citation or preprint of record.

## Scope at a glance

| Raw family | Paper label | Environment | H-site records |
|---|---|---|---:|
| `SC` | `SCR` | bar(3)10-oriented rotated single crystal | 31 |
| `S3` | Sigma 3 | symmetric-tilt grain boundary | 113 |
| `S5` | Sigma 5 | symmetric-tilt grain boundary | 425 |
| `S27` | Sigma 27 | symmetric-tilt grain boundary | 471 |
| `S5m` | Sigma 5m | mixed/near-CSL grain boundary | 608 |
| **Total** |  | **25 structure-chemistry combinations** | **1,648** |

Each family is provided for:

- `PureFe`
- `FeCr`
- `FeCu`
- `FeTi`
- `FeVac`

All 1,648 candidate sites are published. The 45 sites excluded from the
surrogate-model training by EPOD screening remain part of the scientific data
release; model screening is not treated as deletion of raw DFT data.

### Important `SC`/`SCR` naming distinction

The raw `SC_*` files describe the paper's 120-atom rotated single crystal,
abbreviated `SCR`. They do **not** describe the separate 128-atom 4x4x4 bulk
single crystal called `SC` in the manuscript's computational table.

Raw names are retained as source identifiers, while paper-facing generated
names use `SCR_*`.

## Files used by this release

The scientific input layer consists of:

- 25 single-frame, LAMMPS-style `*_allH` aggregate site maps;
- `Solution_Energy_Descriptors.txt`, a comma-delimited descriptor table despite
  its `.txt` extension.

The scalar generator uses the two input categories above. In addition, 25
`CHGCAR_*` files are distributed in five downloadable release archives for
SC/SCR, S3, S5, S27, and S5m. Version 1.1.0 adds the five S5 files. Their
inventory and checksums are in `charge_density_manifest.json`. They are
additional raw resources, not inputs to the scalar converter: reported `CD`
values are imported from the descriptor table rather than recomputed.

After generation, `generated_metadata/` contains:

```text
generated_metadata/
  README.md
  dataset.conceptual.json
  computational_settings.json
  field_dictionary.json
  manifest_index.json
  source_manifest.json
  validation_report.json
  schemas/
    hydrogen_site_collection.schema.json
    manifest.schema.json
  <paper-structure>/
    <paper-structure>.conceptual.json
    <paper-structure>.hydrogen_sites.json
```

The generated files separate three data layers:

1. **Raw/imported:** dump coordinates and raw site columns, descriptor-table
   properties, and Voronoi topology.
2. **Derived:** geometric descriptors recomputed from the canonical
   `Vertices`, `Edges`, `Faces`, and `Voronoi_Volume` fields.
3. **Normalized:** explicit JSON `null` values restored from known missing-value
   sentinels, without numerical imputation.

## Critical aggregate-H semantics

A `*_allH` dump is an **aggregate candidate-site map**, not a physical
multi-hydrogen structure. Its type-2 rows are mutually alternative initial H
positions. Every type-2 row corresponds to a separate calculation containing
exactly one H atom.

The dump header's total atom count must therefore not be interpreted as the
atom count of a physical Fe-many-H sample.

### Atom-type interpretation

| Chemistry | Type 1 | Type 2 | Type 3 |
|---|---|---|---|
| PureFe | Fe | alternative H candidate | absent |
| FeVac | Fe | alternative H candidate | absent |
| FeCr | Fe | alternative H candidate | Cr |
| FeCu | Fe | alternative H candidate | Cu |
| FeTi | Fe | alternative H candidate | Ti |

Vacancies are represented by absent Fe atoms, not by a separate atom type.
Relative to the corresponding PureFe host:

- `SC_FeVac` contains one vacancy;
- each GB-family FeVac host contains two vacancies.

Similarly, SCR alloy hosts contain one substitutional solute, while the GB
alloy hosts contain two substitutional solutes.

Site-property columns are scientifically meaningful only on type-2 rows. Any
apparent site-property value on a non-H row is ignored by the generator.

## Reconstructing one physical calculation

To reconstruct the initial configuration for H site `i`:

1. Read the box and atomic coordinates from the relevant `*_allH` dump.
2. Retain every row whose `type != 2`; these atoms define the H-free host.
3. Select exactly one row satisfying `type == 2` and `id == i`.
4. Add that single H at the selected row's Cartesian position.
5. Preserve the dump cell and periodic boundary conditions.
6. Apply either the frozen or non-frozen ionic-relaxation protocol below.

In compact notation:

```text
physical sample i = {all rows with type != 2}
                  + {the one type-2 row with id i}
```

Dump coordinates and dump cells are authoritative.

Some host Cartesian coordinates are stored outside the nominal box bounds but
are periodic images of valid positions. The raw unwrapped coordinates are
preserved unchanged. A wrapped coordinate, if made by a downstream user, is a
derived representation and must not overwrite the source value.

## Frozen and non-frozen calculations

| State | Starting structure | Ionic constraints | Cell |
|---|---|---|---|
| Frozen (`F`) | reconstructed one-H structure | only H relaxes; all non-H atoms remain at their H-free dump positions | fixed |
| Non-frozen (`NF`) | same reconstructed one-H structure | all atoms relax | fixed |

"Frozen" does not mean that H is fixed. It means that the H-free host
environment is fixed while H is allowed to relax.

For the non-frozen calculation, all ions are free to relax from the same
one-H starting configuration. `ISIF=2` keeps cell shape and volume fixed. No
additional volume expansion is applied after H insertion.

The non-H dump positions are the corresponding relaxed, pre-insertion H-free
positions. They are not described as ideal lattice positions.

The release reconstructs the starting configuration and reported protocol. It
does not contain the relaxed H-containing output coordinates for every site.

## Solution-energy definitions

For branch `C`, where `C` is frozen or non-frozen:

```text
E_sol,C = E_host+H,C - E_host - 0.5 E_H2
```

The executable source uses:

- `E_H2 = -6.759828 eV`
- `0.5 E_H2 = -3.379914 eV`

The stored solution energy is a thermodynamic insertion/segregation quantity,
not a migration barrier.

### Zero-point energy

```text
ZPE = 0.5 * sum_i(h * nu_i)
```

The H-site ZPE uses the three local H vibrational modes. Frequencies were
evaluated with harmonic finite differences by displacing only H while keeping
all host atoms fixed. For a non-frozen site, the local-mode calculation is
performed about the non-frozen relaxed geometry. The H2 contribution contains
its stretching mode.

The corrected solution energy is:

```text
E_sol,C^ZPE = E_sol,C + ZPE_H,C - 0.5 ZPE_H2
```

The numerical dataset reference is:

- `0.5 ZPE_H2 = 0.1334825 eV`
- `ZPE_H2 = 0.266965 eV`

This numerical value is verified from executable source and the stored
corrected values; it is not printed numerically in the manuscript.

Corrected energies are imported as reported values and independently checked
with this expression. They are not silently replaced by recomputed values.

### Energy-field mapping

| Raw/source field | Release meaning | Unit |
|---|---|---|
| `SE`, `Atom_SE` | frozen solution energy | eV |
| `SE_NF`, `Atom_SE_NF` | non-frozen solution energy | eV |
| `ZPE`, `Atom_ZPE` | frozen local-H ZPE | eV |
| `ZPE_NF`, `Atom_ZPE_NF` | non-frozen local-H ZPE | eV |
| `Atom_SE_ZPEC` | `SE_ZPEC_F` | eV |
| `Atom_SE_NF_ZPEC` | `SE_ZPEC_NF` | eV |

## Charge density and Bader values

### Charge density (`CD`)

`CD` is evaluated for the candidate H region in the **H-free** host:

```text
CD = total electronic charge integrated over the candidate-H Voronoi region
     / volume of that candidate-H Voronoi region
```

Its unit is `e/Å³` (`e/angstrom^3` in ASCII-only contexts). Although stored on
an H candidate row, it is a pre-insertion host descriptor, not an output of the
H-containing relaxation. It is a reported precomputed scalar. Twenty-five
volumetric charge-density files are available as release assets, but the
integration workflow is not included and correspondence to every stored `CD`
value has not been independently verified.

### Bader values (`BC`, `BC_NF`)

`BC` and `BC_NF` are the reported H Bader charge-transfer values for the
frozen and non-frozen calculations, respectively, in units of `e`.

Values are preserved exactly as reported. The release performs no sign
reversal, subtraction of one electron, or conversion to a nominal ionic
charge.

## Geometric descriptors

The canonical descriptor table imports:

- candidate position;
- minimum H-host nearest-neighbor distance;
- closest-neighbor element;
- Voronoi volume `V`;
- original Voronoi face count;
- serialized vertices, edges, and triangulated hull faces.

The generator recomputes the paper's derived descriptors:

```text
V^(1/3)
A/V
sphericity = pi^(1/3) * (6V)^(2/3) / A
lambda1/lambda3
lambda2/lambda3
Rout/Rin
```

Here:

- `A` is the sum of the triangulated convex-hull face areas;
- the PCA eigenvalues satisfy `lambda1 >= lambda2 >= lambda3`;
- `Rout` and `Rin` are the maximum and minimum vertex radii relative
  to the vertex centroid.

`Num_Faces` is the original Voronoi face count. The serialized `Faces` field
contains triangulated hull facets and therefore has a different count.

The external `H_Voronoi_selected_descriptors_4SE.csv` was used only as a
numerical cross-check. It cannot override canonical dump coordinates, raw
topology, or descriptor-table properties. The broader
`H_Voronoi_descriptors.csv` is not used as a production source because its
charge-density values are stale for all 32 `S3_FeTi` sites. Neither of these
external CSV files is included here. This is historical provenance from the
supplied export, not a claim that either CSV is available in the repository.

## Crystallography and reference GB energies

| Raw family | Plane/orientation | Axis | Explicit angle | Angle status | Reported mixed components | Reference GB energy |
|---|---|---|---:|---|---|---:|
| `SC` / `SCR` | `[-3,1,0]` orientation | `[0,0,1]` | 18.434948823 deg | derived orientation rotation; not a GB misorientation | not applicable | not applicable |
| `S3` | `[1,1,-2]` | `[1,-1,0]` | 70.528779366 deg | derived ideal CSL axis-angle | none | 0.41 J/m2 |
| `S5` | `[-3,1,0]` | `[0,0,1]` | 36.869897646 deg | derived ideal CSL axis-angle | none | 1.54 J/m2 |
| `S27` | `[5,-5,2]` | `[1,1,0]` | 31.586338097 deg | derived ideal CSL axis-angle | none | 1.68 J/m2 |
| `S5m` | `[1,1,-2]` | `[1,0,0]` | 36.869897646 deg | derived ideal Sigma-5 CSL axis-angle | tilt 33.6 deg; twist 0.27 deg | 2.09 J/m2 |

The SCR cell orientation is represented as:

- cell x parallel to `[-3,1,0]`;
- cell y parallel to `[0,0,1]`;
- cell z parallel to `[1,3,0]`.

This is a single-crystal orientation descriptor, not a GB misorientation.

The explicit derived angles use the following crystallographic relations:

```text
SCR  = atan(1/3)              = 18.434948823 deg
S3   = 2 atan(sqrt(2)/2)      = 70.528779366 deg
S5   = 2 atan(1/3)            = 36.869897646 deg
S27  = 2 atan(sqrt(2)/5)      = 31.586338097 deg
S5m  = ideal Sigma-5 angle    = 36.869897646 deg
```

For SCR, the stored signed angle is `-18.434948823 deg` from `[-1,0,0]`
toward `[-3,1,0]` using the right-hand rule about positive `[0,0,1]`. The
manuscript reports the S5m tilt and twist components numerically. It does not
print the ideal CSL angles for S3, S5, S27, or S5m, nor the SCR orientation
angle; those values are therefore marked `derived`. The S5m tilt and twist
components must not be added together or substituted for the CSL axis-angle.
The CSL numbers in the table use the lower-angle representative in the
`[0,90] deg` interval for each listed ideal plane/axis construction;
symmetry-equivalent complementary conventions may appear in other sources.
The S5m value `0.27 deg` and its degree unit are retained exactly as printed in
the manuscript and are not reinterpreted or converted.

The GB energies are family-level reference values reported for the
corresponding boundary models. They are not independently recalculated GB
energies for every alloy, vacancy, H site, or relaxed one-H configuration.

## DFT and structure-preparation settings

### Global settings

| Setting | Value |
|---|---|
| Code | VASP; version not recorded |
| Core-valence treatment | PAW |
| Exchange-correlation | GGA-PBE |
| Spin treatment | collinear spin-polarized |
| `ENCUT` | 400 eV |
| `EDIFF` | `1e-6 eV` |
| `EDIFFG` | `-0.01 eV/angstrom` |
| `ISTART` | 0 |
| `ALGO` | Fast |
| `ISMEAR` | 1 |
| `SIGMA` | 0.1 eV |
| `IBRION` | 2 |
| `POTIM` | 0.2 |
| `NSW` | 1000 |
| `ISIF` | 2 |
| `LWAVE` | false |
| `LCHARG` | true |
| `LVTOT` | true |
| `LELF` | true |
| `LORBIT` | 10 |
| `NBANDS` | VASP default; the supplied line is commented out |
| Equilibrium bcc Fe lattice constant | 2.833 angstrom |

The representative INCAR contains `EDIFF=1e-5`, while the manuscript and the
author's explicit correction specify `1e-6 eV` for this dataset. The generated
metadata uses `1e-6 eV` as authoritative and records the discrepancy in
provenance.

Frozen/non-frozen atom constraints are not encoded by the shared INCAR tags.
The frozen branch requires atom-specific selective-dynamics constraints or an
equivalent workflow rule.

### Initial magnetic moments

| Element | Initial moment |
|---|---:|
| Fe | +3 mu_B |
| Cr | -3 mu_B |
| H | 0 mu_B |
| Ti | 0 mu_B |
| Cu | 0 mu_B |

These are initial VASP moments, not claims about converged local magnetic
moments. A per-calculation `MAGMOM` list must follow the reconstructed one-H
composition and atom ordering.

### Family-specific settings

| Family | Manuscript initial cell (angstrom) | k mesh | Normal-volume optimization | RBT |
|---|---|---|---|---|
| SCR | 26.87 x 5.66 x 8.95 | 1 x 6 x 4 | yes | no/not applicable |
| S3 | 27.75 x 8.01 x 4.90 | 1 x 4 x 6 | yes | yes |
| S5 | 26.87 x 5.66 x 8.95 | 1 x 6 x 4 | yes | yes |
| S27 | 41.63 x 8.01 x 7.36 | 1 x 4 x 4 | yes | no |
| S5m | 20.81 x 6.33 x 15.51 | 1 x 3 x 1 | yes | no |

These are manuscript-reported **initial** dimensions. Reconstruction uses the
actual per-file dump cell, which is authoritative after structure preparation
and volume optimization.

The manuscript's `Single Site` segregation-optimization entry for SCR/S5m
does not mean that only one H candidate is present in this release. It refers
to structure/solute or vacancy optimization; all identified H sites are
published.

## Missing-value policy

Missing scientific values are serialized as JSON `null`. They are never
emitted as IEEE `NaN`, the strings `"nan"`/`"NA"`, or an artificial numerical
sentinel. No value is imputed.

Known normalization:

- 6 missing frozen Bader values become `null`;
- 2 missing non-frozen ZPE values and the corresponding corrected energies
  become `null`;
- 52 S27 non-frozen Bader values that were upstream `NA` and later zero-filled
  are restored to `null`.

Every affected `(Structure, Atom_id)` key is listed in
`generated_metadata/validation_report.json`.

## Known source anomalies

The release preserves and reports rather than hiding these anomalies:

1. For all six `SC_FeCr` H sites (atom IDs 120-125), the dump reports a
   non-zero `ZPE_NF`, while `Solution_Energy_Descriptors.txt` stores
   `Atom_ZPE_NF=0.0`. Its reported `Atom_SE_NF_ZPEC` was calculated with the
   table's zero. The site JSON preserves the dump ZPE, preserves the reported
   corrected energy, stores both source values, and adds a conflict flag.
2. Thirteen dumps contain a non-zero `BC` on non-H atom ID 2. Site-property
   columns on all non-H rows are ignored.
3. Some host coordinates are unwrapped outside nominal box bounds. They are
   periodic images and remain unchanged.
4. `H_Voronoi_descriptors.csv` contains stale charge density for the 32
   `S3_FeTi` rows and is not used as a production source.

## Source authority and provenance

When sources overlap, the release applies this order:

1. Dump: cell, atom ID, LAMMPS type, Cartesian coordinates, and primary raw
   `ZPE_NF` where it conflicts with the table.
2. `Solution_Energy_Descriptors.txt`: corrected energies, raw Voronoi topology,
   and tabulated site properties.
3. Recomputed geometry: derived descriptors calculated from canonical topology.
4. Manuscript: crystallographic plane/axis, reported S5m tilt/twist components,
   family-level GB energy, and computational description.
5. Derived crystallography: explicit ideal CSL axis-angles and the SCR
   orientation rotation, stored with formulas and provenance status.
6. Author clarification: `EDIFF`, elemental initial moments, Bader convention,
   charge-density definition/unit, coordinate authority, and missing policy.

Every source artifact is assigned a stable artifact ID, byte size, and SHA-256
checksum in `source_manifest.json`. Every site links to its source dump,
descriptor-table row, atom ID, and descriptor index.

Paths for included artifacts are relative to the canonical folder. For the
unpublished manuscript, representative INCAR, and external descriptor
cross-check, the manifest deliberately records `path: null` while retaining
the filename, role, size, and checksum. This avoids embedding machine-local
directory layouts in the release.

Derived data never silently overwrite raw data.

## Conceptual-dictionary representation

`dataset.conceptual.json` and each per-structure `*.conceptual.json` follow the
current `OCDO/conceptual_dictionary` top-level shape:

```text
dataset
computational_sample
workflow
operation
math_operation
```

The native conceptual layer describes:

- the 25 H-free computational hosts;
- material composition and ideal bcc reference phase;
- dump-derived simulation cells and host positions;
- native vacancy/substitutional and GB blocks;
- family-level GB-energy properties;
- DFT workflows using `method: DensityFunctionalTheory`, `algorithm: null`,
  and `xc_functional: PBE`.

The aggregate raw dump is not used as `atom_attribute.file_path`, because a
generic LAMMPS reader would load every candidate H as simultaneous occupancy.
Host positions/species are therefore embedded inline after filtering out all
type-2 rows.

Information that is not fully expressible by the generic conceptual template
is placed under the documented `x_hydrogen_landscape` extension namespace:

```text
urn:ocdo:hydrogen-solution-landscape:1
```

This extension carries, among other things:

- the exact one-H selection rule;
- inserted H species, source ID, and Cartesian position;
- atom-specific frozen/non-frozen constraints;
- explicit angle metadata for SCR, S3, S5, S27, and S5m, including reported
  versus derived status;
- separate reported S5m tilt/twist components and the SCR cell orientation;
- composite charge-density unit;
- missing-value status and detailed provenance.

The generic native `AddAtom` operation stores input/output sample IDs but not
the inserted species and position, which is why these details remain in the
extension sidecar.

## Validation

The generator performs and records the following checks:

- exactly 25 supported single-frame dumps;
- exactly 1,648 unique type-2 site records;
- complete one-to-one join by `(Structure, Atom_id)`;
- descriptor indices unique and continuous;
- dump/table coordinate agreement within `1e-6 angstrom`;
- raw scalar agreement within `1e-6` in native units, except the documented
  `SC_FeCr` NF-ZPE conflict;
- corrected-energy equation checks within `1e-6 eV` source precision;
- geometry recomputation and selected-descriptor cross-check;
- valid vertex/edge/face indices;
- host position/species/ID/type array consistency;
- expected solute/vacancy counts and type maps;
- complete and internally consistent angle metadata for all 25 structures and
  all five families;
- JSON `null` handling with no NaN/Infinity;
- source and output SHA-256 checksums;
- strict controlled-vocabulary validation of all 26 conceptual documents with
  `ConceptualDict.validate` when the reference package is available.

The generated validation report distinguishes errors from expected warnings.
A successful export may therefore have status `passed_with_warnings` when it
contains documented raw-source conflicts or preserved coordinate conditions.

## Running the generator

Python 3.10 or later and NumPy are required. The existing JSON needs no
regeneration. To regenerate in a separate copy of the repository:

```bash
python -m pip install -r requirements.txt
python convert_h_solution_landscape.py --overwrite
```

The default input root is the directory containing the script. The generator
reads only the 25 canonical filenames in that directory; it does not scan
recursively and cannot double-count mirrored project copies.

An existing `generated_metadata/` requires explicit `--overwrite`; the safety
check restricts replacement to this named folder directly inside the input
root. Timestamps and artifact hashes change on regeneration.

Optional `--paper`, `--incar`, and `--selected-descriptors` paths refer to
external files not distributed here. The first two record artifact provenance,
not automatic extraction of calculation settings. `--selected-descriptors`
adds numerical cross-checks. `--conceptual-dictionary-source` can point to a
separately obtained OCDO/conceptual_dictionary checkout for strict vocabulary
validation (the author report records commit
`4a41fda37a46e25dad2f70e273a67a517147409e`). Without it, that check is skipped.
The original successful external checks are historical results, not repeatable
from the public files alone. See `provenance/README.md` for release-time checks.

## Reproducibility and limitations

This release reconstructs the intended one-H starting structures and associates
each candidate with its reported scalar results. It does not contain all:

- relaxed H-containing output geometries;
- raw total energies for every frozen/non-frozen branch;
- force histories;
- vibrational Hessians;
- the integration scripts needed to independently reproduce every stored
  `CD` value from the supplied grids.

Exact numerical reproduction may additionally depend on details not recorded
in the available sources, including:

- VASP version;
- exact PAW/POTCAR labels and releases;
- k-point centering;
- isolated-H2 cell and geometry;
- finite-displacement amplitude;
- Bader implementation and grid settings.

The release therefore says that it reconstructs the intended starting
configuration and reported protocol. It does not promise bitwise-identical
energies from an independent rerun.

The ZPE approximation contains local H modes, not the full host phonon
spectrum, and neglects anharmonicity and explicit H-host vibrational coupling.

The dataset covers selected ferritic-Fe GB and dilute solute/vacancy
environments. It is not a complete sampling of all ferritic-steel grain
boundaries.

## Citation and DOI

The associated manuscript is not yet published:

- Citation: `null`
- DOI: `null`

These fields can be updated after publication. No provisional journal citation,
publication year or DOI is invented by this release. The author-approved
licenses are CC BY 4.0 for data/documentation and MIT for code; see LICENSE.md.
