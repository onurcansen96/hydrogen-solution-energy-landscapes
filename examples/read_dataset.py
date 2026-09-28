"""Read energies and reconstruct one initial one-H configuration without writing files."""

import argparse
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def reconstruct_one_h(dump_path, atom_id):
    """Return source orthogonal bounds and one H plus host, without wrapping.

    This is an initial configuration, not the relaxed output of DFT.
    """
    lines = Path(dump_path).read_text(encoding="utf-8").splitlines()
    if lines[4] != "ITEM: BOX BOUNDS pp pp pp":
        raise ValueError("This example accepts only the released orthogonal cells")
    bounds = [[float(v) for v in row.split()] for row in lines[5:8]]
    columns = lines[8].split()[2:]
    rows = [dict(zip(columns, row.split())) for row in lines[9:]]
    if len(rows) != int(lines[3]):
        raise ValueError("Dump atom-count mismatch")
    selected = [r for r in rows if int(r["type"]) == 2 and int(r["id"]) == atom_id]
    if len(selected) != 1:
        raise ValueError(f"Expected exactly one type-2 row with id {atom_id}")
    chemistry = Path(dump_path).name.split("_")[1]
    species = {1: "Fe", 2: "H"}
    if chemistry in ("FeCr", "FeCu", "FeTi"):
        species[3] = chemistry[2:]
    retained = [r for r in rows if int(r["type"]) != 2] + selected
    return {
        "bounds_angstrom": bounds,
        "periodic": [True, True, True],
        "atoms": [{"id": int(r["id"]), "species": species[int(r["type"])],
                   "position_angstrom": [float(r[k]) for k in ("x", "y", "z")]}
                  for r in retained],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--structure", default="SCR_PureFe_allH", help="Paper-facing ID")
    parser.add_argument("--atom-id", type=int, help="Source H atom ID; default first site")
    args = parser.parse_args()
    manifest = read_json(ROOT / "generated_metadata" / "manifest_index.json")
    entry = next((e for e in manifest["structures"]
                  if e["paper_structure_id"] == args.structure), None)
    if entry is None:
        parser.error("Unknown structure; consult generated_metadata/manifest_index.json")
    collection = read_json(ROOT / "generated_metadata" / entry["site_file"])
    site = next((s for s in collection["sites"]
                 if args.atom_id is None or s["source_key"]["atom_id"] == args.atom_id), None)
    if site is None:
        parser.error("No H site has this atom ID in the selected structure")
    sample = reconstruct_one_h(ROOT / entry["source_structure_id"],
                               site["source_key"]["atom_id"])
    assert sum(a["species"] == "H" for a in sample["atoms"]) == 1
    assert len(sample["atoms"]) == entry["host_atom_count"] + 1
    print(f"Dataset: {manifest['structure_count']} hosts, {manifest['site_count']} H sites")
    print(f"Selected: {site['site_id']}; initial atoms: {len(sample['atoms'])}")
    for branch in ("frozen", "non_frozen"):
        p = site["calculations"][branch]["properties"]
        print(f"{branch}: SE={p['solution_energy']['value']} eV; "
              f"SE_ZPEC={p['solution_energy_zpe_corrected']['value']} eV")
    print("None means missing; review each property's provenance and conflict flags.")


if __name__ == "__main__":
    main()
