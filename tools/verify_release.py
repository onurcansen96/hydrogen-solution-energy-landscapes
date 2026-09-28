"""Read-only checksum, site-key, and single-H reconstruction checks (stdlib only)."""

import argparse
import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "examples"))
from read_dataset import reconstruct_one_h


def load(path):
    def invalid(value):
        raise ValueError(f"Non-finite JSON token in {path}: {value}")
    return json.loads(path.read_text(encoding="utf-8"), parse_constant=invalid)


def digest(path):
    result = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def check_file(path, record):
    if not path.is_file():
        raise ValueError(f"Missing file: {path}")
    if path.stat().st_size != record["byte_size"] or digest(path) != record["sha256"]:
        raise ValueError(f"Size/checksum mismatch: {path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chgcar-dir", type=Path)
    args = parser.parse_args()
    meta = ROOT / "generated_metadata"
    manifest = load(meta / "manifest_index.json")
    source = load(meta / "source_manifest.json")
    checked = 0
    for item in source["artifacts"]:
        if item["included_in_package"]:
            check_file(ROOT / item["path"], item)
            checked += 1
    for item in manifest["output_artifacts"]:
        check_file(meta / item["path"], item)
        checked += 1
    keys = set()
    family_counts = {}
    for entry in manifest["structures"]:
        collection = load(meta / entry["site_file"])
        sites = collection["sites"]
        if len(sites) != entry["hydrogen_site_count"]:
            raise ValueError(f"Site count mismatch: {entry['paper_structure_id']}")
        for site in sites:
            key = (site["source_key"]["structure"], site["source_key"]["atom_id"])
            if key in keys:
                raise ValueError(f"Duplicate site key: {key}")
            keys.add(key)
            sample = reconstruct_one_h(ROOT / key[0], key[1])
            if sum(a["species"] == "H" for a in sample["atoms"]) != 1:
                raise ValueError(f"Single-H constraint failed: {key}")
            if len(sample["atoms"]) != entry["host_atom_count"] + 1:
                raise ValueError(f"Atom count mismatch: {key}")
            if sample["atoms"][-1]["position_angstrom"] != site["candidate_atom"]["position"]["value"]:
                raise ValueError(f"H coordinate mismatch: {key}")
        family = entry["raw_family"]
        family_counts[family] = family_counts.get(family, 0) + len(sites)
    if len(keys) != 1648 or len(manifest["structures"]) != 25:
        raise ValueError("Expected 25 hosts and 1,648 sites in this release")
    if family_counts != {"SC": 31, "S3": 113, "S5": 425, "S27": 471, "S5m": 608}:
        raise ValueError(f"Family count mismatch: {family_counts}")
    catalog = load(ROOT / "charge_density_manifest.json")
    if len(catalog["files"]) != 20:
        raise ValueError("Expected 20 charge-density resources")
    resource_ref = source.get("charge_density_resources")
    if resource_ref and digest(ROOT / resource_ref["catalog_path"]) != resource_ref["catalog_sha256"]:
        raise ValueError("Charge-density catalog checksum mismatch")
    if args.chgcar_dir:
        for item in catalog["files"]:
            check_file(args.chgcar_dir / item["filename"], item)
    print(json.dumps({"status": "passed", "checked_source_and_output_files": checked,
                      "host_count": 25, "site_count": len(keys),
                      "one_H_reconstructions_checked": len(keys),
                      "family_site_counts": family_counts,
                      "CHGCAR_files_verified": 20 if args.chgcar_dir else 0,
                      "CHGCAR_note": "Pass --chgcar-dir to verify downloaded raw grids"}, indent=2))


if __name__ == "__main__":
    main()
