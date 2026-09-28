#!/usr/bin/env python3
"""Build the hydrogen solution-energy metadata release.

The source ``*_allH`` files are aggregate candidate-site maps.  They are not
physical many-H configurations.  This generator therefore writes:

* conceptual-dictionary-compatible JSON containing only the H-free hosts;
* one lossless site collection per structure, with one record per alternative
  one-H calculation;
* shared workflow, field-dictionary, source-manifest and validation files.

The script deliberately preserves reported values and reports source conflicts
instead of silently resolving them.  Dump cells and coordinates are
authoritative.  ZPE-corrected energies are imported from the descriptor table.
"""

from __future__ import annotations

import argparse
import ast
import csv
import hashlib
import json
import math
import os
import shutil
import statistics
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import numpy as np


GENERATOR_NAME = "convert_h_solution_landscape.py"
GENERATOR_VERSION = "1.1.1"
SCHEMA_VERSION = "1.1.0"
EXTENSION_KEY = "x_hydrogen_landscape"
EXTENSION_NAMESPACE = "urn:ocdo:hydrogen-solution-landscape:1"

DESCRIPTOR_FILENAME = "Solution_Energy_Descriptors.txt"
SELECTED_DESCRIPTOR_FILENAME = "H_Voronoi_selected_descriptors_4SE.csv"
FULL_DESCRIPTOR_FILENAME = "H_Voronoi_descriptors.csv"

HALF_H2_ZPE_EV = 0.1334825
H2_ZPE_EV = 2.0 * HALF_H2_ZPE_EV
H2_ELECTRONIC_ENERGY_EV = -6.759828

ENERGY_TOLERANCE_EV = 1.0e-6
COORDINATE_TOLERANCE_ANGSTROM = 1.0e-6
SCALAR_SOURCE_TOLERANCE = 1.0e-6
DERIVED_DESCRIPTOR_TOLERANCE = 1.0e-9

CREATORS = [
    {"id": None, "name": "Onur Can Şen"},
    {"id": None, "name": "Farnoosh Ghaderi"},
    {"id": None, "name": "Santiago Benito"},
    {"id": None, "name": "Sebastian Weber"},
    {"id": None, "name": "Rebecca Janisch"},
]

CHEMISTRIES = ("PureFe", "FeCr", "FeCu", "FeTi", "FeVac")
RAW_FAMILIES = ("SC", "S3", "S5", "S27", "S5m")

FAMILY_METADATA: dict[str, dict[str, Any]] = {
    "SC": {
        "paper_abbreviation": "SCR",
        "paper_name": "bar(3)10-oriented rotated single crystal",
        "structure_class": "rotated_single_crystal",
        "paper_initial_cell_angstrom": [26.87, 5.66, 8.95],
        "paper_base_host_atoms": 120,
        "k_mesh": [1, 6, 4],
        "segregation_optimization": "single_site",
        "normal_volume_optimization": True,
        "rigid_body_translation": False,
        "orientation_plane": [-3, 1, 0],
        "cell_orientation": {
            "x_parallel": [-3, 1, 0],
            "y_parallel": [0, 0, 1],
            "z_parallel": [1, 3, 0],
        },
        "orientation_rotation": {
            "axis": [0, 0, 1],
            "reference_direction": [-1, 0, 0],
            "target_direction": [-3, 1, 0],
            "angle_magnitude_degrees": 18.43494882292201,
            "signed_angle_degrees": -18.43494882292201,
            "right_hand_rule_axis": [0, 0, 1],
            "status": "derived_from_orientation_vectors",
            "derivation": "angle magnitude = atan(1/3); signed angle is clockwise from [-1,0,0] to [-3,1,0] about +[0,0,1]",
            "not_a_grain_boundary_misorientation": True,
        },
        "grain_boundary": None,
    },
    "S3": {
        "paper_abbreviation": "S3",
        "paper_name": "Sigma 3 symmetric-tilt grain boundary",
        "structure_class": "symmetric_tilt_grain_boundary",
        "paper_initial_cell_angstrom": [27.75, 8.01, 4.90],
        "paper_base_host_atoms": 96,
        "k_mesh": [1, 4, 6],
        "segregation_optimization": "all_sites",
        "normal_volume_optimization": True,
        "rigid_body_translation": True,
        "grain_boundary": {
            "native_key": "symmetric_tilt_grain_boundary",
            "sigma": 3,
            "plane": [1, 1, -2],
            "rotation_axis": [1, -1, 0],
            "misorientation_angle_degrees": 70.52877936550931,
            "misorientation_angle_status": "derived_ideal_CSL_axis_angle",
            "misorientation_angle_derivation": "2*atan(sqrt(2)/2)",
            "reference_energy_j_per_m2": 0.41,
        },
    },
    "S5": {
        "paper_abbreviation": "S5",
        "paper_name": "Sigma 5 symmetric-tilt grain boundary",
        "structure_class": "symmetric_tilt_grain_boundary",
        "paper_initial_cell_angstrom": [26.87, 5.66, 8.95],
        "paper_base_host_atoms": 120,
        "k_mesh": [1, 6, 4],
        "segregation_optimization": "all_sites",
        "normal_volume_optimization": True,
        "rigid_body_translation": True,
        "grain_boundary": {
            "native_key": "symmetric_tilt_grain_boundary",
            "sigma": 5,
            "plane": [-3, 1, 0],
            "rotation_axis": [0, 0, 1],
            "misorientation_angle_degrees": 36.86989764584402,
            "misorientation_angle_status": "derived_ideal_CSL_axis_angle",
            "misorientation_angle_derivation": "2*atan(1/3)",
            "reference_energy_j_per_m2": 1.54,
        },
    },
    "S27": {
        "paper_abbreviation": "S27",
        "paper_name": "Sigma 27 symmetric-tilt grain boundary",
        "structure_class": "symmetric_tilt_grain_boundary",
        "paper_initial_cell_angstrom": [41.63, 8.01, 7.36],
        "paper_base_host_atoms": 216,
        "k_mesh": [1, 4, 4],
        "segregation_optimization": "all_sites",
        "normal_volume_optimization": True,
        "rigid_body_translation": False,
        "grain_boundary": {
            "native_key": "symmetric_tilt_grain_boundary",
            "sigma": 27,
            "plane": [5, -5, 2],
            "rotation_axis": [1, 1, 0],
            "misorientation_angle_degrees": 31.586338096527925,
            "misorientation_angle_status": "derived_ideal_CSL_axis_angle",
            "misorientation_angle_derivation": "2*atan(sqrt(2)/5)",
            "reference_energy_j_per_m2": 1.68,
        },
    },
    "S5m": {
        "paper_abbreviation": "S5m",
        "paper_name": "mixed near-CSL Sigma 5m grain boundary",
        "structure_class": "mixed_grain_boundary",
        "paper_initial_cell_angstrom": [20.81, 6.33, 15.51],
        "paper_base_host_atoms": 180,
        "k_mesh": [1, 3, 1],
        "segregation_optimization": "single_site",
        "normal_volume_optimization": True,
        "rigid_body_translation": False,
        "grain_boundary": {
            "native_key": "mixed_grain_boundary",
            "sigma": 5,
            "plane": [1, 1, -2],
            "rotation_axis": [1, 0, 0],
            "misorientation_angle_degrees": 36.86989764584402,
            "misorientation_angle_status": "derived_ideal_CSL_axis_angle",
            "misorientation_angle_derivation": "2*atan(1/3)",
            "tilt_angle_degrees": 33.6,
            "tilt_angle_status": "reported_in_manuscript",
            "twist_angle_degrees": 0.27,
            "twist_angle_status": "reported_in_manuscript",
            "reference_energy_j_per_m2": 2.09,
        },
    },
}

ELEMENT_INITIAL_MAGNETIC_MOMENT_MU_B = {
    "Fe": 3.0,
    "Cr": -3.0,
    "H": 0.0,
    "Ti": 0.0,
    "Cu": 0.0,
}

# These keys are NA in the upstream S27 NF Bader summary but were converted to
# literal 0.0 in the dump/descriptor table.  They are restored to JSON null.
BC_NF_ZERO_PLACEHOLDER_KEYS = {
    *(('S27_PureFe_allH', i) for i in (235, 240, 242, 260, 276, 280, 286, 300)),
    *(('S27_FeCr_allH', i) for i in (227, 229, 232, 237, 239, 250, 260)),
    *(('S27_FeCu_allH', i) for i in (226, 232, 234, 237, 238, 255, 267, 270, 276, 281, 289)),
    *(('S27_FeTi_allH', i) for i in (218, 221, 223, 227, 230, 232, 236, 240, 248, 257, 278, 279, 281, 282)),
    *(('S27_FeVac_allH', i) for i in (217, 229, 236, 247, 248, 250, 253, 254, 264, 272, 280, 302)),
}

EXPECTED_MISSING_FROZEN_BC_KEYS = {
    ("S27_FeCr_allH", 226),
    ("S27_FeTi_allH", 238),
    ("S27_FeTi_allH", 279),
    ("S27_FeVac_allH", 252),
    ("S27_FeVac_allH", 292),
    ("S27_FeVac_allH", 297),
}

EXPECTED_MISSING_NF_ZPE_KEYS = {
    ("S5_FeTi_allH", 261),
    ("S27_FeTi_allH", 257),
}

# At these six sites the dump reports a non-zero NF ZPE, whereas the descriptor
# table stores zero and its corrected energy was calculated with that zero.
SC_FECR_NF_ZPE_CONFLICT_KEYS = {
    ("SC_FeCr_allH", i) for i in range(120, 126)
}

DUMP_PROPERTY_TO_DESCRIPTOR = {
    "SE": "Atom_SE",
    "CD": "Atom_CD",
    "SE_NF": "Atom_SE_NF",
    "BC": "Atom_BC",
    "ZPE": "Atom_ZPE",
    "ZPE_NF": "Atom_ZPE_NF",
    "BC_NF": "Atom_BC_NF",
}


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def parse_optional_float(value: Any) -> float | None:
    if value is None:
        return None
    text = str(value).strip()
    if text == "" or text.lower() in {"nan", "na", "none", "null"}:
        return None
    number = float(text)
    if not math.isfinite(number):
        return None
    return number


def json_number(value: Any) -> int | float | None:
    if value is None:
        return None
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        value = float(value)
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def write_json(path: Path, data: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(data, handle, indent=2, ensure_ascii=False, allow_nan=False)
        handle.write("\n")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def portable_relpath(path: Path, start: Path) -> str:
    return Path(os.path.relpath(path.resolve(), start.resolve())).as_posix()


def artifact_record(
    artifact_id: str,
    path: Path,
    path_base: Path,
    role: str,
    included_in_package: bool,
    media_type: str = "text/plain",
) -> dict[str, Any]:
    return {
        "artifact_id": artifact_id,
        "path": portable_relpath(path, path_base) if included_in_package else None,
        "filename": path.name,
        "role": role,
        "media_type": media_type,
        "included_in_package": included_in_package,
        "location_status": (
            "relative_to_canonical_source_root"
            if included_in_package
            else "external_reference_not_included"
        ),
        "byte_size": path.stat().st_size,
        "sha256": sha256_file(path),
    }


def split_structure_name(source_structure: str) -> tuple[str, str]:
    parts = source_structure.split("_")
    if len(parts) < 3 or parts[-1] != "allH":
        raise ValueError(f"Unexpected structure name: {source_structure}")
    family, chemistry = parts[0], parts[1]
    if family not in RAW_FAMILIES or chemistry not in CHEMISTRIES:
        raise ValueError(f"Unsupported family/chemistry: {source_structure}")
    return family, chemistry


def paper_structure_name(source_structure: str) -> str:
    family, chemistry = split_structure_name(source_structure)
    paper_family = FAMILY_METADATA[family]["paper_abbreviation"]
    return f"{paper_family}_{chemistry}_allH"


def type_species_map(family: str, chemistry: str) -> dict[int, str]:
    mapping = {1: "Fe", 2: "H"}
    if chemistry in {"FeCr", "FeCu", "FeTi"}:
        mapping[3] = chemistry[2:]
    return mapping


def build_angle_metadata(family: str) -> dict[str, Any]:
    """Return explicit family-level angle metadata with source status."""
    metadata = FAMILY_METADATA[family]
    gb = metadata.get("grain_boundary")
    if gb is None:
        rotation = metadata["orientation_rotation"]
        return {
            "angle_scope": "family-level rotated-single-crystal orientation",
            "angle_type": "single_crystal_orientation_rotation",
            "orientation_rotation": {
                "angle_magnitude": {
                    "value": rotation["angle_magnitude_degrees"],
                    "unit": "degree",
                },
                "signed_angle": {
                    "value": rotation["signed_angle_degrees"],
                    "unit": "degree",
                },
                "axis": list(rotation["axis"]),
                "reference_direction": list(rotation["reference_direction"]),
                "target_direction": list(rotation["target_direction"]),
                "right_hand_rule_axis": list(rotation["right_hand_rule_axis"]),
                "status": rotation["status"],
                "derivation": rotation["derivation"],
                "not_a_grain_boundary_misorientation": True,
            },
        }

    result: dict[str, Any] = {
        "angle_scope": "family-level grain-boundary crystallography; identical for every chemistry and H site in the family",
        "angle_type": "grain_boundary_misorientation",
        "misorientation": {
            "value": gb["misorientation_angle_degrees"],
            "unit": "degree",
            "rotation_axis": list(gb["rotation_axis"]),
            "boundary_plane": list(gb["plane"]),
            "status": gb["misorientation_angle_status"],
            "derivation": gb["misorientation_angle_derivation"],
            "angle_convention": "lower-angle representative in [0,90] degrees for the listed ideal cubic CSL plane/axis construction; symmetry-equivalent complementary conventions may exist",
            "definition": "Ideal cubic CSL axis-angle for the listed Sigma representation; the manuscript specifies the plane/axis but does not print this number.",
        },
    }
    if gb["native_key"] == "mixed_grain_boundary":
        result["mixed_character_components"] = {
            "tilt": {
                "value": gb["tilt_angle_degrees"],
                "unit": "degree",
                "status": gb["tilt_angle_status"],
            },
            "twist": {
                "value": gb["twist_angle_degrees"],
                "unit": "degree",
                "status": gb["twist_angle_status"],
            },
            "interpretation": "Reported mixed-boundary character components; values and degree units are preserved exactly as printed, without reinterpretation. Do not add them together or substitute their sum for the CSL axis-angle.",
        }
    return result


def read_dump(path: Path) -> dict[str, Any]:
    lines = path.read_text(encoding="utf-8").splitlines()
    if len(lines) < 10 or lines[0] != "ITEM: TIMESTEP":
        raise ValueError(f"{path.name}: not a supported LAMMPS dump")
    timestep = int(lines[1].strip())
    if lines[2] != "ITEM: NUMBER OF ATOMS":
        raise ValueError(f"{path.name}: missing atom-count header")
    declared_count = int(lines[3].strip())
    box_header = lines[4]
    if not box_header.startswith("ITEM: BOX BOUNDS"):
        raise ValueError(f"{path.name}: missing box-bounds header")
    bounds = []
    for line in lines[5:8]:
        values = [float(v) for v in line.split()]
        if len(values) < 2:
            raise ValueError(f"{path.name}: malformed box bound")
        bounds.append(values[:2])
    atoms_header = lines[8]
    if not atoms_header.startswith("ITEM: ATOMS "):
        raise ValueError(f"{path.name}: missing atom header")
    columns = atoms_header[len("ITEM: ATOMS ") :].split()
    expected_columns = ["type", "id", "x", "y", "z", "SE", "CD", "SE_NF", "BC", "ZPE", "ZPE_NF", "BC_NF"]
    if columns != expected_columns:
        raise ValueError(f"{path.name}: unexpected columns {columns}")
    atom_lines = lines[9:]
    if len(atom_lines) != declared_count:
        raise ValueError(
            f"{path.name}: declared {declared_count} atoms but found {len(atom_lines)} rows"
        )
    atoms = []
    for row_number, line in enumerate(atom_lines, start=10):
        tokens = line.split()
        if len(tokens) != len(columns):
            raise ValueError(f"{path.name}:{row_number}: expected {len(columns)} fields")
        raw = dict(zip(columns, tokens))
        atom = {
            "type": int(raw["type"]),
            "id": int(raw["id"]),
            "position": [float(raw["x"]), float(raw["y"]), float(raw["z"])],
            "properties": {name: parse_optional_float(raw[name]) for name in expected_columns[5:]},
            "raw_properties": {name: raw[name] for name in expected_columns[5:]},
        }
        atoms.append(atom)
    ids = [atom["id"] for atom in atoms]
    if len(ids) != len(set(ids)):
        raise ValueError(f"{path.name}: duplicate atom IDs")
    if sorted(ids) != list(range(1, declared_count + 1)):
        raise ValueError(f"{path.name}: atom IDs are not contiguous 1..N")
    lengths = [upper - lower for lower, upper in bounds]
    return {
        "path": path,
        "name": path.name,
        "timestep": timestep,
        "declared_count": declared_count,
        "box_header": box_header,
        "bounds": bounds,
        "lengths": lengths,
        "columns": columns,
        "atoms": atoms,
    }


def read_descriptor_table(path: Path) -> tuple[list[dict[str, str]], dict[tuple[str, int], dict[str, str]]]:
    csv.field_size_limit(min(sys.maxsize, 2**31 - 1))
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    keyed: dict[tuple[str, int], dict[str, str]] = {}
    for row in rows:
        key = (row["Structure"], int(row["Atom_id"]))
        if key in keyed:
            raise ValueError(f"Duplicate descriptor key: {key}")
        keyed[key] = row
    return rows, keyed


def read_selected_descriptor_reference(path: Path | None) -> dict[tuple[str, int], dict[str, str]]:
    if path is None or not path.exists():
        return {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    keyed: dict[tuple[str, int], dict[str, str]] = {}
    for row in rows:
        key = (row["Structure"], int(row["Index"]))
        if key in keyed:
            raise ValueError(f"Duplicate selected-descriptor key: {key}")
        keyed[key] = row
    return keyed


def parse_literal_list(row: dict[str, str], field: str) -> list[Any]:
    try:
        value = ast.literal_eval(row[field])
    except (SyntaxError, ValueError) as exc:
        raise ValueError(f"Cannot parse {field} for descriptor index {row.get('Index')}") from exc
    if not isinstance(value, list):
        raise ValueError(f"{field} is not a list for descriptor index {row.get('Index')}")
    return value


def recompute_geometry(row: dict[str, str]) -> tuple[dict[str, Any], dict[str, Any]]:
    vertices_list = parse_literal_list(row, "Vertices")
    edges = parse_literal_list(row, "Edges")
    faces = parse_literal_list(row, "Faces")
    vertices = np.asarray(vertices_list, dtype=float)
    if vertices.ndim != 2 or vertices.shape[1] != 3 or len(vertices) < 4:
        raise ValueError(f"Invalid vertex array at descriptor index {row['Index']}")
    for face in faces:
        if len(face) != 3 or any((not isinstance(i, int) or i < 0 or i >= len(vertices)) for i in face):
            raise ValueError(f"Invalid triangular face at descriptor index {row['Index']}: {face}")
    for edge in edges:
        if len(edge) != 2 or any((not isinstance(i, int) or i < 0 or i >= len(vertices)) for i in edge):
            raise ValueError(f"Invalid edge at descriptor index {row['Index']}: {edge}")

    area = 0.0
    for i, j, k in faces:
        area += float(np.linalg.norm(np.cross(vertices[j] - vertices[i], vertices[k] - vertices[i])) / 2.0)

    centroid = vertices.mean(axis=0)
    centered = vertices - centroid
    covariance = np.cov(centered, rowvar=False)
    eigenvalues_ascending = np.linalg.eigvalsh(covariance)
    if eigenvalues_ascending[0] <= 0.0:
        raise ValueError(f"Non-positive PCA eigenvalue at descriptor index {row['Index']}")
    lambda3, lambda2, lambda1 = [float(v) for v in eigenvalues_ascending]
    radii = np.linalg.norm(centered, axis=1)
    r_in = float(radii.min())
    r_out = float(radii.max())
    if r_in <= 0.0:
        raise ValueError(f"Zero centroid-to-vertex radius at descriptor index {row['Index']}")

    volume = float(row["Voronoi_Volume"])
    if volume <= 0.0 or area <= 0.0:
        raise ValueError(f"Non-positive geometry at descriptor index {row['Index']}")
    derived = {
        "surface_area": area,
        "volume_cuberoot": volume ** (1.0 / 3.0),
        "area_over_volume": area / volume,
        "sphericity": math.pi ** (1.0 / 3.0) * (6.0 * volume) ** (2.0 / 3.0) / area,
        "elongation_lambda1_over_lambda3": lambda1 / lambda3,
        "anisotropy_lambda2_over_lambda3": lambda2 / lambda3,
        "rout_over_rin": r_out / r_in,
    }
    details = {
        "vertex_centroid": [float(v) for v in centroid],
        "pca_eigenvalues_descending": [lambda1, lambda2, lambda3],
        "rout": r_out,
        "rin": r_in,
    }
    topology = {
        "vertices": vertices_list,
        "edges": edges,
        "faces": faces,
    }
    return {**derived, **details}, topology


def make_quantity(
    value: float | int | None,
    unit: str | None,
    status: str,
    source_artifact: str | None,
    source_column: str | None,
    definition: str,
    **extra: Any,
) -> dict[str, Any]:
    quantity = {
        "value": json_number(value),
        "unit": unit,
        "status": status,
        "definition": definition,
        "source": {
            "artifact_id": source_artifact,
            "column": source_column,
        },
    }
    quantity.update(extra)
    return quantity


def property_status(value: Any, missing_status: str = "not_calculated") -> str:
    return "reported" if value is not None else missing_status


def build_computational_settings() -> dict[str, Any]:
    return {
        "schema_version": SCHEMA_VERSION,
        "method": "DensityFunctionalTheory",
        "software": {"name": "VASP", "version": None},
        "electronic_structure": {
            "core_valence_method": "PAW",
            "exchange_correlation": "GGA-PBE",
            "spin_polarized": True,
            "equilibrium_bcc_fe_lattice_constant": {"value": 2.833, "unit": "angstrom"},
        },
        "effective_incar": {
            "ISTART": 0,
            "ISIF": 2,
            "ENCUT": 400,
            "ALGO": "Fast",
            "EDIFF": 1.0e-6,
            "ISPIN": 2,
            "EDIFFG": -0.01,
            "NSW": 1000,
            "ISMEAR": 1,
            "SIGMA": 0.1,
            "IBRION": 2,
            "POTIM": 0.2,
            "LWAVE": False,
            "LCHARG": True,
            "LVTOT": True,
            "LELF": True,
            "LORBIT": 10,
            "NBANDS": None,
        },
        "incar_units": {
            "ENCUT": "eV",
            "EDIFF": "eV",
            "EDIFFG": "eV/angstrom",
            "SIGMA": "eV",
        },
        "author_confirmed_initial_magnetic_moments": {
            element: {"value": value, "unit": "mu_B"}
            for element, value in ELEMENT_INITIAL_MAGNETIC_MOMENT_MU_B.items()
        },
        "magmom_note": "Initial collinear moments, not constrained or converged local moments; repetition counts follow each reconstructed one-H composition and atom ordering.",
        "relaxation_protocols": {
            "frozen": {
                "initial_structure": "H-free dump host plus exactly one selected candidate H",
                "cell": "fixed",
                "movable_atoms": "selected H only",
                "fixed_atoms": "all non-H atoms at their H-free pre-insertion dump positions",
                "selective_dynamics_required": True,
            },
            "non_frozen": {
                "initial_structure": "same reconstructed one-H structure as frozen branch",
                "cell": "fixed",
                "movable_atoms": "all atoms",
                "fixed_atoms": "none",
                "selective_dynamics_required": False,
            },
        },
        "solution_energy": {
            "electronic_expression": "E_host_plus_H - E_host - 0.5 * E_H2",
            "h2_electronic_energy": {"value": H2_ELECTRONIC_ENERGY_EV, "unit": "eV", "provenance": "upstream SE_maker notebook"},
            "half_h2_electronic_energy": {"value": H2_ELECTRONIC_ENERGY_EV / 2.0, "unit": "eV"},
        },
        "zero_point_energy": {
            "expression": "0.5 * sum_i(h * nu_i)",
            "site_protocol": "harmonic finite differences; displace only H while holding all host atoms fixed; use three local H modes",
            "h2_protocol": "one H2 stretching mode",
            "zpe_corrected_solution_energy_expression": "SE + ZPE_site - 0.5 * ZPE_H2",
            "h2_zpe_reference": {"value": H2_ZPE_EV, "unit": "eV"},
            "half_h2_zpe_reference": {"value": HALF_H2_ZPE_EV, "unit": "eV"},
            "reference_role": "dataset numerical reference verified from executable source and stored corrected values; not printed numerically in the manuscript",
        },
        "post_processing": {
            "bader": {
                "quantity": "H Bader charge-transfer value",
                "unit": "e",
                "transformation": "preserved exactly as reported; no sign reversal or subtraction",
                "software_version": None,
            },
            "charge_density": {
                "expression": "total_charge_in_candidate_H_Voronoi_region / candidate_H_Voronoi_volume",
                "unit": "e/Å³",
                "unit_ascii": "e/angstrom^3",
                "geometry_state": "H-free pre-insertion host",
            },
        },
        "family_settings": {
            FAMILY_METADATA[family]["paper_abbreviation"]: {
                "raw_family": family,
                "paper_initial_cell_angstrom": FAMILY_METADATA[family]["paper_initial_cell_angstrom"],
                "k_mesh": FAMILY_METADATA[family]["k_mesh"],
                "segregation_optimization": FAMILY_METADATA[family]["segregation_optimization"],
                "normal_volume_optimization": FAMILY_METADATA[family]["normal_volume_optimization"],
                "rigid_body_translation": FAMILY_METADATA[family]["rigid_body_translation"],
                "angle_metadata": build_angle_metadata(family),
            }
            for family in RAW_FAMILIES
        },
        "known_provenance_limitations": {
            "vasp_version": None,
            "potcar_labels_and_releases": None,
            "k_point_centering": None,
            "finite_displacement_amplitude": None,
            "isolated_h2_cell_and_bond_length": None,
            "bader_software_and_grid_settings": None,
            "relaxed_h_containing_output_coordinates_in_release": False,
        },
        "ediff_provenance_note": "The supplied representative INCAR contains 1e-5 eV, while the manuscript and the author's explicit correction specify 1e-6 eV. This release uses 1e-6 eV as authoritative.",
    }


def build_workflow(
    host_sample_id: str,
    state: str,
    computational_settings_ref: str,
) -> dict[str, Any]:
    if state == "frozen":
        movable = "selected H only"
        fixed = "all non-H atoms"
    elif state == "non_frozen":
        movable = "all atoms"
        fixed = "none"
    else:
        raise ValueError(state)
    return {
        "id": f"workflow_{host_sample_id}_{state}",
        "algorithm": None,
        "method": "DensityFunctionalTheory",
        "xc_functional": "PBE",
        "input_parameter": [],
        "input_sample": [host_sample_id],
        "output_sample": [],
        "output_parameter": [],
        "calculated_property": [],
        "degrees_of_freedom": ["AtomicPositionRelaxation"],
        "interatomic_potential": {"potential_type": None, "uri": None},
        "software": [{"uri": "https://www.vasp.at/", "version": None, "label": "VASP"}],
        "workflow_manager": {"uri": None, "version": None, "label": None},
        "thermodynamic_ensemble": None,
        EXTENSION_KEY: {
            "namespace": EXTENSION_NAMESPACE,
            "protocol_role": f"one-H {state.replace('_', '-')} relaxation",
            "cell": "fixed",
            "movable_atoms": movable,
            "fixed_atoms": fixed,
            "output_geometry_available": False,
            "computational_settings_ref": computational_settings_ref,
            "site_results_location": "per-structure hydrogen-site sidecar",
        },
    }


def expected_defect_count(family: str, chemistry: str) -> int:
    if chemistry == "PureFe":
        return 0
    return 1 if family == "SC" else 2


def build_host_sample(dump: dict[str, Any], family: str, chemistry: str) -> tuple[dict[str, Any], dict[str, Any]]:
    metadata = FAMILY_METADATA[family]
    species_map = type_species_map(family, chemistry)
    host_atoms = [atom for atom in dump["atoms"] if atom["type"] != 2]
    host_species = [species_map[atom["type"]] for atom in host_atoms]
    host_positions = [atom["position"] for atom in host_atoms]
    host_ids = [atom["id"] for atom in host_atoms]
    host_types = [atom["type"] for atom in host_atoms]
    composition = Counter(host_species)
    host_count = len(host_atoms)
    element_ratio = {element: count / host_count for element, count in sorted(composition.items())}

    source_name = dump["name"]
    paper_name = paper_structure_name(source_name)
    paper_family = metadata["paper_abbreviation"]
    host_sample_id = f"host_{paper_family}_{chemistry}"
    volume = math.prod(dump["lengths"])
    calculated_property = []
    gb = metadata.get("grain_boundary")
    if gb is not None:
        calculated_property.append(
            {
                "id": f"gb_energy_{paper_family}_{chemistry}",
                "label": "Paper reference grain-boundary energy",
                "basename": "GrainBoundaryEnergy",
                "value": gb["reference_energy_j_per_m2"],
                "unit": "J-PER-M2",
                "associate_to_sample": [host_sample_id],
            }
        )

    sample: dict[str, Any] = {
        "id": host_sample_id,
        "material": {
            "element_ratio": element_ratio,
            "crystal_structure": {
                "spacegroup_symbol": "Im-3m",
                "spacegroup_number": 229,
                "unit_cell": {
                    "bravais_lattice": "https://www.wikidata.org/wiki/Q851536",
                    "lattice_parameter": [2.833, 2.833, 2.833],
                    "angle": [90.0, 90.0, 90.0],
                },
            },
        },
        "simulation_cell": {
            "volume": {"value": volume},
            "number_of_atoms": host_count,
            "length": dump["lengths"],
            "vector": [
                [dump["lengths"][0], 0.0, 0.0],
                [0.0, dump["lengths"][1], 0.0],
                [0.0, 0.0, dump["lengths"][2]],
            ],
            "angle": [90.0, 90.0, 90.0],
            "repetitions": [],
            "grain_size": None,
            "number_of_grains": 1 if family == "SC" else None,
        },
        "atom_attribute": {
            "position": host_positions,
            "species": host_species,
            "file_path": None,
            "file_format": None,
            "file_species": None,
        },
        "calculated_property": calculated_property,
    }

    defect_count = expected_defect_count(family, chemistry)
    if chemistry in {"FeCr", "FeCu", "FeTi"}:
        sample["substitutional"] = {
            "concentration": defect_count / host_count,
            "number": defect_count,
        }
    elif chemistry == "FeVac":
        reference_lattice_sites = host_count + defect_count
        sample["vacancy"] = {
            "concentration": defect_count / reference_lattice_sites,
            "number": defect_count,
        }

    if gb is not None:
        sample[gb["native_key"]] = {
            "sigma": gb["sigma"],
            "plane": gb["plane"],
            "misorientation_angle": gb["misorientation_angle_degrees"],
            "rotation_axis": gb["rotation_axis"],
        }

    extension = {
        "namespace": EXTENSION_NAMESPACE,
        "schema_version": SCHEMA_VERSION,
        "source_structure_id": source_name,
        "paper_structure_id": paper_name,
        "raw_family": family,
        "paper_abbreviation": paper_family,
        "chemistry": chemistry,
        "structure_class": metadata["structure_class"],
        "host_definition": "all source dump atoms whose LAMMPS type is not 2",
        "host_atom_ids": host_ids,
        "host_source_types": host_types,
        "host_species_counts": dict(sorted(composition.items())),
        "type_species_map": {str(key): value for key, value in species_map.items()},
        "source_box_bounds": dump["bounds"],
        "source_box_bounds_header": dump["box_header"],
        "periodic_boundary_conditions": [True, True, True],
        "source_coordinates_preserved_unwrapped": True,
        "paper_initial_cell_angstrom": metadata["paper_initial_cell_angstrom"],
        "actual_dump_cell_is_authoritative": True,
        "k_mesh": metadata["k_mesh"],
        "segregation_optimization": metadata["segregation_optimization"],
        "normal_volume_optimization": metadata["normal_volume_optimization"],
        "rigid_body_translation": metadata["rigid_body_translation"],
        "solute_species": chemistry[2:] if chemistry in {"FeCr", "FeCu", "FeTi"} else None,
        "solute_or_vacancy_count": defect_count,
        "reference_grain_boundary_energy_scope": "family-level paper reference; not independently recomputed for every chemistry or H site" if gb is not None else None,
        "orientation_plane": metadata.get("orientation_plane"),
        "cell_orientation": metadata.get("cell_orientation"),
        "angle_metadata": build_angle_metadata(family),
    }
    sample[EXTENSION_KEY] = extension
    return sample, extension


def quantity_summary(values: Iterable[float | None]) -> dict[str, Any]:
    materialized = list(values)
    complete = [float(value) for value in materialized if value is not None]
    return {
        "count": len(complete),
        "missing": sum(value is None for value in materialized),
        "minimum": min(complete) if complete else None,
        "maximum": max(complete) if complete else None,
        "mean": statistics.fmean(complete) if complete else None,
    }


def build_site_record(
    dump: dict[str, Any],
    h_atom: dict[str, Any],
    descriptor_row: dict[str, str],
    family: str,
    chemistry: str,
    site_sequence: int,
    derived_geometry: dict[str, Any],
    topology: dict[str, Any],
    selected_reference_row: dict[str, str] | None,
) -> tuple[dict[str, Any], dict[str, Any]]:
    source_structure = dump["name"]
    paper_structure = paper_structure_name(source_structure)
    paper_family = FAMILY_METADATA[family]["paper_abbreviation"]
    host_sample_id = f"host_{paper_family}_{chemistry}"
    atom_id = h_atom["id"]
    key = (source_structure, atom_id)
    site_id = f"site_{paper_family}_{chemistry}_H_{atom_id}"
    initial_sample_id = f"initial_{paper_family}_{chemistry}_H_{atom_id}"
    dump_artifact_id = f"dump_{source_structure}"
    descriptor_artifact_id = "descriptor_table"

    se_f = parse_optional_float(descriptor_row["Atom_SE"])
    se_nf = parse_optional_float(descriptor_row["Atom_SE_NF"])
    zpe_f = parse_optional_float(descriptor_row["Atom_ZPE"])
    zpe_nf_descriptor = parse_optional_float(descriptor_row["Atom_ZPE_NF"])
    zpe_nf_dump = h_atom["properties"]["ZPE_NF"]
    zpe_nf_primary = zpe_nf_dump if zpe_nf_dump is not None else zpe_nf_descriptor
    zpec_f = parse_optional_float(descriptor_row["Atom_SE_ZPEC"])
    zpec_nf = parse_optional_float(descriptor_row["Atom_SE_NF_ZPEC"])
    bc_f = parse_optional_float(descriptor_row["Atom_BC"])
    bc_nf_raw = parse_optional_float(descriptor_row["Atom_BC_NF"])
    cd = parse_optional_float(descriptor_row["Atom_CD"])
    bc_nf_is_placeholder = key in BC_NF_ZERO_PLACEHOLDER_KEYS
    bc_nf = None if bc_nf_is_placeholder else bc_nf_raw

    flags: list[dict[str, Any]] = []
    if bc_f is None:
        flags.append({"code": "missing_frozen_bader", "severity": "warning"})
    if zpe_nf_primary is None:
        flags.append({"code": "missing_non_frozen_zpe", "severity": "warning"})
    if zpec_nf is None:
        flags.append({"code": "missing_non_frozen_zpe_corrected_energy", "severity": "warning"})
    if bc_nf_is_placeholder:
        flags.append(
            {
                "code": "restored_missing_non_frozen_bader",
                "severity": "warning",
                "detail": "Upstream NA was converted to 0.0 in dump/descriptor; release value restored to null.",
            }
        )
    if key in SC_FECR_NF_ZPE_CONFLICT_KEYS:
        flags.append(
            {
                "code": "dump_descriptor_non_frozen_zpe_conflict",
                "severity": "warning",
                "detail": "Dump ZPE_NF is non-zero; descriptor table ZPE_NF is 0.0 and its reported corrected energy uses 0.0. Both are preserved.",
            }
        )

    calc_f = se_f + zpe_f - HALF_H2_ZPE_EV if None not in (se_f, zpe_f) else None
    calc_nf_primary = se_nf + zpe_nf_primary - HALF_H2_ZPE_EV if None not in (se_nf, zpe_nf_primary) else None
    calc_nf_descriptor = se_nf + zpe_nf_descriptor - HALF_H2_ZPE_EV if None not in (se_nf, zpe_nf_descriptor) else None

    calculations = {
        "frozen": {
            "workflow_id": f"workflow_{host_sample_id}_frozen",
            "constraints": {
                "cell": "fixed",
                "movable_atoms": "selected H only",
                "fixed_atoms": "all non-H atoms",
            },
            "properties": {
                "solution_energy": make_quantity(
                    se_f,
                    "eV",
                    property_status(se_f),
                    descriptor_artifact_id,
                    "Atom_SE",
                    "E_host+H,frozen - E_host - 0.5*E_H2",
                ),
                "zpe_contribution": make_quantity(
                    zpe_f,
                    "eV",
                    property_status(zpe_f),
                    descriptor_artifact_id,
                    "Atom_ZPE",
                    "Three-mode local H zero-point energy in the frozen calculation geometry.",
                ),
                "solution_energy_zpe_corrected": make_quantity(
                    zpec_f,
                    "eV",
                    "reported_derived" if zpec_f is not None else "not_calculated",
                    descriptor_artifact_id,
                    "Atom_SE_ZPEC",
                    "Reported frozen solution energy corrected by the site-H and H2 zero-point-energy terms.",
                    derivation={
                        "expression": "SE + ZPE_site - 0.5*ZPE_H2",
                        "half_zpe_h2_eV": HALF_H2_ZPE_EV,
                        "recomputed_value": calc_f,
                        "difference_from_reported": (calc_f - zpec_f) if None not in (calc_f, zpec_f) else None,
                    },
                ),
                "bader_charge_transfer_as_reported": make_quantity(
                    bc_f,
                    "e",
                    property_status(bc_f),
                    descriptor_artifact_id,
                    "Atom_BC",
                    "Reported H Bader charge-transfer value, preserved without sign or reference-state conversion.",
                ),
            },
        },
        "non_frozen": {
            "workflow_id": f"workflow_{host_sample_id}_non_frozen",
            "constraints": {
                "cell": "fixed",
                "movable_atoms": "all atoms",
                "fixed_atoms": "none",
            },
            "properties": {
                "solution_energy": make_quantity(
                    se_nf,
                    "eV",
                    property_status(se_nf),
                    descriptor_artifact_id,
                    "Atom_SE_NF",
                    "E_host+H,non-frozen - E_host - 0.5*E_H2",
                ),
                "zpe_contribution": make_quantity(
                    zpe_nf_primary,
                    "eV",
                    property_status(zpe_nf_primary),
                    dump_artifact_id if zpe_nf_dump is not None else descriptor_artifact_id,
                    "ZPE_NF" if zpe_nf_dump is not None else "Atom_ZPE_NF",
                    "Three-mode local H zero-point energy evaluated about the non-frozen relaxed geometry.",
                    source_comparison={
                        "dump_value": zpe_nf_dump,
                        "descriptor_table_value": zpe_nf_descriptor,
                        "agreement": (
                            None
                            if zpe_nf_dump is None or zpe_nf_descriptor is None
                            else abs(zpe_nf_dump - zpe_nf_descriptor) <= SCALAR_SOURCE_TOLERANCE
                        ),
                    },
                ),
                "solution_energy_zpe_corrected": make_quantity(
                    zpec_nf,
                    "eV",
                    "reported_derived" if zpec_nf is not None else "not_calculated",
                    descriptor_artifact_id,
                    "Atom_SE_NF_ZPEC",
                    "Reported non-frozen solution energy corrected by the site-H and H2 zero-point-energy terms.",
                    derivation={
                        "expression": "SE_NF + ZPE_NF - 0.5*ZPE_H2",
                        "half_zpe_h2_eV": HALF_H2_ZPE_EV,
                        "recomputed_from_primary_dump_zpe": calc_nf_primary,
                        "difference_using_primary_dump_zpe": (calc_nf_primary - zpec_nf) if None not in (calc_nf_primary, zpec_nf) else None,
                        "recomputed_from_descriptor_table_zpe": calc_nf_descriptor,
                        "difference_using_descriptor_table_zpe": (calc_nf_descriptor - zpec_nf) if None not in (calc_nf_descriptor, zpec_nf) else None,
                    },
                ),
                "bader_charge_transfer_as_reported": make_quantity(
                    bc_nf,
                    "e",
                    "restored_missing" if bc_nf_is_placeholder else property_status(bc_nf),
                    descriptor_artifact_id,
                    "Atom_BC_NF",
                    "Reported H Bader charge-transfer value, preserved without sign or reference-state conversion.",
                    raw_stored_value=bc_nf_raw if bc_nf_is_placeholder else None,
                    normalization_reason="upstream_NA_zero_filled_then_restored_to_null" if bc_nf_is_placeholder else None,
                ),
            },
        },
    }

    geometry_descriptors = {
        "minimum_nearest_neighbor_distance": make_quantity(
            parse_optional_float(descriptor_row["Min_NN_Distance"]),
            "angstrom",
            "reported",
            descriptor_artifact_id,
            "Min_NN_Distance",
            "Shortest distance from the candidate H position to a host atom.",
        ),
        "closest_neighbor_element": {
            "value": descriptor_row["Closest_NN_Element"],
            "unit": None,
            "status": "reported",
            "definition": "Element identity of the closest host atom.",
            "source": {"artifact_id": descriptor_artifact_id, "column": "Closest_NN_Element"},
        },
        "voronoi_volume": make_quantity(
            parse_optional_float(descriptor_row["Voronoi_Volume"]),
            "angstrom^3",
            "reported",
            descriptor_artifact_id,
            "Voronoi_Volume",
            "Volume of the candidate-H Voronoi region in the H-free host.",
        ),
        "number_of_faces": make_quantity(
            int(float(descriptor_row["Num_Faces"])),
            "1",
            "reported",
            descriptor_artifact_id,
            "Num_Faces",
            "Reported number of Voronoi-polyhedron faces.",
        ),
        "surface_area": make_quantity(
            derived_geometry["surface_area"],
            "angstrom^2",
            "recomputed",
            descriptor_artifact_id,
            "Faces+Vertices",
            "Sum of triangular convex-hull face areas.",
        ),
        "volume_cuberoot": make_quantity(
            derived_geometry["volume_cuberoot"],
            "angstrom",
            "recomputed",
            descriptor_artifact_id,
            "Voronoi_Volume",
            "Cube root of the Voronoi volume.",
        ),
        "area_over_volume": make_quantity(
            derived_geometry["area_over_volume"],
            "angstrom^-1",
            "recomputed",
            descriptor_artifact_id,
            "Faces+Vertices+Voronoi_Volume",
            "Convex-hull surface area divided by Voronoi volume.",
        ),
        "sphericity": make_quantity(
            derived_geometry["sphericity"],
            "1",
            "recomputed",
            descriptor_artifact_id,
            "Faces+Vertices+Voronoi_Volume",
            "pi^(1/3)*(6V)^(2/3)/A.",
        ),
        "elongation_lambda1_over_lambda3": make_quantity(
            derived_geometry["elongation_lambda1_over_lambda3"],
            "1",
            "recomputed",
            descriptor_artifact_id,
            "Vertices",
            "Ratio of largest to smallest PCA eigenvalue of the void vertices.",
        ),
        "anisotropy_lambda2_over_lambda3": make_quantity(
            derived_geometry["anisotropy_lambda2_over_lambda3"],
            "1",
            "recomputed",
            descriptor_artifact_id,
            "Vertices",
            "Ratio of middle to smallest PCA eigenvalue of the void vertices.",
        ),
        "rout_over_rin": make_quantity(
            derived_geometry["rout_over_rin"],
            "1",
            "recomputed",
            descriptor_artifact_id,
            "Vertices",
            "Maximum divided by minimum vertex radius relative to the vertex centroid.",
        ),
    }
    if selected_reference_row is not None:
        geometry_descriptors["cross_check"] = {
            "artifact_id": "selected_descriptor_crosscheck",
            "role": "numerical validation only; never authoritative",
            "matched_by": {"Structure": source_structure, "Index": int(descriptor_row["Index"])},
        }

    record = {
        "site_id": site_id,
        "host_sample_id": host_sample_id,
        "initial_sample_id": initial_sample_id,
        "source_key": {
            "structure": source_structure,
            "atom_id": atom_id,
            "descriptor_index": int(descriptor_row["Index"]),
            "site_sequence": site_sequence,
        },
        "candidate_atom": {
            "species": "H",
            "source_type": 2,
            "source_atom_id": atom_id,
            "position": {
                "value": h_atom["position"],
                "unit": "angstrom",
                "coordinate_frame": "source_dump_cartesian",
                "wrapped": False,
                "authoritative_source": dump_artifact_id,
            },
        },
        "construction": {
            "operation": {
                "method": "AddAtom",
                "input_sample": host_sample_id,
                "output_sample": initial_sample_id,
            },
            "host_selection": "retain every source atom with type != 2",
            "added_atom_selection": f"retain only source type-2 atom with id {atom_id}",
            "physical_atom_count": sum(atom["type"] != 2 for atom in dump["atoms"]) + 1,
            "cell_source": dump_artifact_id,
            "periodic_boundary_conditions": [True, True, True],
        },
        "calculations": calculations,
        "pre_insertion_descriptors": {
            "charge_density": make_quantity(
                cd,
                None,
                property_status(cd),
                descriptor_artifact_id,
                "Atom_CD",
                "Total charge integrated over the candidate-H Voronoi region in the H-free host, divided by that region's Voronoi volume.",
                unit_symbol="e/Å³",
                composite_unit={"numerator": "e", "denominator": "angstrom^3"},
            ),
            "geometry": geometry_descriptors,
        },
        "voronoi_topology": {
            "coordinate_unit": "angstrom",
            "vertices": topology["vertices"],
            "edges": topology["edges"],
            "faces": topology["faces"],
            "derivation_details": {
                "vertex_centroid": derived_geometry["vertex_centroid"],
                "pca_eigenvalues_descending": derived_geometry["pca_eigenvalues_descending"],
                "rout": derived_geometry["rout"],
                "rin": derived_geometry["rin"],
            },
        },
        "quality_flags": flags,
        "provenance": {
            "dump_artifact_id": dump_artifact_id,
            "descriptor_artifact_id": descriptor_artifact_id,
            "coordinate_authority": "dump",
            "reported_corrected_energy_authority": "descriptor_table",
            "raw_non_frozen_zpe_authority": "dump",
        },
    }

    validation = {
        "key": {"structure": source_structure, "atom_id": atom_id},
        "frozen_zpec_difference": (calc_f - zpec_f) if None not in (calc_f, zpec_f) else None,
        "non_frozen_zpec_difference_using_primary_dump_zpe": (calc_nf_primary - zpec_nf) if None not in (calc_nf_primary, zpec_nf) else None,
        "non_frozen_zpec_difference_using_descriptor_zpe": (calc_nf_descriptor - zpec_nf) if None not in (calc_nf_descriptor, zpec_nf) else None,
        "quality_flag_codes": [flag["code"] for flag in flags],
    }
    return record, validation


def build_field_dictionary() -> dict[str, Any]:
    return {
        "schema_version": SCHEMA_VERSION,
        "extension_namespace": EXTENSION_NAMESPACE,
        "source_authority_order": [
            "dump: cell, atom ID, LAMMPS type and Cartesian coordinates; primary raw ZPE_NF when sources conflict",
            "Solution_Energy_Descriptors.txt: corrected energies, raw Voronoi topology and tabulated site properties",
            "recomputed geometry: derived paper descriptors from canonical topology",
            "manuscript: crystallographic plane/axis, reported S5m tilt/twist components, family-level GB energy and computational description",
            "derived crystallography: explicit ideal CSL axis-angles and the SCR orientation angle, with formulas and reference directions",
            "author clarification: EDIFF, initial magnetic moments, BC convention, CD definition/unit and missing-value policy",
        ],
        "aggregate_dump_semantics": {
            "type_2_role": "mutually alternative candidate H positions",
            "physical_configuration": "all atoms with type != 2 plus exactly one selected type-2 row",
            "forbidden_interpretation": "all type-2 rows simultaneously occupied by H",
        },
        "crystallographic_angles": {
            "native_json_path": "computational_sample[*].<grain_boundary_block>.misorientation_angle",
            "extended_json_path": "computational_sample[*].x_hydrogen_landscape.angle_metadata",
            "site_collection_json_path": "angle_metadata",
            "manifest_json_path": "structures[*].angle_metadata",
            "status_policy": "S5m tilt/twist are manuscript-reported exactly as printed; lower-angle ideal CSL representatives and the SCR orientation rotation are explicitly marked derived.",
            "families": {
                FAMILY_METADATA[family]["paper_abbreviation"]: build_angle_metadata(family)
                for family in RAW_FAMILIES
            },
        },
        "field_object_shape": "Quantitative JSON paths identify property objects; the numerical datum is stored in their value member.",
        "fields": {
            "Structure": {"json_path": "source_key.structure", "unit": None, "provenance_class": "imported_identifier", "definition": "Raw aggregate-dump structure key, including the _allH suffix."},
            "Index": {"json_path": "source_key.descriptor_index", "unit": None, "provenance_class": "imported_identifier", "definition": "Unique zero-based row index in Solution_Energy_Descriptors.txt."},
            "Atom_id": {"json_path": "source_key.atom_id", "unit": None, "provenance_class": "imported_identifier", "definition": "LAMMPS atom ID of the selected candidate-H row."},
            "Atom_type": {"json_path": "candidate_atom.source_type", "unit": None, "provenance_class": "imported_identifier", "definition": "LAMMPS type of the selected atom; candidate H is always type 2."},
            "Coordinates": {"json_path": "candidate_atom.position.value", "unit": "angstrom", "provenance_class": "imported_authoritative", "definition": "Cartesian x, y, z position copied from the authoritative aggregate dump."},
            "SE": {"json_path": "calculations.frozen.properties.solution_energy", "unit": "eV", "provenance_class": "reported", "source_aliases": ["SE", "Atom_SE"], "definition": "Frozen-host H solution energy E_host+H,F - E_host - 0.5 E_H2."},
            "SE_NF": {"json_path": "calculations.non_frozen.properties.solution_energy", "unit": "eV", "provenance_class": "reported", "source_aliases": ["SE_NF", "Atom_SE_NF"], "definition": "All-ion-relaxed H solution energy at fixed cell: E_host+H,NF - E_host - 0.5 E_H2."},
            "ZPE": {"json_path": "calculations.frozen.properties.zpe_contribution", "unit": "eV", "provenance_class": "reported", "source_aliases": ["ZPE", "Atom_ZPE"], "definition": "Three-mode local-H zero-point energy for the frozen branch."},
            "ZPE_NF": {"json_path": "calculations.non_frozen.properties.zpe_contribution", "unit": "eV", "provenance_class": "reported", "source_aliases": ["ZPE_NF", "Atom_ZPE_NF"], "definition": "Three-mode local-H zero-point energy about the non-frozen relaxed geometry; the dump is primary where sources conflict."},
            "Atom_SE_ZPEC": {"json_path": "calculations.frozen.properties.solution_energy_zpe_corrected", "unit": "eV", "provenance_class": "reported_derived", "release_name": "SE_ZPEC_F", "definition": "Reported frozen ZPE-corrected solution energy: SE + ZPE - 0.5 ZPE_H2."},
            "Atom_SE_NF_ZPEC": {"json_path": "calculations.non_frozen.properties.solution_energy_zpe_corrected", "unit": "eV", "provenance_class": "reported_derived", "release_name": "SE_ZPEC_NF", "definition": "Reported non-frozen ZPE-corrected solution energy: SE_NF + ZPE_NF - 0.5 ZPE_H2."},
            "BC": {"json_path": "calculations.frozen.properties.bader_charge_transfer_as_reported", "unit": "e", "provenance_class": "reported", "definition": "H Bader charge-transfer value preserved exactly as reported, without sign or reference-state conversion."},
            "BC_NF": {"json_path": "calculations.non_frozen.properties.bader_charge_transfer_as_reported", "unit": "e", "provenance_class": "reported_or_restored_missing", "definition": "Non-frozen H Bader charge-transfer value preserved exactly as reported; known upstream NA-to-zero placeholders are restored to null."},
            "CD": {"json_path": "pre_insertion_descriptors.charge_density", "unit": "e/Å³", "unit_ascii": "e/angstrom^3", "provenance_class": "reported_pre_insertion_descriptor", "definition": "Precomputed total charge in the candidate-H Voronoi region of the H-free host divided by the region volume; imported from the descriptor table. Separately distributed CHGCAR resources are listed in charge_density_manifest.json and are not read by this scalar converter."},
            "Voronoi_Volume": {"json_path": "pre_insertion_descriptors.geometry.voronoi_volume", "unit": "angstrom^3", "provenance_class": "reported_pre_insertion_descriptor", "definition": "Volume of the candidate-H Voronoi region in the H-free host."},
            "Num_Faces": {"json_path": "pre_insertion_descriptors.geometry.number_of_faces", "unit": "1", "provenance_class": "reported_pre_insertion_descriptor", "definition": "Number of faces of the candidate-site Voronoi polyhedron."},
            "Min_NN_Distance": {"json_path": "pre_insertion_descriptors.geometry.minimum_nearest_neighbor_distance", "unit": "angstrom", "provenance_class": "reported_pre_insertion_descriptor", "definition": "Shortest distance from the candidate position to a host atom."},
            "Closest_NN_Element": {"json_path": "pre_insertion_descriptors.geometry.closest_neighbor_element", "unit": None, "provenance_class": "reported_pre_insertion_descriptor", "definition": "Chemical element of the closest host atom."},
            "Vertices": {"json_path": "voronoi_topology.vertices", "unit": "angstrom", "provenance_class": "imported_raw_topology", "definition": "Cartesian coordinates of the serialized candidate-site Voronoi vertices."},
            "Edges": {"json_path": "voronoi_topology.edges", "unit": None, "provenance_class": "imported_raw_topology", "definition": "Pairs of zero-based indices into the Vertices array."},
            "Faces": {"json_path": "voronoi_topology.faces", "unit": None, "provenance_class": "imported_raw_topology", "definition": "Triangular faces represented by triples of zero-based indices into the Vertices array."},
        },
        "derived_descriptors": {
            "surface_area": {"json_path": "pre_insertion_descriptors.geometry.surface_area", "unit": "angstrom^2", "provenance_class": "recomputed", "definition": "Sum of the triangular face areas defined by Faces and Vertices."},
            "volume_cuberoot": {"json_path": "pre_insertion_descriptors.geometry.volume_cuberoot", "unit": "angstrom", "provenance_class": "recomputed", "definition": "V^(1/3), where V is the reported Voronoi volume."},
            "area_over_volume": {"json_path": "pre_insertion_descriptors.geometry.area_over_volume", "unit": "angstrom^-1", "provenance_class": "recomputed", "definition": "A/V, where A is the recomputed triangular-face surface area."},
            "sphericity": {"json_path": "pre_insertion_descriptors.geometry.sphericity", "unit": "1", "provenance_class": "recomputed", "definition": "pi^(1/3)*(6V)^(2/3)/A."},
            "elongation_lambda1_over_lambda3": {"json_path": "pre_insertion_descriptors.geometry.elongation_lambda1_over_lambda3", "unit": "1", "provenance_class": "recomputed", "definition": "Largest divided by smallest PCA eigenvalue of the Voronoi vertices."},
            "anisotropy_lambda2_over_lambda3": {"json_path": "pre_insertion_descriptors.geometry.anisotropy_lambda2_over_lambda3", "unit": "1", "provenance_class": "recomputed", "definition": "Middle divided by smallest PCA eigenvalue of the Voronoi vertices."},
            "rout_over_rin": {"json_path": "pre_insertion_descriptors.geometry.rout_over_rin", "unit": "1", "provenance_class": "recomputed", "definition": "Maximum divided by minimum vertex distance from the vertex centroid."},
        },
        "missing_value_policy": {
            "json_representation": "JSON null",
            "no_imputation": True,
            "forbidden_sentinels": ["NaN", "Infinity", "NA", "nan", "zero-filled upstream NA"],
            "status_values": ["reported", "reported_derived", "recomputed", "not_calculated", "restored_missing"],
        },
    }


def build_site_schema() -> dict[str, Any]:
    return {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "$id": "urn:ocdo:hydrogen-solution-landscape:site-collection:1",
        "title": "Hydrogen solution-energy site collection",
        "type": "object",
        "required": ["schema_version", "structure_id", "host_sample_id", "site_count", "angle_metadata", "sites"],
        "properties": {
            "schema_version": {"type": "string"},
            "structure_id": {"type": "string"},
            "host_sample_id": {"type": "string"},
            "site_count": {"type": "integer", "minimum": 1},
            "angle_metadata": {"type": "object"},
            "sites": {
                "type": "array",
                "items": {
                    "type": "object",
                    "required": ["site_id", "source_key", "candidate_atom", "construction", "calculations", "pre_insertion_descriptors", "provenance"],
                    "properties": {
                        "site_id": {"type": "string"},
                        "host_sample_id": {"type": "string"},
                        "initial_sample_id": {"type": "string"},
                        "quality_flags": {"type": "array"},
                    },
                },
            },
        },
    }


def build_manifest_schema() -> dict[str, Any]:
    return {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "$id": "urn:ocdo:hydrogen-solution-landscape:manifest:1",
        "title": "Hydrogen solution-energy metadata manifest",
        "type": "object",
        "required": ["schema_version", "structure_count", "site_count", "structures"],
        "properties": {
            "schema_version": {"type": "string"},
            "structure_count": {"const": 25},
            "site_count": {"const": 1648},
            "structures": {
                "type": "array",
                "minItems": 25,
                "maxItems": 25,
                "items": {
                    "type": "object",
                    "required": ["paper_structure_id", "raw_family", "angle_metadata"],
                    "properties": {"angle_metadata": {"type": "object"}},
                },
            },
        },
    }


def ensure_no_nonfinite(value: Any, path: str = "$") -> None:
    if isinstance(value, float) and not math.isfinite(value):
        raise ValueError(f"Non-finite JSON number at {path}: {value}")
    if isinstance(value, dict):
        for key, child in value.items():
            ensure_no_nonfinite(child, f"{path}.{key}")
    elif isinstance(value, list):
        for index, child in enumerate(value):
            ensure_no_nonfinite(child, f"{path}[{index}]")


def clean_output_directory(output_root: Path, input_root: Path) -> None:
    resolved_output = output_root.resolve()
    resolved_input = input_root.resolve()
    if resolved_output.parent != resolved_input or resolved_output.name != "generated_metadata":
        raise ValueError(f"Refusing to clean unexpected output directory: {resolved_output}")
    if resolved_output.exists():
        shutil.rmtree(resolved_output)
    resolved_output.mkdir(parents=True, exist_ok=True)


def validate_conceptual_dictionary(
    document_paths: Iterable[Path],
    conceptual_source: Path | None,
    document_root: Path,
) -> dict[str, Any]:
    if conceptual_source is None or not conceptual_source.exists():
        return {
            "status": "not_run",
            "reason": "conceptual_dictionary source was not supplied or found",
        }
    paths = list(document_paths)
    sys.path.insert(0, str(conceptual_source.resolve()))
    current_document = None
    try:
        from conceptual_dictionary import ConceptualDict  # type: ignore

        validated_documents = []
        for path in paths:
            current_document = portable_relpath(path, document_root)
            document = ConceptualDict.from_json(path)
            document.validate(strict=True)
            validated_documents.append(current_document)
        return {
            "status": "passed",
            "source_repository": "https://github.com/OCDO/conceptual_dictionary",
            "source_commit": _git_commit(conceptual_source),
            "validation": "ConceptualDict.validate(strict=True)",
            "document_count": len(validated_documents),
            "documents": validated_documents,
        }
    except Exception as exc:  # noqa: BLE001 - validation result must be reported
        return {
            "status": "failed",
            "document": current_document,
            "error": f"{type(exc).__name__}: {exc}",
        }
    finally:
        try:
            sys.path.remove(str(conceptual_source.resolve()))
        except ValueError:
            pass


def _git_commit(path: Path) -> str | None:
    head = path / ".git" / "HEAD"
    if not head.exists():
        return None
    text = head.read_text(encoding="utf-8").strip()
    if text.startswith("ref: "):
        ref = path / ".git" / text[5:]
        return ref.read_text(encoding="utf-8").strip() if ref.exists() else None
    return text


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--descriptor-table", type=Path, default=None)
    parser.add_argument("--selected-descriptors", type=Path, default=None)
    parser.add_argument("--paper", type=Path, default=None)
    parser.add_argument("--incar", type=Path, default=None)
    parser.add_argument("--conceptual-dictionary-source", type=Path, default=None)
    parser.add_argument("--overwrite", action="store_true", help="Allow replacement of an existing generated_metadata directory")
    args = parser.parse_args()

    input_root = args.input_root.resolve()
    output_root = (args.output or (input_root / "generated_metadata")).resolve()
    descriptor_path = (args.descriptor_table or (input_root / DESCRIPTOR_FILENAME)).resolve()
    selected_path = args.selected_descriptors
    if selected_path is None:
        candidate = input_root.parent / SELECTED_DESCRIPTOR_FILENAME
        selected_path = candidate if candidate.exists() else None
    elif selected_path is not None:
        selected_path = selected_path.resolve()
    conceptual_source = args.conceptual_dictionary_source
    if conceptual_source is None:
        candidate = input_root / "tmp" / "conceptual_dictionary_reference"
        conceptual_source = candidate if candidate.exists() else None
    elif conceptual_source is not None:
        conceptual_source = conceptual_source.resolve()

    if not descriptor_path.exists():
        raise FileNotFoundError(descriptor_path)
    if output_root.exists() and not args.overwrite:
        parser.error("Output exists; use a separate input copy or pass --overwrite explicitly")
    generated_at = utc_now()

    descriptor_rows, descriptors_by_key = read_descriptor_table(descriptor_path)
    selected_reference = read_selected_descriptor_reference(selected_path)

    dump_paths = [input_root / f"{family}_{chemistry}_allH" for family in RAW_FAMILIES for chemistry in CHEMISTRIES]
    missing_dump_paths = [path for path in dump_paths if not path.exists()]
    if missing_dump_paths:
        raise FileNotFoundError(f"Missing canonical dumps: {missing_dump_paths}")
    dumps = {path.name: read_dump(path) for path in dump_paths}

    # Parse required inputs before replacing a previous export.
    clean_output_directory(output_root, input_root)

    settings = build_computational_settings()
    write_json(output_root / "computational_settings.json", settings)
    write_json(output_root / "field_dictionary.json", build_field_dictionary())
    write_json(output_root / "schemas" / "hydrogen_site_collection.schema.json", build_site_schema())
    write_json(output_root / "schemas" / "manifest.schema.json", build_manifest_schema())

    source_artifacts: list[dict[str, Any]] = []
    generator_path = Path(__file__).resolve()
    source_artifacts.append(
        artifact_record(
            "generator_script",
            generator_path,
            input_root,
            "deterministic metadata generator used for this release",
            generator_path.parent == input_root,
            "text/x-python",
        )
    )
    readme_path = input_root / "README.md"
    if readme_path.exists():
        source_artifacts.append(
            artifact_record(
                "release_readme",
                readme_path,
                input_root,
                "human-readable release documentation maintained with the generator",
                True,
                "text/markdown",
            )
        )
    source_artifacts.append(
        artifact_record(
            "descriptor_table",
            descriptor_path,
            input_root,
            "canonical per-site descriptor table and raw Voronoi topology",
            descriptor_path.parent == input_root,
            "text/csv",
        )
    )
    if selected_path is not None and selected_path.exists():
        source_artifacts.append(
            artifact_record(
                "selected_descriptor_crosscheck",
                selected_path,
                input_root,
                "external numerical cross-check only; never authoritative",
                False,
                "text/csv",
            )
        )
    if args.paper is not None and args.paper.exists():
        source_artifacts.append(
            artifact_record(
                "unpublished_manuscript",
                args.paper.resolve(),
                input_root,
                "unpublished methodology source; not a formal publication citation",
                False,
                "application/pdf",
            )
        )
    if args.incar is not None and args.incar.exists():
        source_artifacts.append(
            artifact_record(
                "representative_incar",
                args.incar.resolve(),
                input_root,
                "representative input settings; EDIFF corrected by author clarification",
                False,
                "text/plain",
            )
        )

    all_samples: list[dict[str, Any]] = []
    all_workflows: list[dict[str, Any]] = []
    manifest_structures: list[dict[str, Any]] = []
    all_site_ids: set[str] = set()
    all_source_keys: set[tuple[str, int]] = set()
    validation_sites: list[dict[str, Any]] = []
    unexpected_property_mismatches: list[dict[str, Any]] = []
    expected_property_conflicts: list[dict[str, Any]] = []
    coordinate_differences: list[float] = []
    derived_crosscheck_max = defaultdict(float)
    derived_crosscheck_mismatches: list[dict[str, Any]] = []
    out_of_bounds_by_structure: dict[str, list[int]] = {}
    stray_host_properties: list[dict[str, Any]] = []
    site_collection_angle_checks: list[dict[str, Any]] = []
    site_record_angle_duplication_count = 0

    descriptor_keys_used: set[tuple[str, int]] = set()
    selected_keys_used: set[tuple[str, int]] = set()

    for source_structure in sorted(dumps):
        dump = dumps[source_structure]
        family, chemistry = split_structure_name(source_structure)
        paper_structure = paper_structure_name(source_structure)
        paper_family = FAMILY_METADATA[family]["paper_abbreviation"]
        structure_dir = output_root / paper_structure
        structure_dir.mkdir(parents=True, exist_ok=True)
        host_sample, host_extension = build_host_sample(dump, family, chemistry)
        host_sample_id = host_sample["id"]
        workflows = [
            build_workflow(host_sample_id, "frozen", "../computational_settings.json"),
            build_workflow(host_sample_id, "non_frozen", "../computational_settings.json"),
        ]
        all_samples.append(host_sample)
        all_workflows.extend(
            [
                build_workflow(host_sample_id, "frozen", "computational_settings.json"),
                build_workflow(host_sample_id, "non_frozen", "computational_settings.json"),
            ]
        )

        dump_artifact_id = f"dump_{source_structure}"
        dump_artifact = artifact_record(
            dump_artifact_id,
            dump["path"],
            input_root,
            "authoritative aggregate candidate-site map for host cell, atom IDs, types, coordinates and raw site columns",
            True,
        )
        source_artifacts.append(dump_artifact)

        h_atoms = [atom for atom in dump["atoms"] if atom["type"] == 2]
        h_atoms.sort(key=lambda atom: atom["id"])
        site_records: list[dict[str, Any]] = []
        structure_validation: list[dict[str, Any]] = []
        structure_out_of_bounds = []
        for atom in dump["atoms"]:
            if atom["type"] == 2:
                continue
            position = atom["position"]
            if any(
                position[axis] < dump["bounds"][axis][0] or position[axis] > dump["bounds"][axis][1]
                for axis in range(3)
            ):
                structure_out_of_bounds.append(atom["id"])
            nonzero_custom = {
                key: value
                for key, value in atom["properties"].items()
                if value is not None and abs(value) > 0.0
            }
            if nonzero_custom:
                stray_host_properties.append(
                    {"structure": source_structure, "atom_id": atom["id"], "values": nonzero_custom}
                )
        if structure_out_of_bounds:
            out_of_bounds_by_structure[source_structure] = structure_out_of_bounds

        for sequence, h_atom in enumerate(h_atoms, start=1):
            key = (source_structure, h_atom["id"])
            descriptor_row = descriptors_by_key.get(key)
            if descriptor_row is None:
                raise ValueError(f"No descriptor row for {key}")
            descriptor_keys_used.add(key)
            descriptor_position = [float(descriptor_row[axis]) for axis in ("x", "y", "z")]
            max_coord_difference = max(
                abs(h_atom["position"][axis] - descriptor_position[axis]) for axis in range(3)
            )
            coordinate_differences.append(max_coord_difference)
            if max_coord_difference > COORDINATE_TOLERANCE_ANGSTROM:
                raise ValueError(f"Coordinate mismatch for {key}: {max_coord_difference}")

            for dump_column, descriptor_column in DUMP_PROPERTY_TO_DESCRIPTOR.items():
                dump_value = h_atom["properties"][dump_column]
                descriptor_value = parse_optional_float(descriptor_row[descriptor_column])
                if dump_value is None and descriptor_value is None:
                    continue
                if dump_value is None or descriptor_value is None:
                    mismatch = True
                    difference = None
                else:
                    difference = dump_value - descriptor_value
                    mismatch = abs(difference) > SCALAR_SOURCE_TOLERANCE
                if mismatch:
                    entry = {
                        "structure": source_structure,
                        "atom_id": h_atom["id"],
                        "dump_column": dump_column,
                        "descriptor_column": descriptor_column,
                        "dump_value": dump_value,
                        "descriptor_value": descriptor_value,
                        "difference": difference,
                    }
                    if key in SC_FECR_NF_ZPE_CONFLICT_KEYS and dump_column == "ZPE_NF":
                        expected_property_conflicts.append(entry)
                    else:
                        unexpected_property_mismatches.append(entry)

            derived, topology = recompute_geometry(descriptor_row)
            selected_row = selected_reference.get((source_structure, int(descriptor_row["Index"])))
            if selected_row is not None:
                selected_keys_used.add((source_structure, int(descriptor_row["Index"])))
                selected_mapping = {
                    "volume_cuberoot": "V13",
                    "area_over_volume": "A_over_V",
                    "elongation_lambda1_over_lambda3": "Elongation_L1_over_L3",
                    "anisotropy_lambda2_over_lambda3": "Anisotropy_L2_over_L3",
                    "rout_over_rin": "Rout_over_Rin",
                    "sphericity": "Sphericity",
                }
                for local_name, reference_name in selected_mapping.items():
                    reference_value = float(selected_row[reference_name])
                    difference = derived[local_name] - reference_value
                    derived_crosscheck_max[local_name] = max(
                        derived_crosscheck_max[local_name], abs(difference)
                    )
                    if abs(difference) > DERIVED_DESCRIPTOR_TOLERANCE:
                        derived_crosscheck_mismatches.append(
                            {
                                "structure": source_structure,
                                "atom_id": h_atom["id"],
                                "descriptor": local_name,
                                "recomputed": derived[local_name],
                                "reference": reference_value,
                                "difference": difference,
                            }
                        )

            record, site_validation = build_site_record(
                dump,
                h_atom,
                descriptor_row,
                family,
                chemistry,
                sequence,
                derived,
                topology,
                selected_row,
            )
            if record["site_id"] in all_site_ids:
                raise ValueError(f"Duplicate site ID: {record['site_id']}")
            if key in all_source_keys:
                raise ValueError(f"Duplicate source key: {key}")
            all_site_ids.add(record["site_id"])
            all_source_keys.add(key)
            site_records.append(record)
            structure_validation.append(site_validation)
            validation_sites.append(site_validation)

        site_collection = {
            "schema_name": "ocdo_hydrogen_solution_landscape_site_collection",
            "schema_version": SCHEMA_VERSION,
            "extension_namespace": EXTENSION_NAMESPACE,
            "structure_id": paper_structure,
            "source_structure_id": source_structure,
            "host_sample_id": host_sample_id,
            "site_count": len(site_records),
            "angle_metadata": build_angle_metadata(family),
            "record_semantics": "Each record is one alternative one-H calculation; records are not simultaneous occupancies.",
            "reconstruction_rule": "Retain all source atoms with type != 2 and add only the selected type-2 H row.",
            "source_artifacts": {
                "dump": dump_artifact_id,
                "descriptor_table": "descriptor_table",
                "selected_descriptor_crosscheck": "selected_descriptor_crosscheck" if selected_reference else None,
            },
            "sites": site_records,
        }
        site_collection_angle_checks.append(
            {
                "structure": paper_structure,
                "matches_expected_family_metadata": (
                    site_collection["angle_metadata"] == build_angle_metadata(family)
                ),
            }
        )
        site_record_angle_duplication_count += sum(
            "angle_metadata" in site for site in site_records
        )
        ensure_no_nonfinite(site_collection)
        site_filename = f"{paper_structure}.hydrogen_sites.json"
        write_json(structure_dir / site_filename, site_collection)

        conceptual_document = {
            "dataset": {
                "identifier": None,
                "title": f"Hydrogen solution-energy landscape host metadata for {paper_structure}",
                "creators": CREATORS,
                "publication": {"id": None, "identifier": None, "title": None},
                "samples": [host_sample_id],
            },
            "computational_sample": [host_sample],
            "workflow": workflows,
            "operation": [],
            "math_operation": [],
            EXTENSION_KEY: {
                "namespace": EXTENSION_NAMESPACE,
                "schema_version": SCHEMA_VERSION,
                "publication_status": "unpublished",
                "source_structure_id": source_structure,
                "paper_structure_id": paper_structure,
                "site_collection": site_filename,
                "site_count": len(site_records),
                "computational_settings": "../computational_settings.json",
                "field_dictionary": "../field_dictionary.json",
                "angle_metadata": build_angle_metadata(family),
                "raw_aggregate_dump_is_not_a_physical_many_H_sample": True,
                "source_artifact_ids": [dump_artifact_id, "descriptor_table"],
            },
        }
        ensure_no_nonfinite(conceptual_document)
        conceptual_filename = f"{paper_structure}.conceptual.json"
        write_json(structure_dir / conceptual_filename, conceptual_document)

        summaries = {
            "solution_energy_frozen_eV": quantity_summary(
                site["calculations"]["frozen"]["properties"]["solution_energy"]["value"] for site in site_records
            ),
            "solution_energy_non_frozen_eV": quantity_summary(
                site["calculations"]["non_frozen"]["properties"]["solution_energy"]["value"] for site in site_records
            ),
            "solution_energy_zpe_corrected_frozen_eV": quantity_summary(
                site["calculations"]["frozen"]["properties"]["solution_energy_zpe_corrected"]["value"] for site in site_records
            ),
            "solution_energy_zpe_corrected_non_frozen_eV": quantity_summary(
                site["calculations"]["non_frozen"]["properties"]["solution_energy_zpe_corrected"]["value"] for site in site_records
            ),
        }
        manifest_structures.append(
            {
                "source_structure_id": source_structure,
                "paper_structure_id": paper_structure,
                "host_sample_id": host_sample_id,
                "raw_family": family,
                "paper_family": paper_family,
                "chemistry": chemistry,
                "structure_class": FAMILY_METADATA[family]["structure_class"],
                "host_atom_count": sum(atom["type"] != 2 for atom in dump["atoms"]),
                "hydrogen_site_count": len(site_records),
                "conceptual_file": f"{paper_structure}/{conceptual_filename}",
                "site_file": f"{paper_structure}/{site_filename}",
                "source_dump_artifact_id": dump_artifact_id,
                "angle_metadata": build_angle_metadata(family),
                "summary": summaries,
            }
        )

    unused_descriptor_keys = set(descriptors_by_key) - descriptor_keys_used
    selected_expected_keys = {
        (row["Structure"], int(row["Index"])) for row in descriptor_rows
    }
    unused_selected_keys = set(selected_reference) - selected_keys_used
    missing_selected_keys = selected_expected_keys - selected_keys_used if selected_reference else set()

    if len(all_site_ids) != 1648:
        raise ValueError(f"Expected 1,648 sites, generated {len(all_site_ids)}")
    if len(manifest_structures) != 25:
        raise ValueError(f"Expected 25 structures, generated {len(manifest_structures)}")

    aggregate_core = {
        "dataset": {
            "identifier": None,
            "title": "Hydrogen solution-energy landscapes in ferritic Fe environments",
            "creators": CREATORS,
            "publication": {"id": None, "identifier": None, "title": None},
            "samples": [sample["id"] for sample in all_samples],
        },
        "computational_sample": all_samples,
        "workflow": all_workflows,
        "operation": [],
        "math_operation": [],
        EXTENSION_KEY: {
            "namespace": EXTENSION_NAMESPACE,
            "schema_version": SCHEMA_VERSION,
            "publication_status": "unpublished",
            "citation": None,
            "doi": None,
            "host_structure_count": 25,
            "one_H_site_record_count": 1648,
            "manifest": "manifest_index.json",
            "site_records": "per-structure *.hydrogen_sites.json files",
            "raw_aggregate_dump_is_not_a_physical_many_H_sample": True,
        },
    }
    ensure_no_nonfinite(aggregate_core)
    aggregate_core_path = output_root / "dataset.conceptual.json"
    write_json(aggregate_core_path, aggregate_core)

    conceptual_validation = validate_conceptual_dictionary(
        sorted(output_root.rglob("*.conceptual.json")),
        conceptual_source,
        output_root,
    )

    expected_frozen_missing = sorted(
        [
            entry["key"]
            for entry in validation_sites
            if "missing_frozen_bader" in entry["quality_flag_codes"]
        ],
        key=lambda item: (item["structure"], item["atom_id"]),
    )
    expected_nf_zpe_missing = sorted(
        [
            entry["key"]
            for entry in validation_sites
            if "missing_non_frozen_zpe" in entry["quality_flag_codes"]
        ],
        key=lambda item: (item["structure"], item["atom_id"]),
    )
    restored_bc_nf = sorted(
        [
            entry["key"]
            for entry in validation_sites
            if "restored_missing_non_frozen_bader" in entry["quality_flag_codes"]
        ],
        key=lambda item: (item["structure"], item["atom_id"]),
    )
    primary_nf_formula_mismatches = [
        entry
        for entry in validation_sites
        if entry["non_frozen_zpec_difference_using_primary_dump_zpe"] is not None
        and abs(entry["non_frozen_zpec_difference_using_primary_dump_zpe"]) > ENERGY_TOLERANCE_EV
    ]
    descriptor_nf_formula_mismatches = [
        entry
        for entry in validation_sites
        if entry["non_frozen_zpec_difference_using_descriptor_zpe"] is not None
        and abs(entry["non_frozen_zpec_difference_using_descriptor_zpe"]) > ENERGY_TOLERANCE_EV
    ]
    frozen_formula_mismatches = [
        entry
        for entry in validation_sites
        if entry["frozen_zpec_difference"] is not None
        and abs(entry["frozen_zpec_difference"]) > ENERGY_TOLERANCE_EV
    ]
    summary_count_mismatches = [
        {
            "structure": structure["paper_structure_id"],
            "quantity": quantity,
            "count": summary["count"],
            "missing": summary["missing"],
            "expected_total": structure["hydrogen_site_count"],
        }
        for structure in manifest_structures
        for quantity, summary in structure["summary"].items()
        if summary["count"] + summary["missing"] != structure["hydrogen_site_count"]
    ]
    angle_metadata_mismatches: list[dict[str, Any]] = []
    for sample in all_samples:
        extension = sample[EXTENSION_KEY]
        family = extension["raw_family"]
        expected = build_angle_metadata(family)
        if extension.get("angle_metadata") != expected:
            angle_metadata_mismatches.append(
                {"location": "computational_sample", "sample_id": sample["id"], "family": family}
            )
        gb = FAMILY_METADATA[family].get("grain_boundary")
        if gb is not None:
            native_angle = sample[gb["native_key"]].get("misorientation_angle")
            if native_angle != gb["misorientation_angle_degrees"]:
                angle_metadata_mismatches.append(
                    {"location": "native_grain_boundary_block", "sample_id": sample["id"], "family": family}
                )
    for structure in manifest_structures:
        family = structure["raw_family"]
        if structure.get("angle_metadata") != build_angle_metadata(family):
            angle_metadata_mismatches.append(
                {"location": "manifest", "structure": structure["paper_structure_id"], "family": family}
            )
    for check in site_collection_angle_checks:
        if not check["matches_expected_family_metadata"]:
            angle_metadata_mismatches.append(
                {"location": "site_collection", "structure": check["structure"]}
            )
    if site_record_angle_duplication_count:
        angle_metadata_mismatches.append(
            {
                "location": "individual_site_records",
                "reason": "family-level angle metadata must not be duplicated inside sites",
                "count": site_record_angle_duplication_count,
            }
        )
    native_gb_misorientation_count = sum(
        FAMILY_METADATA[sample[EXTENSION_KEY]["raw_family"]].get("grain_boundary") is not None
        and sample[
            FAMILY_METADATA[sample[EXTENSION_KEY]["raw_family"]]["grain_boundary"]["native_key"]
        ].get("misorientation_angle") is not None
        for sample in all_samples
    )
    scr_orientation_rotation_count = sum(
        sample[EXTENSION_KEY]["raw_family"] == "SC"
        and sample[EXTENSION_KEY]["angle_metadata"]["orientation_rotation"]["angle_magnitude"]["value"] is not None
        for sample in all_samples
    )
    s5m_mixed_component_count = sum(
        sample[EXTENSION_KEY]["raw_family"] == "S5m"
        and "mixed_character_components" in sample[EXTENSION_KEY]["angle_metadata"]
        for sample in all_samples
    )
    expected_angle_coverage = {
        "host_structures": 25,
        "site_collections": 25,
        "manifest_structures": 25,
        "native_gb_misorientations": 20,
        "scr_orientation_rotations": 5,
        "s5m_mixed_component_sets": 5,
    }
    actual_angle_coverage = {
        "host_structures": len(all_samples),
        "site_collections": len(site_collection_angle_checks),
        "manifest_structures": len(manifest_structures),
        "native_gb_misorientations": native_gb_misorientation_count,
        "scr_orientation_rotations": scr_orientation_rotation_count,
        "s5m_mixed_component_sets": s5m_mixed_component_count,
    }
    if actual_angle_coverage != expected_angle_coverage:
        angle_metadata_mismatches.append(
            {
                "location": "coverage",
                "expected": expected_angle_coverage,
                "actual": actual_angle_coverage,
            }
        )

    errors: list[dict[str, Any]] = []
    warnings: list[dict[str, Any]] = []
    if unexpected_property_mismatches:
        errors.append({"code": "unexpected_dump_descriptor_property_mismatch", "entries": unexpected_property_mismatches})
    if unused_descriptor_keys:
        errors.append({"code": "unused_descriptor_rows", "keys": sorted(unused_descriptor_keys)})
    if derived_crosscheck_mismatches:
        errors.append({"code": "derived_descriptor_crosscheck_mismatch", "entries": derived_crosscheck_mismatches})
    if frozen_formula_mismatches:
        errors.append({"code": "frozen_zpec_formula_mismatch", "entries": frozen_formula_mismatches})
    if descriptor_nf_formula_mismatches:
        errors.append({"code": "descriptor_internal_non_frozen_zpec_formula_mismatch", "entries": descriptor_nf_formula_mismatches})
    if summary_count_mismatches:
        errors.append({"code": "manifest_summary_count_mismatch", "entries": summary_count_mismatches})
    if angle_metadata_mismatches:
        errors.append({"code": "crystallographic_angle_metadata_mismatch", "entries": angle_metadata_mismatches})
    if conceptual_validation["status"] == "failed":
        errors.append({"code": "conceptual_dictionary_validation_failed", "detail": conceptual_validation})

    if expected_property_conflicts:
        warnings.append({"code": "known_SC_FeCr_non_frozen_ZPE_conflict", "entries": expected_property_conflicts})
    if primary_nf_formula_mismatches:
        warnings.append({"code": "reported_NF_ZPEC_inconsistent_with_primary_dump_ZPE", "entries": primary_nf_formula_mismatches})
    if out_of_bounds_by_structure:
        warnings.append({"code": "unwrapped_host_coordinates_outside_nominal_bounds", "entries": out_of_bounds_by_structure})
    if stray_host_properties:
        warnings.append({"code": "nonzero_custom_property_on_non_H_row_ignored", "entries": stray_host_properties})
    if conceptual_validation["status"] == "not_run":
        warnings.append({"code": "conceptual_dictionary_validation_not_run", "detail": conceptual_validation})
    if unused_selected_keys:
        warnings.append({"code": "unused_selected_descriptor_crosscheck_rows", "count": len(unused_selected_keys)})
    if missing_selected_keys:
        warnings.append({"code": "missing_selected_descriptor_crosscheck_rows", "count": len(missing_selected_keys)})

    validation_report = {
        "schema_version": SCHEMA_VERSION,
        "generated_at_utc": generated_at,
        "status": "passed_with_warnings" if not errors and warnings else ("passed" if not errors else "failed"),
        "summary": {
            "errors": len(errors),
            "warnings": len(warnings),
            "structure_count": len(manifest_structures),
            "site_count": len(all_site_ids),
            "unique_source_key_count": len(all_source_keys),
            "descriptor_row_count": len(descriptor_rows),
            "selected_crosscheck_row_count": len(selected_reference),
        },
        "checks": {
            "all_25_canonical_dumps_parsed": len(dumps) == 25,
            "all_1648_site_keys_unique": len(all_source_keys) == 1648,
            "all_descriptor_rows_joined": not unused_descriptor_keys and len(descriptor_keys_used) == 1648,
            "maximum_dump_descriptor_coordinate_difference_angstrom": max(coordinate_differences),
            "coordinate_tolerance_angstrom": COORDINATE_TOLERANCE_ANGSTROM,
            "unexpected_property_mismatch_count": len(unexpected_property_mismatches),
            "known_SC_FeCr_ZPE_NF_conflict_count": len(expected_property_conflicts),
            "derived_descriptor_crosscheck_max_abs_difference": dict(derived_crosscheck_max),
            "derived_descriptor_crosscheck_tolerance": DERIVED_DESCRIPTOR_TOLERANCE,
            "frozen_zpec_formula_mismatch_count": len(frozen_formula_mismatches),
            "descriptor_internal_non_frozen_zpec_formula_mismatch_count": len(descriptor_nf_formula_mismatches),
            "primary_dump_non_frozen_zpec_formula_mismatch_count": len(primary_nf_formula_mismatches),
            "manifest_summary_count_mismatch_count": len(summary_count_mismatches),
            "crystallographic_angles": {
                "status": "passed" if not angle_metadata_mismatches else "failed",
                "host_structure_count": len(all_samples),
                "site_collection_count": len(site_collection_angle_checks),
                "manifest_structure_count": len(manifest_structures),
                "family_count": len(RAW_FAMILIES),
                "native_gb_misorientation_count": native_gb_misorientation_count,
                "scr_orientation_rotation_count": scr_orientation_rotation_count,
                "s5m_mixed_component_set_count": s5m_mixed_component_count,
                "site_record_duplication_count": site_record_angle_duplication_count,
                "mismatch_count": len(angle_metadata_mismatches),
                "families": {
                    FAMILY_METADATA[family]["paper_abbreviation"]: build_angle_metadata(family)
                    for family in RAW_FAMILIES
                },
            },
            "conceptual_dictionary": conceptual_validation,
        },
        "missing_values": {
            "frozen_bader": {"count": len(expected_frozen_missing), "keys": expected_frozen_missing},
            "non_frozen_zpe_and_zpec": {"count": len(expected_nf_zpe_missing), "keys": expected_nf_zpe_missing},
            "non_frozen_bader_restored_from_zero_placeholder": {"count": len(restored_bc_nf), "keys": restored_bc_nf},
        },
        "errors": errors,
        "warnings": warnings,
    }
    ensure_no_nonfinite(validation_report)
    validation_path = output_root / "validation_report.json"
    write_json(validation_path, validation_report)

    source_manifest = {
        "schema_version": SCHEMA_VERSION,
        "generated_at_utc": generated_at,
        "canonical_source_root": "..",
        "artifact_path_base": "canonical_source_root",
        "external_artifact_path_policy": "path is null; filename, role, size and checksum are retained without a machine-local location",
        "publication_status": "unpublished",
        "citation": None,
        "doi": None,
        "authority_order": build_field_dictionary()["source_authority_order"],
        "artifacts": sorted(source_artifacts, key=lambda item: item["artifact_id"]),
        "not_in_canonical_package": {
            "relaxed_H_containing_geometries": True,
            "raw_total_energies_for_every_branch": True,
            "vibrational_hessians": True,
        },
    }
    resource_catalog = input_root / "charge_density_manifest.json"
    if resource_catalog.exists():
        source_manifest["charge_density_resources"] = {
            "catalog_path": "charge_density_manifest.json",
            "path_base": "canonical_source_root",
            "storage": "separate downloadable release assets; not scalar-converter inputs",
            "catalog_sha256": sha256_file(resource_catalog),
        }
    write_json(output_root / "source_manifest.json", source_manifest)

    (output_root / "README.md").write_text(
        "# Generated Metadata\n\n"
        "Start with [the repository README](../README.md) and "
        "[the scientific data guide](../DATA_GUIDE.md).\n\n"
        "Use manifest_index.json to locate host and H-site records. "
        "Each site represents an alternative single-H calculation, not "
        "simultaneous H occupancy.\n",
        encoding="utf-8",
    )

    output_files = sorted(
        path for path in output_root.rglob("*") if path.is_file() and path.name != "manifest_index.json"
    )
    output_artifacts = [
        {
            "path": portable_relpath(path, output_root),
            "byte_size": path.stat().st_size,
            "sha256": sha256_file(path),
        }
        for path in output_files
    ]
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "generator": {"name": GENERATOR_NAME, "version": GENERATOR_VERSION},
        "generated_at_utc": generated_at,
        "dataset_title": "Hydrogen solution-energy landscapes in ferritic Fe environments",
        "publication_status": "unpublished",
        "citation": None,
        "doi": None,
        "structure_count": len(manifest_structures),
        "site_count": len(all_site_ids),
        "conceptual_dictionary_core": "dataset.conceptual.json",
        "computational_settings": "computational_settings.json",
        "field_dictionary": "field_dictionary.json",
        "source_manifest": "source_manifest.json",
        "validation_report": "validation_report.json",
        "structures": sorted(manifest_structures, key=lambda item: item["paper_structure_id"]),
        "output_artifacts": output_artifacts,
    }
    write_json(output_root / "manifest_index.json", manifest)

    if errors:
        print(json.dumps({"status": "failed", "errors": len(errors), "validation_report": str(validation_path)}, indent=2))
        return 1
    print(
        json.dumps(
            {
                "status": validation_report["status"],
                "structures": len(manifest_structures),
                "sites": len(all_site_ids),
                "warnings": len(warnings),
                "output": str(output_root),
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
