#!/usr/bin/env python3
"""Calculate all pairwise RdRp structural alignments with US-align.

The script writes a tidy table containing both directional TM-scores, a
symmetric average-length-normalized TM-score, the aligned-residue RMSD, and
alignment coverage.  It also writes square RMSD/TM-score matrices and a
dependency-free SVG heatmap for inclusion in the repository README.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import html
import itertools
import json
import math
import re
import statistics
import subprocess
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path


STRUCTURES = [
    ("AhPV1", "ahpv1_rdrp_relax_m3_p0_plddt-90.pdb"),
    ("ALV1", "alv1_p1_1140_1610.pdb"),
    ("CaMV", "CaMV_P5_8R0S.pdb"),
    ("CMV", "CMV_2a_273-750.pdb"),
    ("HIV-1 RT", "HIV1RT_3DLK.pdb"),
    ("R2 RT", "nonLTR_RT_8gh6_1_924.pdb"),
    ("PlAMV", "PlAMV_RdRp_q07518_895_1385.pdb"),
    ("PLrV", "PLrV_Polerovirus_P11623_relax_m4_p0_plddt-78.pdb"),
    ("RYMV", "rymv_rdrp_1_464.pdb"),
    ("TCV", "tcv_rdrp_relax_m1_p0_plddt-92.pdb"),
    ("TMV", "tmv_rdrp_1117_relax_m3_p0_plddt-91.pdb"),
    ("TRoV", "trov_p2ab_428.pdb"),
    ("TRV", "trv_134k_1206_1707.pdb"),
    ("TuMV", "TuMV_NIb_m2_plddt-93.pdb"),
    ("TuYV", "TuYV_RdRP_p09507_relax_m1_p0_plddt-80.pdb"),
    ("TYMV", "tymv_206k_1298.pdb"),
    ("YoMV", "YoMV_RdRP_q66220_1120_1597.pdb"),
]


@dataclass
class Cluster:
    """A node in an average-linkage UPGMA dendrogram."""

    name: str
    members: tuple[int, ...]
    height: float
    left: "Cluster | None" = None
    right: "Cluster | None" = None
    merge_distance: float = 0.0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--usalign", required=True, type=Path,
                        help="Path to the official USalign executable")
    parser.add_argument("--input-dir", default=Path("RdRps"), type=Path)
    parser.add_argument(
        "--output-dir",
        default=Path("Results/pairwise_structure_validation"),
        type=Path,
    )
    parser.add_argument(
        "--figure",
        default=Path("Figures/pairwise_RMSD_TMscore_heatmaps.svg"),
        type=Path,
    )
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def first_protein_chain(path: Path) -> str:
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("ATOM  "):
                return line[21].strip() or "_"
    raise ValueError(f"No ATOM records found in {path}")


def extract(pattern: str, text: str, name: str, cast=float):
    match = re.search(pattern, text, flags=re.MULTILINE)
    if not match:
        raise ValueError(f"Could not parse {name} from US-align output")
    return cast(match.group(1))


def align_pair(executable: Path, path_1: Path, path_2: Path) -> dict[str, float | int]:
    command = [
        str(executable), str(path_1), str(path_2),
        "-mol", "prot", "-a", "T", "-outfmt", "-1",
    ]
    completed = subprocess.run(
        command, check=True, capture_output=True, text=True, errors="replace"
    )
    output = completed.stdout
    return {
        "length_1": extract(r"Length of Structure_1:\s*(\d+)", output, "length_1", int),
        "length_2": extract(r"Length of Structure_2:\s*(\d+)", output, "length_2", int),
        "aligned_length": extract(r"Aligned length=\s*(\d+)", output, "aligned_length", int),
        "rmsd_angstrom": extract(r"Aligned length=.*?RMSD=\s*([0-9.]+)", output, "RMSD"),
        "sequence_identity": extract(r"Seq_ID=n_identical/n_aligned=\s*([0-9.]+)", output, "identity"),
        "tm_score_norm_structure_1": extract(
            r"TM-score=\s*([0-9.]+)\s*\(normalized by length of Structure_1", output, "TM1"
        ),
        "tm_score_norm_structure_2": extract(
            r"TM-score=\s*([0-9.]+)\s*\(normalized by length of Structure_2", output, "TM2"
        ),
        "tm_score_norm_average_length": extract(
            r"TM-score=\s*([0-9.]+)\s*\(if normalized by average length", output, "TMavg"
        ),
    }


def write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_matrix(path: Path, labels: list[str], matrix: list[list[float]], digits: int) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["structure", *labels])
        for label, row in zip(labels, matrix):
            writer.writerow([label, *[f"{value:.{digits}f}" for value in row]])


def upgma(labels: list[str], matrix: list[list[float]]) -> tuple[Cluster, list[dict]]:
    """Perform deterministic unweighted pair-group clustering by arithmetic mean."""
    active = [Cluster(label, (index,), 0.0) for index, label in enumerate(labels)]
    merges: list[dict] = []
    step = 0
    while len(active) > 1:
        candidates = []
        for first, second in itertools.combinations(active, 2):
            distances = [matrix[i][j] for i in first.members for j in second.members]
            distance = sum(distances) / len(distances)
            tie_breaker = tuple(sorted((first.name, second.name)))
            candidates.append((distance, tie_breaker, first, second))
        distance, _, left, right = min(candidates, key=lambda item: (item[0], item[1]))
        if right.name < left.name:
            left, right = right, left
        step += 1
        merged = Cluster(
            name=f"node_{step}",
            members=tuple(sorted(left.members + right.members)),
            height=distance / 2.0,
            left=left,
            right=right,
            merge_distance=distance,
        )
        merges.append({
            "step": step,
            "left_cluster": ";".join(labels[index] for index in left.members),
            "right_cluster": ";".join(labels[index] for index in right.members),
            "cluster_size": len(merged.members),
            "average_linkage_distance": distance,
            "node_height": merged.height,
        })
        active = [cluster for cluster in active if cluster is not left and cluster is not right]
        active.append(merged)
    return active[0], merges


def newick_label(label: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", label).strip("_")


def to_newick(node: Cluster) -> str:
    def render(current: Cluster, parent_height: float | None) -> str:
        branch_length = 0.0 if parent_height is None else max(0.0, parent_height - current.height)
        suffix = "" if parent_height is None else f":{branch_length:.6f}"
        if current.left is None or current.right is None:
            return f"{newick_label(current.name)}{suffix}"
        children = f"{render(current.left, current.height)},{render(current.right, current.height)}"
        return f"({children}){suffix}"
    return render(node, None) + ";\n"


def interpolate_color(value: float, minimum: float, maximum: float,
                      low: tuple[int, int, int], high: tuple[int, int, int]) -> str:
    fraction = 0.5 if maximum == minimum else (value - minimum) / (maximum - minimum)
    fraction = max(0.0, min(1.0, fraction))
    rgb = tuple(round(a + fraction * (b - a)) for a, b in zip(low, high))
    return "#{:02x}{:02x}{:02x}".format(*rgb)


def luminance(hex_color: str) -> float:
    red, green, blue = (int(hex_color[i:i + 2], 16) for i in (1, 3, 5))
    return 0.2126 * red + 0.7152 * green + 0.0722 * blue


def heatmap_svg(labels: list[str], rmsd: list[list[float]], tm: list[list[float]]) -> str:
    cell = 42
    left = 112
    top = 175
    panel = cell * len(labels)
    gap = 120
    width = left * 2 + panel * 2 + gap
    height = top + panel + 110
    rmsd_max = max(max(row) for row in rmsd)
    pieces = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<style>text{font-family:Arial,sans-serif}.title{font-size:24px;font-weight:700}.axis{font-size:14px}.cell{font-size:11px;font-weight:600}</style>',
        f'<text x="{width / 2}" y="34" text-anchor="middle" class="title">All-pairwise RdRp structural validation (US-align)</text>',
    ]

    def panel_svg(x0: int, matrix: list[list[float]], title: str,
                  minimum: float, maximum: float, low, high, digits: int) -> None:
        pieces.append(f'<text x="{x0 + panel / 2}" y="73" text-anchor="middle" class="title">{html.escape(title)}</text>')
        for index, label in enumerate(labels):
            x = x0 + index * cell + cell / 2
            y = top - 12
            pieces.append(f'<text x="{x}" y="{y}" transform="rotate(-55 {x} {y})" text-anchor="start" class="axis">{html.escape(label)}</text>')
            pieces.append(f'<text x="{x0 - 9}" y="{top + index * cell + cell * 0.68}" text-anchor="end" class="axis">{html.escape(label)}</text>')
        for row_index, row in enumerate(matrix):
            for column_index, value in enumerate(row):
                fill = interpolate_color(value, minimum, maximum, low, high)
                text_color = "white" if luminance(fill) < 125 else "#111111"
                x = x0 + column_index * cell
                y = top + row_index * cell
                pieces.append(f'<rect x="{x}" y="{y}" width="{cell}" height="{cell}" fill="{fill}" stroke="white" stroke-width="1"/>')
                pieces.append(f'<text x="{x + cell / 2}" y="{y + cell * 0.64}" text-anchor="middle" fill="{text_color}" class="cell">{value:.{digits}f}</text>')
        pieces.append(f'<text x="{x0 + panel / 2}" y="{top + panel + 40}" text-anchor="middle" class="axis">Structure 2</text>')
        pieces.append(f'<text x="{x0 - 84}" y="{top + panel / 2}" transform="rotate(-90 {x0 - 84} {top + panel / 2})" text-anchor="middle" class="axis">Structure 1</text>')

    panel_svg(left, rmsd, "Aligned-residue RMSD (Å)", 0.0, rmsd_max,
              (238, 247, 255), (8, 48, 107), 1)
    panel_svg(left + panel + gap, tm, "TM-score (average-length normalized)", 0.0, 1.0,
              (255, 247, 236), (127, 39, 4), 2)
    pieces.append(f'<text x="{width / 2}" y="{height - 22}" text-anchor="middle" class="axis">Diagonal values represent self-comparisons; all off-diagonal cells correspond to the 136 unique pairwise alignments.</text>')
    pieces.append("</svg>")
    return "\n".join(pieces)


def dendrogram_svg(rmsd_tree: Cluster, tm_tree: Cluster) -> str:
    """Create paired horizontal UPGMA dendrograms without plotting dependencies."""
    width = 1700
    height = 760
    panel_width = 650
    tree_width = 470
    top = 110
    row_spacing = 30
    pieces = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<style>text{font-family:Arial,sans-serif}.title{font-size:24px;font-weight:700}.subtitle{font-size:17px;font-weight:700}.label{font-size:14px}.value{font-size:11px;fill:#7f0000}.branch{stroke:#222;stroke-width:2;fill:none}</style>',
        f'<text x="{width / 2}" y="35" text-anchor="middle" class="title">UPGMA clustering of pairwise structural-alignment metrics</text>',
    ]

    def leaves(node: Cluster) -> list[Cluster]:
        if node.left is None or node.right is None:
            return [node]
        return leaves(node.left) + leaves(node.right)

    def draw_panel(node: Cluster, x0: int, title: str, unit: str, digits: int) -> None:
        leaf_nodes = leaves(node)
        y_map = {id(leaf): top + index * row_spacing for index, leaf in enumerate(leaf_nodes)}
        root_distance = node.merge_distance or 1.0
        pieces.append(f'<text x="{x0 + panel_width / 2}" y="73" text-anchor="middle" class="subtitle">{html.escape(title)}</text>')

        def coordinates(current: Cluster) -> tuple[float, float]:
            linkage_distance = 2.0 * current.height
            x = x0 + tree_width * (1.0 - linkage_distance / root_distance)
            if current.left is None or current.right is None:
                y = y_map[id(current)]
                pieces.append(f'<text x="{x0 + tree_width + 12}" y="{y + 5}" class="label">{html.escape(current.name)}</text>')
                return x, y
            left_x, left_y = coordinates(current.left)
            right_x, right_y = coordinates(current.right)
            y = (left_y + right_y) / 2.0
            pieces.append(f'<line x1="{x}" y1="{left_y}" x2="{x}" y2="{right_y}" class="branch"/>')
            pieces.append(f'<line x1="{x}" y1="{left_y}" x2="{left_x}" y2="{left_y}" class="branch"/>')
            pieces.append(f'<line x1="{x}" y1="{right_y}" x2="{right_x}" y2="{right_y}" class="branch"/>')
            pieces.append(f'<text x="{x + 4}" y="{y - 4}" class="value">{current.merge_distance:.{digits}f}</text>')
            return x, y

        draw_x, _ = coordinates(node)
        scale_y = top + len(leaf_nodes) * row_spacing + 25
        pieces.append(f'<line x1="{draw_x}" y1="{scale_y}" x2="{x0 + tree_width}" y2="{scale_y}" class="branch"/>')
        pieces.append(f'<text x="{draw_x}" y="{scale_y + 20}" text-anchor="start" class="label">{root_distance:.{digits}f}</text>')
        pieces.append(f'<text x="{x0 + tree_width}" y="{scale_y + 20}" text-anchor="end" class="label">0</text>')
        pieces.append(f'<text x="{x0 + tree_width / 2}" y="{scale_y + 46}" text-anchor="middle" class="label">Average-linkage distance ({html.escape(unit)})</text>')

    draw_panel(rmsd_tree, 90, "RMSD distance", "Å", 2)
    draw_panel(tm_tree, 925, "TM distance = 1 − TM-score", "unitless", 3)
    pieces.append('<text x="850" y="742" text-anchor="middle" class="label">Numbers at internal nodes are the UPGMA average-linkage distances.</text>')
    pieces.append("</svg>")
    return "\n".join(pieces)


def main() -> None:
    args = parse_args()
    executable = args.usalign.resolve()
    input_dir = args.input_dir.resolve()
    output_dir = args.output_dir.resolve()
    figure = args.figure.resolve()
    if not executable.is_file():
        raise FileNotFoundError(executable)
    missing = [filename for _, filename in STRUCTURES if not (input_dir / filename).is_file()]
    if missing:
        raise FileNotFoundError(f"Missing input structures: {', '.join(missing)}")
    output_dir.mkdir(parents=True, exist_ok=True)
    figure.parent.mkdir(parents=True, exist_ok=True)

    # US-align prints its version correctly but some Windows builds return 1
    # for the informational ``-v`` invocation, so do not treat that code as an
    # alignment failure.
    version_output = subprocess.run(
        [str(executable), "-v"], check=False, capture_output=True, text=True, errors="replace"
    ).stdout
    version_match = re.search(r"Version\s+([0-9]+)", version_output)
    version = version_match.group(1) if version_match else version_output.strip()

    results = []
    total_pairs = math.comb(len(STRUCTURES), 2)
    for pair_index, ((label_1, filename_1), (label_2, filename_2)) in enumerate(
        itertools.combinations(STRUCTURES, 2), start=1
    ):
        print(f"[{pair_index:3d}/{total_pairs}] {label_1} vs {label_2}", flush=True)
        path_1 = input_dir / filename_1
        path_2 = input_dir / filename_2
        metrics = align_pair(executable, path_1, path_2)
        results.append({
            "structure_1": label_1,
            "file_1": filename_1,
            "chain_1": first_protein_chain(path_1),
            "structure_2": label_2,
            "file_2": filename_2,
            "chain_2": first_protein_chain(path_2),
            **metrics,
            "coverage_structure_1": metrics["aligned_length"] / metrics["length_1"],
            "coverage_structure_2": metrics["aligned_length"] / metrics["length_2"],
        })

    fieldnames = [
        "structure_1", "file_1", "chain_1", "length_1",
        "structure_2", "file_2", "chain_2", "length_2",
        "aligned_length", "rmsd_angstrom", "sequence_identity",
        "tm_score_norm_structure_1", "tm_score_norm_structure_2",
        "tm_score_norm_average_length", "coverage_structure_1",
        "coverage_structure_2",
    ]
    formatted_results = []
    for result in results:
        formatted = dict(result)
        for key in (
            "rmsd_angstrom", "sequence_identity", "tm_score_norm_structure_1",
            "tm_score_norm_structure_2", "tm_score_norm_average_length",
            "coverage_structure_1", "coverage_structure_2",
        ):
            formatted[key] = f"{result[key]:.5f}"
        formatted_results.append(formatted)
    write_csv(output_dir / "pairwise_alignments.csv", formatted_results, fieldnames)

    labels = [label for label, _ in STRUCTURES]
    index = {label: position for position, label in enumerate(labels)}
    rmsd_matrix = [[0.0 for _ in labels] for _ in labels]
    tm_matrix = [[1.0 if row == column else 0.0 for column in range(len(labels))]
                 for row in range(len(labels))]
    for result in results:
        row = index[result["structure_1"]]
        column = index[result["structure_2"]]
        rmsd_matrix[row][column] = rmsd_matrix[column][row] = result["rmsd_angstrom"]
        tm_matrix[row][column] = tm_matrix[column][row] = result["tm_score_norm_average_length"]
    write_matrix(output_dir / "pairwise_rmsd_matrix.csv", labels, rmsd_matrix, 2)
    write_matrix(output_dir / "pairwise_tm_score_matrix.csv", labels, tm_matrix, 5)
    figure.write_text(heatmap_svg(labels, rmsd_matrix, tm_matrix), encoding="utf-8")

    tm_distance_matrix = [
        [0.0 if row == column else 1.0 - tm_matrix[row][column]
         for column in range(len(labels))]
        for row in range(len(labels))
    ]
    rmsd_tree, rmsd_merges = upgma(labels, rmsd_matrix)
    tm_tree, tm_merges = upgma(labels, tm_distance_matrix)
    (output_dir / "upgma_rmsd.newick").write_text(to_newick(rmsd_tree), encoding="utf-8")
    (output_dir / "upgma_tm_distance.newick").write_text(to_newick(tm_tree), encoding="utf-8")
    merge_fields = [
        "step", "left_cluster", "right_cluster", "cluster_size",
        "average_linkage_distance", "node_height",
    ]
    for filename, merges in (
        ("upgma_rmsd_merges.csv", rmsd_merges),
        ("upgma_tm_distance_merges.csv", tm_merges),
    ):
        formatted_merges = []
        for merge in merges:
            formatted = dict(merge)
            formatted["average_linkage_distance"] = f'{merge["average_linkage_distance"]:.6f}'
            formatted["node_height"] = f'{merge["node_height"]:.6f}'
            formatted_merges.append(formatted)
        write_csv(output_dir / filename, formatted_merges, merge_fields)
    dendrogram_path = figure.parent / "pairwise_UPGMA_dendrograms.svg"
    dendrogram_path.write_text(dendrogram_svg(rmsd_tree, tm_tree), encoding="utf-8")

    rmsd_values = [result["rmsd_angstrom"] for result in results]
    tm_values = [result["tm_score_norm_average_length"] for result in results]
    same_fold_count = sum(value >= 0.5 for value in tm_values)
    summary = {
        "number_of_structures": len(labels),
        "number_of_unique_pairs": total_pairs,
        "rmsd_angstrom": {
            "minimum": min(rmsd_values),
            "median": statistics.median(rmsd_values),
            "maximum": max(rmsd_values),
        },
        "tm_score_average_length_normalized": {
            "minimum": min(tm_values),
            "median": statistics.median(tm_values),
            "maximum": max(tm_values),
            "pairs_at_least_0_5": same_fold_count,
            "fraction_at_least_0_5": same_fold_count / total_pairs,
        },
        "top_ten_pairs_by_tm_score": [
            {
                "structure_1": result["structure_1"],
                "structure_2": result["structure_2"],
                "tm_score": result["tm_score_norm_average_length"],
                "rmsd_angstrom": result["rmsd_angstrom"],
                "aligned_length": result["aligned_length"],
            }
            for result in sorted(
                results, key=lambda item: item["tm_score_norm_average_length"], reverse=True
            )[:10]
        ],
    }
    (output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    metadata = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "usalign_version": version,
        "usalign_executable_sha256": sha256(executable),
        "command_template": "USalign structure_1.pdb structure_2.pdb -mol prot -a T -outfmt -1",
        "tm_score_used_in_symmetric_matrix": "normalized by the average length of the two structures (-a T)",
        "rmsd_definition": "US-align RMSD over the residues in the final structural alignment",
        "upgma_rmsd_distance": "aligned-residue RMSD in Angstrom",
        "upgma_tm_distance": "1 - average-length-normalized TM-score",
        "upgma_linkage": "unweighted arithmetic mean of all inter-cluster pairwise distances",
    }
    (output_dir / "run_metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
