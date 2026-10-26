#!/usr/bin/env python3
"""Verify saved k=5 certificates using only the Python standard library.

This checker imports no audit/search code and never computes a fractional
maximum. A feasible residual plus an explicit blocking inequality proves
maximality. It re-enumerates every weight vector in each supplied graph's
reduced domain and requires exactly one certificate for every feasible vector.
Graph-isomorphism completeness is supplied by audit.cpp's exhaustive generator.

Run: python verify.py
"""

import itertools
import json
import time
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent
K = 5


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def connected(n, edges, omitted=None):
    vertices = set(range(n))
    if omitted is not None:
        vertices.remove(omitted)
    if not vertices:
        return True
    reached = {min(vertices)}
    while True:
        updated = set(reached)
        for u, v in edges:
            if u not in vertices or v not in vertices:
                continue
            if u in reached:
                updated.add(v)
            if v in reached:
                updated.add(u)
        if updated == reached:
            return reached == vertices
        reached = updated


def main():
    started = time.monotonic()
    graph_data = json.loads((ROOT / "graphs.json").read_text())
    records = defaultdict(dict)
    record_count = 0
    with (ROOT / "witnesses.jsonl").open() as stream:
        for line in stream:
            record = json.loads(line)
            graph_id = record["graph_id"]
            weights = record["weight_three_edge_mask"]
            require(weights not in records[graph_id], "Duplicate certificate")
            records[graph_id][weights] = record
            record_count += 1

    totals = Counter()
    by_vertices = defaultdict(Counter)
    coefficients = Counter()
    graph_masks = set()
    for expected_id, graph in enumerate(graph_data["graphs"]):
        require(graph["id"] == expected_id, "Nonsequential graph identifier")
        n = graph["n"]
        edges = [tuple(edge) for edge in graph["edges"]]
        m = len(edges)
        q = K * (n - 1) - 2 * m
        require(5 <= n <= 7 and 0 <= q <= m, "Graph outside reduced domain")
        require(edges == sorted(set(edges)), "Unsorted/duplicate edges")
        require(all(0 <= u < v < n for u, v in edges), "Invalid edge")
        require(graph["m"] == m and graph["q"] == q, "Wrong graph metadata")
        degrees = Counter(v for edge in edges for v in edge)
        require(all(degrees[v] >= 3 for v in range(n)), "Degree below three")
        require(connected(n, edges), "Disconnected graph")
        require(all(connected(n, edges, v) for v in range(n)), "Articulation vertex")
        graph_mask = sum(1 << i for i, edge in enumerate(itertools.combinations(range(n), 2)) if edge in edges)
        require(graph["canonical_mask"] == graph_mask, "Wrong graph mask")
        require((n, graph_mask) not in graph_masks, "Duplicate graph")
        graph_masks.add((n, graph_mask))

        subset_data = []
        for size in range(2, n):
            for vertices0 in itertools.combinations(range(n), size):
                vertices = set(vertices0)
                indices = tuple(i for i, (u, v) in enumerate(edges) if u in vertices and v in vertices)
                vertex_mask = sum(1 << v for v in vertices)
                subset_data.append((vertex_mask, size - 1, indices))

        local = Counter(graphs=1)
        certificates = records.pop(expected_id, {})
        for weight_three_indices in itertools.combinations(range(m), q):
            local["candidates"] += 1
            weights = [2] * m
            weight_mask = 0
            for i in weight_three_indices:
                weights[i] = 3
                weight_mask |= 1 << i
            require(sum(weights) == K * (n - 1), "Wrong input total")
            input_sums = [sum(weights[i] for i in indices) for _, _, indices in subset_data]
            feasible = all(value <= K * rank for value, (_, rank, _) in zip(input_sums, subset_data))
            if not feasible:
                require(weight_mask not in certificates, "Certificate for infeasible input")
                continue
            local["feasible"] += 1
            require(weight_mask in certificates, "Missing certificate")
            record = certificates.pop(weight_mask)
            coefficient = record["coefficient"]
            tree_mask = record["tree_edge_mask"]
            require(isinstance(coefficient, int) and 1 <= coefficient <= K, "Invalid coefficient")
            require(isinstance(tree_mask, int) and 0 <= tree_mask < (1 << m), "Invalid tree mask")
            selected = [bool(tree_mask & (1 << i)) for i in range(m)]
            tree_edges = [edge for i, edge in enumerate(edges) if selected[i]]
            require(len(tree_edges) == n - 1 and connected(n, tree_edges), "Not a spanning tree")

            scale = K - coefficient
            residual = [value - coefficient * used for value, used in zip(weights, selected)]
            require(all(z >= 0 for z in residual), "Negative residual")
            require(sum(residual) == scale * (n - 1), "Wrong residual total")
            blocker_kind = record["blocking_kind"]
            blocker_index = record["blocking_index"]
            stored_blocker_verified = blocker_kind == 0 and coefficient == K
            if blocker_kind == 1:
                require(0 <= blocker_index < m, "Invalid blocking edge index")
                stored_blocker_verified = selected[blocker_index] and residual[blocker_index] == 0
            any_blocker = coefficient == K or any(used and z == 0 for used, z in zip(selected, residual))
            for vertex_mask, rank, indices in subset_data:
                value = sum(residual[i] for i in indices)
                bound = scale * rank
                require(value <= bound, "Residual subset violation")
                blocks = value == bound and sum(selected[i] for i in indices) < rank
                any_blocker = any_blocker or blocks
                if blocker_kind == 2 and blocker_index == vertex_mask:
                    stored_blocker_verified = blocks
            require(any_blocker and stored_blocker_verified, "No valid maximality certificate")
            local["certified_witnesses"] += 1
            coefficients[coefficient] += 1

        require(not certificates, "Unused certificate records")
        for key in ("graphs", "candidates", "feasible", "certified_witnesses"):
            require(graph["counts"][key] == local[key], "Wrong per-graph reported count")
        totals.update(local)
        by_vertices[n].update(local)

    require(not records, "Certificates refer to unknown graphs")
    require(totals["certified_witnesses"] == record_count, "Record count mismatch")
    reported = json.loads((ROOT / "results.json").read_text())
    for key, value in totals.items():
        require(reported["total"][key] == value, "Aggregate count mismatch")
    for n, counts in by_vertices.items():
        for key, value in counts.items():
            require(reported["by_vertices"][str(n)][key] == value, "Vertex count mismatch")
    require(reported["total"]["coefficient_counts"] == {str(k): v for k, v in sorted(coefficients.items())},
            "Coefficient histogram mismatch")

    result = {
        "status": "all saved certificates and complete weight domains verified",
        "checker": "Python standard library; direct residual inequalities and blocking constraints",
        "ratio_search_code_imported": False,
        "graph_isomorphism_completeness": "provided by the separate exhaustive C++ graph generator",
        "by_vertices": {str(n): dict(counts) for n, counts in by_vertices.items()},
        "total": dict(totals),
        "coefficient_counts": dict(sorted(coefficients.items())),
        "runtime_seconds": time.monotonic() - started,
    }
    (ROOT / "verification.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
