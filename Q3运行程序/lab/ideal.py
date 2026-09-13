"""Offline reference: exact open-route TSP over (census ring + source sites).

This is NOT part of the delivered program. It estimates how much travel a
joint census/source route could plausibly need, so we can tell whether a
candidate planner change is worth implementing.

Usage:
    python3 lab/ideal.py --scenario uniform --sources 10 --seed 2026091900
    python3 lab/ideal.py --cases 6 --seed-base 2026091900
"""

import argparse
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import simulator  # noqa: E402


def covering_radius(a, R=1800.0):
    if a >= math.sqrt(3) * R:
        return R
    return max(a / math.sqrt(3), math.sqrt(R * R + a * a - math.sqrt(3) * R * a))


def tsp_exact(start, points):
    """Exact open route: start fixed, visit all points, free end."""
    n = len(points)
    if n == 0:
        return 0.0
    size = 1 << n
    inf = float('inf')
    dp = [[inf] * n for _ in range(size)]
    for j in range(n):
        dp[1 << j][j] = math.dist(start, points[j])
    for mask in range(size):
        row = dp[mask]
        for j in range(n):
            base = row[j]
            if base == inf:
                continue
            pj = points[j]
            for k in range(n):
                if mask >> k & 1:
                    continue
                new = base + math.dist(pj, points[k])
                if new < dp[mask | 1 << k][k]:
                    dp[mask | 1 << k][k] = new
    full = size - 1
    return min(dp[full])


def tsp_heuristic(start, points, restarts=3):
    """Nearest neighbour + 2-opt + or-opt, several starting choices."""
    points = list(points)
    n = len(points)
    if n == 0:
        return 0.0
    best = float('inf')
    for r in range(restarts):
        remaining = set(range(n))
        order = []
        cur = start
        first = min(remaining, key=lambda j: math.dist(cur, points[j])) if r == 0 else list(remaining)[r % n]
        order.append(first)
        remaining.remove(first)
        while remaining:
            nxt = min(remaining, key=lambda j: math.dist(points[order[-1]], points[j]))
            order.append(nxt)
            remaining.remove(nxt)

        def length(route):
            seq = [start] + [points[j] for j in route]
            return sum(math.dist(a, b) for a, b in zip(seq, seq[1:]))

        improved = True
        while improved:
            improved = False
            for i in range(n):
                for j in range(i + 1, n):
                    trial = order[:i] + order[i:j + 1][::-1] + order[j + 1:]
                    if length(trial) < length(order) - 1e-9:
                        order = trial
                        improved = True
            for i in range(n):
                for j in range(n):
                    if i == j or i + 1 == j:
                        continue
                    node = order[i]
                    rest = order[:i] + order[i + 1:]
                    insert = rest[:j] + [node] + rest[j:]
                    if length(insert) < length(order) - 1e-9:
                        order = insert
                        improved = True
        best = min(best, length(order))
    return best


def tsp(start, points):
    return tsp_exact(start, points) if len(points) <= 14 else tsp_heuristic(start, points)


def best_ring_route(start, sources, radii, steps=12):
    best = (float('inf'), None)
    for a in radii:
        if covering_radius(a) > 999.9:
            continue
        for k in range(steps):
            phi = math.pi * k / (3 * steps)
            ring = [(a * math.cos(phi + i * math.pi / 3), a * math.sin(phi + i * math.pi / 3)) for i in range(6)]
            length = tsp(start, ring + list(sources))
            if length < best[0]:
                best = (length, dict(radius=a, phi_deg=math.degrees(phi)))
    return best


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('--scenario', default='uniform')
    parser.add_argument('--sources', type=int, default=10)
    parser.add_argument('--seed', type=int, default=2026091900)
    parser.add_argument('--cases', type=int, default=0)
    parser.add_argument('--seed-base', type=int, default=2026091900)
    args = parser.parse_args(argv)

    radii = (1123.2, 1200.0, 1400.0, 1600.0)
    cases = []
    if args.cases:
        for si, scenario in enumerate(('uniform', 'boundary', 'cluster')):
            for i in range(args.cases):
                cases.append((scenario, (10, 13, 16)[i % 3], args.seed_base + si * 1000 + i))
    else:
        cases.append((args.scenario, args.sources, args.seed))

    out = []
    for scenario, n, seed in cases:
        raw = simulator.sources_for(seed, scenario, n)
        sources = [tuple(s['position']) for s in raw]
        source_only = tsp((0.0, 0.0), sources)
        ring_length, ring_info = best_ring_route((0.0, 0.0), sources, radii)
        row = dict(scenario=scenario, seed=seed, n=len(sources),
                   source_only_tsp_m=round(source_only, 1),
                   best_ring_route_m=round(ring_length, 1),
                   ring=ring_info)
        out.append(row)
        print(json.dumps(row, ensure_ascii=False))

    print(json.dumps(dict(
        source_only_avg_s=round(sum(r['source_only_tsp_m'] for r in out) / len(out) / 5, 1),
        best_ring_avg_s=round(sum(r['best_ring_route_m'] for r in out) / len(out) / 5, 1)), ensure_ascii=False))
    return 0


if __name__ == '__main__':
    sys.exit(main())
