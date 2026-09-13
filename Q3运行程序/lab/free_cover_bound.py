"""Lower bound: best open route over {origin} + 6 free census stations + sources.

The stations are chosen freely (not restricted to a ring) subject to covering
the arena with 999.9 m disks.  This is still optimistic: it ignores the fact
that a source can only be visited after it has been detected.

Usage:
    python3 lab/free_cover_bound.py --scenario boundary --sources 10 --seed 2026092903
"""

import argparse
import math
import random
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import robot_dog  # noqa: E402
import simulator  # noqa: E402


def make_grid():
    return robot_dog.StationGrid(step=40.0)


def covered_by(grid, points, radius):
    """Distance-based coverage on the sample grid (faster than bitmasks here)."""
    radius2 = radius * radius
    for x, y in grid.points:
        ok = False
        for px, py in points:
            if (x - px) ** 2 + (y - py) ** 2 <= radius2:
                ok = True
                break
        if not ok:
            return False
    return True


def route_length(start, points):
    return robot_dog.fast_route_length(start, points)


def optimise(sources, rng, restarts=6, iterations=400):
    grid = make_grid()
    best = None
    for restart in range(restarts):
        if restart == 0:
            stations = [(1123.2 * math.cos(k * math.pi / 3), 1123.2 * math.sin(k * math.pi / 3)) for k in range(6)]
        else:
            stations = [(rng.uniform(-1800, 1800), rng.uniform(-1800, 1800)) for _ in range(6)]
        current = route_length((0.0, 0.0), stations + list(sources))
        for _ in range(iterations):
            i = rng.randrange(6)
            old = stations[i]
            scale = 900.0 * (0.5 ** (rng.random() * 4))
            candidate = (old[0] + rng.gauss(0, scale), old[1] + rng.gauss(0, scale))
            if math.hypot(*candidate) > 1900:
                continue
            stations[i] = candidate
            if not covered_by(grid, [(0.0, 0.0)] + stations, 999.9):
                stations[i] = old
                continue
            value = route_length((0.0, 0.0), stations + list(sources))
            if value <= current:
                current = value
            else:
                stations[i] = old
        if best is None or current < best[0]:
            best = (current, list(stations))
    return best


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('--scenario', default='boundary')
    parser.add_argument('--sources', type=int, default=10)
    parser.add_argument('--seed', type=int, default=2026092903)
    parser.add_argument('--cases', type=int, default=0)
    parser.add_argument('--seed-base', type=int, default=2026091900)
    args = parser.parse_args(argv)
    rng = random.Random(12345)
    cases = []
    if args.cases:
        for si, scenario in enumerate(('uniform', 'boundary', 'cluster')):
            for i in range(args.cases):
                cases.append((scenario, (10, 13, 16)[i % 3], args.seed_base + si * 1000 + i))
    else:
        cases.append((args.scenario, args.sources, args.seed))
    for scenario, n, seed in cases:
        sources = [tuple(s['position']) for s in simulator.sources_for(seed, scenario, n)]
        length, stations = optimise(sources, rng)
        print('%s n=%d seed=%d  bounded route %.0f m  (%.1f s travel)  stations %s' % (
            scenario, n, seed, length, length / 5.0,
            ' '.join('(%.0f,%.0f)' % s for s in stations)))
    return 0


if __name__ == '__main__':
    sys.exit(main())
