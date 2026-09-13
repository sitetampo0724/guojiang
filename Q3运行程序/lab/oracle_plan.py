"""Diagnostic: how good would the station plan be if every source were known
up front?  (Oracle - not a deliverable policy.)  Compares the ring plan with
the source-visit plan on identical information.
"""

import argparse
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import robot_dog  # noqa: E402
import simulator  # noqa: E402


class FakeCtrl:
    def __init__(self, position):
        self.position = position
        self.cleared = set()
        self.source_positions = {}


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('--cases', type=int, default=6)
    parser.add_argument('--seed-base', type=int, default=2026091900)
    args = parser.parse_args(argv)
    for si, scenario in enumerate(('uniform', 'boundary', 'cluster')):
        for i in range(args.cases):
            n = (10, 13, 16)[i % 3]
            seed = args.seed_base + si * 1000 + i
            raw = simulator.sources_for(seed, scenario, n)
            sources = {s['channel']: tuple(s['position']) for s in raw}
            cover = robot_dog.AreaCover()
            cover.add((0.0, 0.0))
            grid = robot_dog.StationGrid()
            covered = grid.mask((0.0, 0.0))
            source_route = [sources[c] for c in robot_dog.warm_route((0.0, 0.0), sources)]
            scan_cost = 6.0 * (20 - n)
            terms = (0.0, None)
            ring = robot_dog.choose_plan_joint(FakeCtrl((0.0, 0.0)), cover, sources, list(range(1, 21 - n)))
            options = robot_dog.ring_station_options(covered, cover.points, grid, source_route)
            best_ring = None
            for option in options:
                stations = [(p, 'ring') for p in option['stations']]
                cost = robot_dog.station_plan_cost((0.0, 0.0), source_route, stations, scan_cost, terms)
                if best_ring is None or cost < best_ring[0]:
                    best_ring = (cost, option)
            print(json.dumps(dict(scenario=scenario, seed=seed, n=n,
                                  best_ring_s=round(best_ring[0], 1) if best_ring else None,
                                  best_ring_radius=best_ring[1]['radius'] if best_ring else None,
                                  joint_s=None if ring is None else round(ring['estimated_s'], 1),
                                  joint_stations=None if ring is None else ring['stations_planned'],
                                  joint_kinds=None if ring is None else ring['station_kinds']),
                             ensure_ascii=False))
    return 0


if __name__ == '__main__':
    sys.exit(main())
