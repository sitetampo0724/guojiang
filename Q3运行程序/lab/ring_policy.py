"""Controlled experiment: force a fixed census ring (radius, rotation) and let a
plain receding-horizon DP choose between scanning a ring point and clearing a
pending source.  Used to measure how the ring geometry affects travel.

Usage:
    python3 lab/ring_policy.py --radius 1150 --cases 4 --seed-base 2026091900
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


def install(radius, rotation_deg=0.0, skip_origin=False):
    rotation = math.radians(rotation_deg)
    ring = [(radius * math.cos(rotation + k * math.pi / 3), radius * math.sin(rotation + k * math.pi / 3)) for k in range(6)]
    if not skip_origin:
        ring = [(0.0, 0.0)] + ring

    def plan(ctrl, cover, original, area_planner):
        infos = {c: robot_dog.adaptive_summary(ctrl, c) for c in ctrl.pending()}
        sources = {c: info['center'] for c, info in infos.items()}
        unknown = [c for c in range(1, 21) if c not in ctrl.cleared and c not in ctrl.observations]
        tasks = dict(sources)
        kinds = {}
        for index, point in enumerate(ring):
            if all(math.dist(point, q) > 1.0 for q in cover.points):
                tasks[21 + index] = point
                kinds[21 + index] = 'ring'
        scan_cost = 6.0 * len(unknown)
        route, stats = robot_dog.solve_route(ctrl.position, tasks, horizon=16)
        return dict(centers=tasks, route=route, dp=stats,
                    estimated_s=stats['route_distance_m'] / 5.0 + len(kinds) * scan_cost,
                    plan_type='static_ring', scan_tasks=list(kinds), scan_task_kinds=kinds,
                    beliefs=infos, unknown_channels=unknown, reference_estimated_s=0.0,
                    alternatives=[], beam_states=0, circle_candidate_count=0,
                    anticipated_scan_sources=[])

    robot_dog.choose_plan = plan


def run_case(scenario, n, seed):
    sources = simulator.sources_for(seed, scenario, n)
    transport = simulator.MockTransport(seed=seed, sources=sources, robot_id='local-test')
    Path('runs').mkdir(exist_ok=True)
    client = simulator.HTTPClient('local-test', transport=transport, log_path=Path('runs') / f'ring-q3-{scenario}-{seed}.jsonl')
    result = robot_dog.execute_q3(client)
    return dict(scenario=scenario, seed=seed, n=n, virtual_s=round(result['virtual_time_s'], 1),
                travel_m=round(result['travel_m'], 1), status=result['status'],
                stations=len(result['census_points'] or []))


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('--radius', type=float, default=1150.0)
    parser.add_argument('--rotation', type=float, default=0.0)
    parser.add_argument('--skip-origin', action='store_true')
    parser.add_argument('--cases', type=int, default=4)
    parser.add_argument('--seed-base', type=int, default=2026091900)
    args = parser.parse_args(argv)
    install(args.radius, args.rotation, args.skip_origin)
    rows = []
    for si, scenario in enumerate(('uniform', 'boundary', 'cluster')):
        for i in range(args.cases):
            rows.append(run_case(scenario, (10, 13, 16)[i % 3], args.seed_base + si * 1000 + i))
    for row in rows:
        print(json.dumps(row, ensure_ascii=False))
    by = {}
    for row in rows:
        by.setdefault(row['scenario'], []).append(row['virtual_s'])
    print(json.dumps({'radius': args.radius, 'rotation': args.rotation,
                      'avg': {k: round(sum(v) / len(v), 1) for k, v in by.items()},
                      'ok': all(r['status'] == 'complete' for r in rows)}, ensure_ascii=False))
    return 0


if __name__ == '__main__':
    sys.exit(main())
