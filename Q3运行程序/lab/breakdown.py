"""Cost breakdown harness for the Q3 policy (local synthetic simulator only).

Usage:
    python3 lab/breakdown.py --cases 6 --seed-base 2026091900
    python3 lab/breakdown.py --scenario uniform --sources 13 --seed 2026091900 --trace
"""

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import simulator  # noqa: E402
from robot_dog import execute_q3  # noqa: E402


def run_case(scenario, n, seed, trace=False):
    sources = simulator.sources_for(seed, scenario, n)
    transport = simulator.MockTransport(seed=seed, sources=sources, robot_id='local-test')
    log_path = Path('runs') / f'lab-q3-{scenario}-{seed}.jsonl'
    client = simulator.HTTPClient('local-test', transport=transport, log_path=log_path)
    result = execute_q3(client)
    try:
        report = transport.report()
    except Exception:
        report = {'source_count': len(sources), 'cleared_count': result['cleared_count']}
    row = dict(
        scenario=scenario,
        seed=seed,
        sources=report['source_count'],
        cleared=result['cleared_count'],
        status=result['status'],
        virtual_s=round(result['virtual_time_s'], 2),
        travel_m=round(result['travel_m'], 1),
        travel_s=round(result['travel_time_s'], 2),
        measure_s=round(result['measurement_time_s'], 2),
        measure_n=result['measure_count'],
        switch_s=round(result['switching_time_s'], 2),
        switch_n=result['switch_count'],
        clear_s=round(result['clearance_time_s'], 2),
        clear_n=result['clear_attempts'],
        actions=result['action_count'],
        census_points=len(result['census_points'] or []),
        cover_complete=result['area_coverage_complete'],
        unresolved_cells=result['unresolved_coverage_cells'],
    )
    if trace:
        row['phases'] = result['phase_time_s']
    return row


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('--cases', type=int, default=6)
    parser.add_argument('--seed-base', type=int, default=2026091900)
    parser.add_argument('--scenario', default=None)
    parser.add_argument('--sources', type=int, default=None)
    parser.add_argument('--seed', type=int, default=2026091900)
    parser.add_argument('--trace', action='store_true')
    args = parser.parse_args(argv)

    if args.scenario:
        rows = [run_case(args.scenario, args.sources, args.seed, args.trace)]
    else:
        rows = []
        for si, scenario in enumerate(('uniform', 'boundary', 'cluster')):
            for i in range(args.cases):
                rows.append(run_case(scenario, (10, 13, 16)[i % 3], args.seed_base + si * 1000 + i, args.trace))

    for row in rows:
        print(json.dumps(row, ensure_ascii=False))

    keys = ('virtual_s', 'travel_s', 'measure_s', 'switch_s', 'clear_s')
    summary = {k: round(sum(r[k] for r in rows) / len(rows), 2) for k in keys}
    summary['travel_share'] = round(sum(r['travel_s'] for r in rows) / sum(r['virtual_s'] for r in rows), 3)
    summary['census_points_avg'] = round(sum(r['census_points'] for r in rows) / len(rows), 2)
    groups = {}
    for row in rows:
        groups.setdefault(row['scenario'], []).append(row)
    print(json.dumps({'overall': summary,
                      'by_scenario': {k: round(sum(r['virtual_s'] for r in v) / len(v), 1) for k, v in groups.items()}},
                     ensure_ascii=False))
    return 0


if __name__ == '__main__':
    sys.exit(main())
