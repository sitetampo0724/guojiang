"""Run q4 cases and keep the full result breakdown.

Usage (from the delivery folder):
    python3 lab/bench.py --cases 6 --seed-base 2026092100 --tag base
    python3 lab/bench.py --scenarios mixed cluster --cases 2 --tag smoke

Writes one JSONL trace per case into runs/ and a summary JSON into results/.
"""

import argparse
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import simulator as S  # noqa: E402
import robot_dog as R  # noqa: E402
from robot_dog import HTTPClient, execute_q4  # noqa: E402


def run_case(scenario, seed, tag, log_dir=Path('runs'), planner=None):
    sources = S.sources_for_q4(seed, scenario, None)
    transport = S.MockTransport4(seed=seed, sources=sources, robot_id='local-test')
    log_path = Path(log_dir) / f'{tag}-{scenario}-{seed}.jsonl'
    client = HTTPClient('local-test', transport=transport, log_path=log_path)
    result = execute_q4(client, planner=planner)
    try:
        report = transport.report()
    except Exception:
        report = {'source_count': len(sources), 'cleared_count': result['cleared_count']}
    result['scenario'] = scenario
    result['seed'] = seed
    result['source_count'] = report['source_count']
    result['cleared_true_count'] = report['cleared_count']
    result['ok'] = (result['status'] == 'complete'
                    and report['cleared_count'] == report['source_count'])
    result.pop('census_points', None)
    return result


def summarize(rows):
    groups = []
    ok = True
    for scenario in sorted({r['scenario'] for r in rows}):
        sub = [r for r in rows if r['scenario'] == scenario]
        good = [r for r in sub if r['ok']]
        ok = ok and len(good) == len(sub)
        n = len(good)

        def mean(key):
            return sum(r[key] for r in good) / n if n else None

        groups.append(dict(
            scenario=scenario,
            cases=len(sub),
            success=len(good),
            avg_virtual_time_s=mean('virtual_time_s'),
            avg_per_source_s=mean('average_time_per_cleared_s'),
            avg_travel_m=mean('travel_m'),
            avg_travel_time_s=mean('travel_time_s'),
            avg_measure_count=mean('measure_count'),
            avg_measurement_time_s=mean('measurement_time_s'),
            avg_switch_count=mean('switch_count'),
            avg_switching_time_s=mean('switching_time_s'),
            avg_clear_attempts=mean('clear_attempts'),
            avg_clearance_time_s=mean('clearance_time_s'),
            avg_action_count=mean('action_count'),
            max_real_runtime_s=max((r['real_runtime_s'] for r in sub), default=None),
        ))
    all_good = [r for r in rows if r['ok']]
    overall = None
    if all_good:
        overall = dict(
            cases=len(rows),
            success=len(all_good),
            avg_virtual_time_s=sum(r['virtual_time_s'] for r in all_good) / len(all_good),
            avg_per_source_s=sum(r['average_time_per_cleared_s'] for r in all_good) / len(all_good),
            avg_travel_m=sum(r['travel_m'] for r in all_good) / len(all_good),
            avg_measure_count=sum(r['measure_count'] for r in all_good) / len(all_good),
            avg_switch_count=sum(r['switch_count'] for r in all_good) / len(all_good),
            avg_clear_attempts=sum(r['clear_attempts'] for r in all_good) / len(all_good),
            max_real_runtime_s=max(r['real_runtime_s'] for r in rows),
        )
    return dict(ok=ok, groups=groups, overall=overall, rows=rows)


def main(argv=None):
    parser = argparse.ArgumentParser(prog='lab/bench.py')
    parser.add_argument('--cases', type=int, default=6)
    parser.add_argument('--seed-base', type=int, default=2026092100)
    parser.add_argument('--scenarios', nargs='*', default=list(S.Q4_SCENARIOS))
    parser.add_argument('--tag', default='bench')
    parser.add_argument('--out', default=None)
    parser.add_argument('--log-dir', default='runs')
    parser.add_argument('--planner', default=None, help='v2 (default) or legacy')
    parser.add_argument('--tune', default=None, help='JSON dict overriding Q4_TUNING')
    args = parser.parse_args(argv)
    if args.tune:
        R.Q4_TUNING.update(json.loads(args.tune))
    planner = args.planner or R.Q4_TUNING['provider']
    rows = []
    started = time.perf_counter()
    for si, scenario in enumerate(args.scenarios):
        for i in range(args.cases):
            seed = args.seed_base + si * 1000 + i
            rows.append(run_case(scenario, seed, args.tag, Path(args.log_dir), planner))
    summary = summarize(rows)
    summary['tag'] = args.tag
    summary['wall_s'] = time.perf_counter() - started
    out = Path(args.out) if args.out else ROOT / 'results' / f'{args.tag}_bench.json'
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')
    for g in summary['groups']:
        print(f"{g['scenario']:16s} ok={g['success']}/{g['cases']} "
              f"avg={g['avg_virtual_time_s']:.1f}s travel={g['avg_travel_m']:.0f}m "
              f"meas={g['avg_measure_count']:.1f} switch={g['avg_switch_count']:.1f}")
    print(json.dumps(summary['overall'], ensure_ascii=False))
    print(f"wall={summary['wall_s']:.1f}s -> {out}")
    return 0 if summary['ok'] else 2


if __name__ == '__main__':
    sys.exit(main())
