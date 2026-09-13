"""Per-action trace of the Q3 policy: where the virtual seconds actually go.

Usage:
    python3 lab/trace.py --scenario uniform --sources 10 --seed 2026091900
    python3 lab/trace.py --scenario boundary --sources 13 --seed 2026092901 --dump results/trace_boundary.json
"""

import argparse
import json
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import simulator  # noqa: E402
import robot_dog  # noqa: E402


def instrument(cls):
    original_measure = cls.measure
    original_clear = cls.clear
    original_account = cls._account

    def _account(self, point, response, fixed_cost):
        before = self.virtual_time
        original_account(self, point, response, fixed_cost)
        self.trace[-1]['cost_s'] = round(self.virtual_time - before, 3)
        self.trace[-1]['virtual_s'] = round(self.virtual_time, 3)

    def measure(self, point, channel):
        distance = math.dist(self.position, point)
        self.trace.append(dict(kind='measure', phase=self.phase, channel=channel,
                               at=[round(v, 1) for v in point], travel_m=round(distance, 1),
                               switched=int(channel != self.channel), t=self.virtual_time))
        response = original_measure(self, point, channel)
        self.trace[-1]['result'] = response['measure_result']
        return response

    def clear(self, point, channel):
        distance = math.dist(self.position, point)
        self.trace.append(dict(kind='clear', phase=self.phase, channel=channel,
                               at=[round(v, 1) for v in point], travel_m=round(distance, 1), t=self.virtual_time))
        result = original_clear(self, point, channel)
        self.trace[-1]['result'] = 'success' if result else 'fail'
        return result

    cls.trace = []
    cls._account = _account
    cls.measure = measure
    cls.clear = clear


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('--scenario', default='uniform')
    parser.add_argument('--sources', type=int, default=10)
    parser.add_argument('--seed', type=int, default=2026091900)
    parser.add_argument('--dump', default=None)
    parser.add_argument('--quiet', action='store_true')
    args = parser.parse_args(argv)

    instrument(robot_dog.Controller)

    sources = simulator.sources_for(args.seed, args.scenario, args.sources)
    transport = simulator.MockTransport(seed=args.seed, sources=sources, robot_id='local-test')
    Path('runs').mkdir(exist_ok=True)
    client = simulator.HTTPClient('local-test', transport=transport, log_path=Path('runs') / f'trace-q3-{args.scenario}-{args.seed}.jsonl')
    result = robot_dog.execute_q3(client)
    report = transport.report()
    trace = robot_dog.Controller.trace

    by_phase = defaultdict(lambda: dict(cost=0.0, travel=0.0, measures=0, clears=0))
    for row in trace:
        slot = by_phase[row['phase']]
        slot['cost'] += row['cost_s']
        slot['travel'] += row['travel_m']
        slot['measures' if row['kind'] == 'measure' else 'clears'] += 1

    print(json.dumps(dict(scenario=args.scenario, seed=args.seed, n=len(sources),
                          virtual_s=round(result['virtual_time_s'], 1),
                          travel_m=round(result['travel_m'], 1),
                          actions=len(trace),
                          by_phase={k: dict(cost=round(v['cost'], 1), travel_m=round(v['travel'], 1),
                                            measures=v['measures'], clears=v['clears'])
                                    for k, v in by_phase.items()},
                          phase_time_s=result['phase_time_s']),
                     ensure_ascii=False, indent=1))

    if not args.quiet:
        for row in trace:
            print(('{:>7.1f}s {:<22} {:>7} ch={:<3} travel={:>7} {:>11} {} -> {}'.format(
                row['t'], row['phase'], row['kind'], row['channel'], row['travel_m'],
                row.get('result', ''), row['at'], row['cost_s'])))
    if args.dump:
        Path(args.dump).write_text(json.dumps(
            dict(trace=trace, sources=[(round(s['position'][0], 1), round(s['position'][1], 1), s['channel'], round(s['radius'], 1)) for s in report['sources']],
                 virtual_s=result['virtual_time_s'], travel_m=result['travel_m']), ensure_ascii=False, indent=1))
    return 0


if __name__ == '__main__':
    sys.exit(main())
