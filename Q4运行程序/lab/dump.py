"""Print one q4 trace action by action, with geometry context.

Usage (from the delivery folder):
    python3 lab/dump.py runs/base-mixed-2026092100.jsonl [--start 0] [--limit 60]

Ground truth from the simulator is used for offline diagnosis only.
"""

import argparse
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import simulator as S  # noqa: E402
from robot_dog import q4_RingCover  # noqa: E402

STATIONS = q4_RingCover().sites


def load_actions(path):
    pending = {}
    actions = []
    for line in Path(path).read_text(encoding='utf-8').splitlines():
        if not line.strip():
            continue
        record = json.loads(line)
        if record.get('event') == 'request':
            pending[record['payload']['request_id']] = record
        elif record.get('event') == 'response':
            request = pending.pop(record['request_id'], None)
            if request is not None:
                actions.append((request, record.get('response', record)))
    return actions


def nearest_station(point):
    index = min(range(len(STATIONS)), key=lambda i: math.dist(STATIONS[i], point))
    return index, math.dist(STATIONS[index], point)


def main(argv=None):
    parser = argparse.ArgumentParser(prog='lab/dump.py')
    parser.add_argument('path')
    parser.add_argument('--start', type=int, default=0)
    parser.add_argument('--limit', type=int, default=45)
    args = parser.parse_args(argv)

    name = Path(args.path).name
    parts = name.split('-')
    scenario, seed = parts[1], int(parts[-1].split('.')[0])
    sources = {int(s['channel']): tuple(s['position']) for s in S.sources_for_q4(seed, scenario, None)}

    previous = (0.0, 0.0)
    cumulative = 0.0
    cleared = set()
    print(f'# {args.path} scenario={scenario} sources={len(sources)}')
    print(f"{'#':>4} {'action':>6} {'ch':>3} {'position':>18} {'leg':>7} {'cum':>8} "
          f"{'result':>12} {'nearest_station':>16} {'dist_to_truth':>13}")
    for number, (request, response) in enumerate(load_actions(args.path)):
        kind = request.get('path')
        if kind not in ('/measure', '/clear'):
            continue
        payload = request['payload']
        point = (float(payload['position']['x']), float(payload['position']['y']))
        channel = int(payload['channel'])
        leg = math.dist(previous, point)
        cumulative += leg / 5 + (5 + int(channel != 1) if kind == '/measure' else 3)
        previous = point
        index, gap = nearest_station(point)
        truth = sources.get(channel)
        truth_gap = math.dist(truth, point) if truth else float('nan')
        if kind == '/measure':
            result = response['measure_result']
        else:
            result = 'clear:' + response['clear_result']
            if response['clear_result'] == 'success':
                cleared.add(channel)
        if number < args.start or number >= args.start + args.limit:
            continue
        print(f'{number:4d} {kind[1:]:>6} {channel:3d} '
              f'({point[0]:8.1f},{point[1]:8.1f}) {leg:7.1f} {cumulative:8.1f} {result:>12} '
              f'{index:5d} ({gap:6.1f}m) {truth_gap:13.1f}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
