"""Classify every action of recorded q4 traces and aggregate the cost model.

Usage (from the delivery folder):
    python3 lab/trace.py runs/base-mixed-2026092100.jsonl
    python3 lab/trace.py 'runs/base-*.jsonl'

Virtual time is split into its four real components (travel, measurement,
channel switching, clearance) and each measurement is labelled by purpose:

    station/certificate  a census station where the channel had never answered
    station/pending      a census station re-measuring a channel already located
    station/cleared      a census station re-measuring a channel already cleared
    detour/*             the same split for actions taken off the 25 census sites
"""

import json
import math
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from robot_dog import q4_RingCover  # noqa: E402

STATIONS = q4_RingCover().sites


def load_actions(path):
    """Pair every request with its response, in order."""
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


def station_key(point, tol=1e-6):
    for index, site in enumerate(STATIONS):
        if math.dist(point, site) < tol:
            return index
    return None


def analyse(path):
    actions = load_actions(path)
    answered = set()      # channel produced a direction/near reading somewhere
    cleared = set()
    position = (0.0, 0.0)
    channel_now = 1
    buckets = defaultdict(lambda: dict(n=0, travel_m=0.0, travel_s=0.0, fixed_s=0.0,
                                       switch=0, measure=0, clear=0))
    station_hits = set()
    ring_travel = 0.0
    total = 0.0
    previous_station = None
    clear_events = []

    for request, response in actions:
        kind = request.get('path') or response.get('path')
        if kind not in ('/measure', '/clear'):
            continue
        payload = request['payload']
        point = (float(payload['position']['x']), float(payload['position']['y']))
        channel = int(payload['channel'])
        leg = math.dist(position, point)
        total += leg / 5.0
        index = station_key(point)
        at_station = index is not None
        if at_station:
            station_hits.add(index)
            if previous_station is not None:
                ring_travel += math.dist(STATIONS[previous_station], point)
            previous_station = index

        if kind == '/measure':
            result = response['measure_result']
            if channel in cleared:
                why = 'cleared'
            elif channel in answered:
                why = 'pending'
            else:
                why = 'certificate'
            switch = int(channel != channel_now)
            channel_now = channel
            bucket = buckets[f"{'station' if at_station else 'detour'}/{why}"]
            bucket['n'] += 1
            bucket['measure'] += 1
            bucket['switch'] += switch
            bucket['fixed_s'] += 5 + switch
            total += 5 + switch
            if result in ('direction', 'near'):
                answered.add(channel)
        else:
            result = response['clear_result']
            success = result == 'success'
            bucket = buckets[f"{'station' if at_station else 'detour'}/clear-{result}"]
            bucket['n'] += 1
            bucket['clear'] += 1
            bucket['fixed_s'] += 3 + 2 * success
            total += 3 + 2 * success
            if success:
                cleared.add(channel)
                answered.discard(channel)
            clear_events.append(dict(channel=channel, ok=success, point=point,
                                     leg=leg, travel_s=leg / 5))
        bucket['travel_m'] += leg
        bucket['travel_s'] += leg / 5.0
        position = point

    travel_m = sum(b['travel_m'] for b in buckets.values())
    fixed_s = sum(b['fixed_s'] for b in buckets.values())
    return dict(
        path=str(path),
        total_s=total,
        travel_m=travel_m,
        travel_s=travel_m / 5.0,
        fixed_s=fixed_s,
        ring_travel_m=ring_travel,
        detour_m=travel_m - ring_travel,
        stations_visited=len(station_hits),
        measures=sum(b['measure'] for b in buckets.values()),
        switches=sum(b['switch'] for b in buckets.values()),
        clear_attempts=len(clear_events),
        clears_ok=sum(1 for e in clear_events if e['ok']),
        buckets={k: v for k, v in sorted(buckets.items())},
        clear_events=clear_events,
    )


def print_trace(info):
    print(f"== {info['path']}")
    print(f"   total={info['total_s']:.0f}s travel={info['travel_m']:.0f}m "
          f"({info['travel_s']:.0f}s) fixed={info['fixed_s']:.0f}s "
          f"stations={info['stations_visited']}/{len(STATIONS)} "
          f"ring_leg={info['ring_travel_m']:.0f}m detour={info['detour_m']:.0f}m "
          f"measures={info['measures']} switches={info['switches']} "
          f"clears={info['clears_ok']}/{info['clear_attempts']}")
    for key, bucket in info['buckets'].items():
        print(f"   {key:26s} n={bucket['n']:4d} t={bucket['fixed_s']:7.1f}s "
              f"sw={bucket['switch']:3d} walk={bucket['travel_m']:7.0f}m")


def aggregate(infos):
    totals = defaultdict(float)
    buckets = defaultdict(lambda: defaultdict(float))
    core = ('total_s', 'travel_m', 'travel_s', 'fixed_s', 'ring_travel_m', 'detour_m',
            'measures', 'switches', 'clears_ok', 'clear_attempts', 'stations_visited')
    for info in infos:
        for key in core:
            totals[key] += info[key]
        for name, bucket in info['buckets'].items():
            for key in ('n', 'travel_m', 'travel_s', 'fixed_s', 'switch', 'measure', 'clear'):
                buckets[name][key] += bucket[key]
    n = len(infos)
    print(f"== aggregate over {n} traces")
    print(f"   avg total={totals['total_s'] / n:.0f}s travel={totals['travel_m'] / n:.0f}m "
          f"({totals['travel_s'] / n:.0f}s) fixed={totals['fixed_s'] / n:.0f}s "
          f"stations={totals['stations_visited'] / n:.1f} "
          f"ring_leg={totals['ring_travel_m'] / n:.0f}m detour={totals['detour_m'] / n:.0f}m "
          f"measures={totals['measures'] / n:.1f} switches={totals['switches'] / n:.1f} "
          f"clears={totals['clears_ok'] / n:.1f}/{totals['clear_attempts'] / n:.1f}")
    for name in sorted(buckets):
        bucket = buckets[name]
        print(f"   {name:26s} n/case={bucket['n'] / n:6.1f} t={bucket['fixed_s'] / n:7.1f}s "
              f"sw={bucket['switch'] / n:5.1f} walk={bucket['travel_m'] / n:7.0f}m")


def main(argv):
    paths = []
    for arg in argv:
        paths.extend(sorted(Path().glob(arg)) if any(c in arg for c in '*?[') else [Path(arg)])
    if not paths:
        print(__doc__)
        return 1
    infos = []
    for path in paths:
        info = analyse(path)
        infos.append(info)
        if len(paths) <= 4:
            print_trace(info)
    if len(infos) > 1:
        aggregate(infos)
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
