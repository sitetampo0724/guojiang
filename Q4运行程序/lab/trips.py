"""Break every trace into maximal off-station detours and price each one.

Usage (from the delivery folder):
    python3 lab/trips.py runs/base-*.jsonl

A trip starts when the robot leaves the last visited census station and ends
when it reaches the next census station.  Its cost is the walked length minus
the straight chord between the two stations, i.e. the travel the station tour
would not have paid anyway.
"""

import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import simulator as S  # noqa: E402
from robot_dog import q4_RingCover  # noqa: E402

STATIONS = q4_RingCover().sites
STATION_INDEX = {site: i for i, site in enumerate(STATIONS)}


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


def trips2(path):
    """Split the trace into off-station trips with entry/exit stations."""
    events = []
    for request, response in load_actions(path):
        if request.get('path') not in ('/measure', '/clear'):
            continue
        payload = request['payload']
        point = (float(payload['position']['x']), float(payload['position']['y']))
        events.append(dict(kind=request['path'], point=point, channel=int(payload['channel']),
                           result=(response.get('measure_result')
                                   or 'clear:' + response.get('clear_result')),
                           station=STATION_INDEX.get(point)))
    position = (0.0, 0.0)
    out = []
    last_station = None
    trip = None
    for event in events:
        leg = math.dist(position, event['point'])
        if event['station'] is not None:
            if trip is not None:
                trip['walk_m'] += leg
                chord = (math.dist(trip['entry_point'], event['point'])
                         if trip['entry_point'] is not None else 0.0)
                trip['chord_m'] = chord
                trip['extra_m'] = trip['walk_m'] - chord
                trip['exit_station'] = event['station']
                out.append(trip)
                trip = None
            last_station = event['station']
            position = event['point']
            continue
        if trip is None:
            trip = dict(entry_station=last_station, entry_point=position if last_station is not None else None,
                        walk_m=0.0, chord_m=0.0, extra_m=0.0, count=0,
                        channels=[], results=[], points=[])
        trip['walk_m'] += leg
        trip['points'].append(event['point'])
        trip['count'] += 1
        trip['channels'].append(event['channel'])
        trip['results'].append(event['result'])
        position = event['point']
    if trip is not None:
        trip['exit_station'] = None
        trip['chord_m'] = 0.0
        trip['extra_m'] = trip['walk_m']
        out.append(trip)
    return out


def main(argv):
    paths = []
    for arg in argv:
        paths.extend(sorted(Path().glob(arg)) if any(c in arg for c in '*?[') else [Path(arg)])
    if not paths:
        print(__doc__)
        return 1
    totals = dict(extra=0.0, walk=0.0, trips=0, sources=0, wasted=0.0, wasted_trips=0)
    per_scenario = {}
    for path in paths:
        name = Path(path).name
        scenario = name.split('-')[1]
        seed = int(name.split('-')[-1].split('.')[0])
        truth = {int(s['channel']): tuple(s['position'])
                 for s in S.sources_for_q4(seed, scenario, None)}
        info = per_scenario.setdefault(scenario, dict(extra=0.0, walk=0.0, trips=0,
                                                      sources=0, wasted=0.0, wasted_trips=0,
                                                      cases=0))
        info['cases'] += 1
        for trip in trips2(path):
            totals['extra'] += trip['extra_m']
            totals['walk'] += trip['walk_m']
            totals['trips'] += 1
            info['extra'] += trip['extra_m']
            info['walk'] += trip['walk_m']
            info['trips'] += 1
            failed = [i for i, r in enumerate(trip['results']) if r != 'clear:success'
                      and not r.startswith('direction')]
            if any(r == 'clear:no_target_in_range' for r in trip['results']):
                totals['wasted'] += trip['walk_m']
                totals['wasted_trips'] += 1
                info['wasted'] += trip['walk_m']
                info['wasted_trips'] += 1
            _ = failed
        info['sources'] += len(truth)
        totals['sources'] += len(truth)
    print(f"{'scenario':16s} {'cases':>5s} {'trips':>6s} {'avg_extra_m':>12s} "
          f"{'avg_walk_m':>11s} {'failed_trips':>12s} {'wasted_avg_m':>12s}")
    for scenario, info in sorted(per_scenario.items()):
        n = info['cases']
        print(f"{scenario:16s} {n:5d} {info['trips'] / n:6.1f} {info['extra'] / n:12.0f} "
              f"{info['walk'] / n:11.0f} {info['wasted_trips'] / n:12.1f} "
              f"{info['wasted'] / n:12.0f}")
    n = len(paths)
    print(f"{'MEAN':16s} {n:5d} {totals['trips'] / n:6.1f} {totals['extra'] / n:12.0f} "
          f"{totals['walk'] / n:11.0f} {totals['wasted_trips'] / n:12.1f} "
          f"{totals['wasted'] / n:12.0f}")
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
