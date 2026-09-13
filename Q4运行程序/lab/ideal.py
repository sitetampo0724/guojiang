"""Compare the recorded station tour and source detours against geometric ideals.

Usage (from the delivery folder):
    python3 lab/ideal.py runs/base-*.jsonl

For every trace we recover the order in which the 25 census stations were
visited and compare that leg length with the best structured tour, then compare
the extra travel spent reaching sources with the cheapest possible insertion of
those same sources into the station tour.
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


def path_legs(path):
    """Return (station_order, positions of every action, total travel)."""
    order = []
    positions = []
    total = 0.0
    previous = (0.0, 0.0)
    for request, _ in load_actions(path):
        if request.get('path') not in ('/measure', '/clear'):
            continue
        payload = request['payload']
        point = (float(payload['position']['x']), float(payload['position']['y']))
        total += math.dist(previous, point)
        positions.append(point)
        index = STATION_INDEX.get(point)
        if index is not None and (not order or order[-1] != index):
            order.append(index)
        previous = point
    return order, positions, total


def tour_length(order):
    points = [(0.0, 0.0)] + [STATIONS[i] for i in order]
    return sum(math.dist(a, b) for a, b in zip(points, points[1:]))


def structured_tours():
    """Candidate hand-built tours: which ring first, which direction."""
    inner = sorted(range(1, 9), key=lambda i: math.atan2(STATIONS[i][1], STATIONS[i][0]))
    outer = sorted(range(9, 25), key=lambda i: math.atan2(STATIONS[i][1], STATIONS[i][0]))
    best = None
    for first, second in ((inner, outer), (outer, inner)):
        for direction in (1, -1):
            seq = first[::direction] + second[::direction]
            for reverse_second in (False, True):
                candidate = first[::direction] + (second[::-direction] if reverse_second else second[::direction])
                length = tour_length(candidate)
                if best is None or length < best[0]:
                    best = (length, candidate)
            _ = seq
    return best


def insertion_cost(source, order):
    """Cheapest way to touch `source` while running the station tour."""
    points = [(0.0, 0.0)] + [STATIONS[i] for i in order]
    best = math.dist(points[-1], source)          # after the last station
    for a, b in zip(points, points[1:]):
        extra = math.dist(a, source) + math.dist(source, b) - math.dist(a, b)
        best = min(best, extra)
    return best


def main(argv):
    paths = []
    for arg in argv:
        paths.extend(sorted(Path().glob(arg)) if any(c in arg for c in '*?[') else [Path(arg)])
    if not paths:
        print(__doc__)
        return 1
    best_structured, structured_order = structured_tours()
    print(f'structured tour length = {best_structured:.0f} m '
          f'(open path from origin over all 25 stations)')
    columns = []
    for path in paths:
        name = Path(path).name
        parts = name.split('-')
        scenario = parts[1] if len(parts) > 2 else '?'
        seed = int(parts[-1].split('.')[0])
        sources = S.sources_for_q4(seed, scenario, None)
        order, positions, total = path_legs(path)
        seen = tour_length(order)
        ideal = sum(insertion_cost(tuple(src['position']), order) for src in sources)
        columns.append((scenario, seen + (len(STATIONS) - len(order)) * 0,
                        seen, ideal, total, len(order)))
    n = len(columns)
    print(f"{'scenario':16s} {'visited':>7s} {'recorded_leg':>12s} "
          f"{'ideal_detour':>12s} {'recorded_total':>14s}")
    for scenario, _, seen, ideal, total, visited in columns:
        print(f'{scenario:16s} {visited:7d} {seen:12.0f} {ideal:12.0f} {total:14.0f}')
    print(f"{'MEAN':16s} {sum(c[5] for c in columns) / n:7.1f} "
          f"{sum(c[2] for c in columns) / n:12.0f} "
          f"{sum(c[3] for c in columns) / n:12.0f} "
          f"{sum(c[4] for c in columns) / n:14.0f}")
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
