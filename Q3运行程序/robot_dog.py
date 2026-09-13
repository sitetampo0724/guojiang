from __future__ import annotations
import argparse
import http.client
import json
import math
import sys
import threading
import time
import unicodedata
import urllib.error
import urllib.request
import uuid
from array import array
from collections import defaultdict
from dataclasses import dataclass
from itertools import combinations
from math import atan2, cos, hypot, inf, isfinite, pi, radians, sin
from pathlib import Path
Point = tuple[float, float]

@dataclass
class Region:
    status: str
    vertices: list[Point]
    witness: Point | None

@dataclass
class Circle:
    center: Point
    radius: float
    support: tuple[Point, ...]

def dot(a, b):
    return a[0] * b[0] + a[1] * b[1]

def cross(a, b):
    return a[0] * b[1] - a[1] * b[0]

def distance(a, b):
    return hypot(a[0] - b[0], a[1] - b[1])

def bearing_halfplanes(stations, bearings, half_width_deg=1.0):
    if len(stations) != len(bearings):
        raise ValueError('One bearing is required per station')
    if not 0 < half_width_deg < 90:
        raise ValueError('Half width must be between 0 and 90 degrees')
    A, b = ([], [])
    for s, theta in zip(stations, bearings):
        if len(s) != 2 or not all((isfinite(v) for v in (*s, theta))):
            raise ValueError('Stations and bearings must be finite')
        lo, hi = (radians(theta - half_width_deg), radians(theta + half_width_deg))
        for normal in [(sin(lo), -cos(lo)), (-sin(hi), cos(hi))]:
            A.append(normal)
            b.append(dot(normal, s))
    return (A, b)

def disk_outer_halfplanes(center, radius, sides=180):
    if radius < 0 or sides < 3:
        raise ValueError('Require radius >= 0 and at least 3 sides')
    A = [(cos(2 * pi * k / sides), sin(2 * pi * k / sides)) for k in range(sides)]
    return (A, [dot(a, center) + radius for a in A])

def convex_hull(points, tol=1e-09):
    pts = []
    for p in sorted(points):
        if not pts or distance(p, pts[-1]) > tol:
            pts.append(p)
    if len(pts) <= 2:
        return pts

    def nonleft(a, b, c):
        u, v = ((b[0] - a[0], b[1] - a[1]), (c[0] - b[0], c[1] - b[1]))
        return cross(u, v) <= tol * max(1.0, hypot(*u), hypot(*v))
    lower, upper = ([], [])
    for p in pts:
        while len(lower) >= 2 and nonleft(lower[-2], lower[-1], p):
            lower.pop()
        lower.append(p)
    for p in reversed(pts):
        while len(upper) >= 2 and nonleft(upper[-2], upper[-1], p):
            upper.pop()
        upper.append(p)
    return lower[:-1] + upper[:-1]

def halfplane_intersection(A, b, tol=1e-07):
    if len(A) != len(b):
        raise ValueError('A and b must have equal lengths')
    rows, rhs = ([], [])
    for a, bi in zip(A, b):
        if len(a) != 2 or not all((isfinite(v) for v in (*a, bi))):
            raise ValueError('Half-plane coefficients must be finite')
        length = hypot(*a)
        if length == 0:
            if bi < -tol:
                return Region('empty', [], None)
            continue
        rows.append((a[0] / length, a[1] / length))
        rhs.append(bi / length)

    def feasible(p):
        return all((dot(a, p) <= bi + tol for a, bi in zip(rows, rhs)))
    candidates = [(0.0, 0.0)]
    candidates.extend(((a[0] * bi, a[1] * bi) for a, bi in zip(rows, rhs)))
    for i, j in combinations(range(len(rows)), 2):
        ai, aj = (rows[i], rows[j])
        det = cross(ai, aj)
        if abs(det) > 1e-14:
            candidates.append(((rhs[i] * aj[1] - ai[1] * rhs[j]) / det, (ai[0] * rhs[j] - rhs[i] * aj[0]) / det))
    valid = [p for p in candidates if feasible(p)]
    if not valid:
        return Region('empty', [], None)
    witness = min(valid, key=lambda p: dot(p, p))
    rays = [(sign * -a[1], sign * a[0]) for a in rows for sign in (-1.0, 1.0)]
    unbounded = not rows or any((all((dot(a, d) <= 1e-12 for a in rows)) for d in rays))
    if unbounded:
        return Region('unbounded', [], witness)
    vertices = convex_hull(valid, tol)
    status = {1: 'point', 2: 'segment'}.get(len(vertices), 'polygon')
    return Region(status, vertices, witness)

def localization_region(stations, bearings, half_width_deg=1.0, extra_A=(), extra_b=(), tol=1e-07):
    A, b = bearing_halfplanes(stations, bearings, half_width_deg)
    return halfplane_intersection(A + list(extra_A), b + list(extra_b), tol)

def diameter(vertices):
    if not vertices:
        raise ValueError('Diameter requires a nonempty bounded region')
    pair = (vertices[0], vertices[0])
    best = 0.0
    for p, q in combinations(vertices, 2):
        d = distance(p, q)
        if d > best:
            best, pair = (d, (p, q))
    return (best, pair)

def _circumcircle(a, b, c):
    u, v = ((b[0] - a[0], b[1] - a[1]), (c[0] - a[0], c[1] - a[1]))
    det = 2 * cross(u, v)
    if abs(det) <= 1e-13 * max(1.0, hypot(*u) * hypot(*v)):
        return None
    uu, vv = (dot(u, u), dot(v, v))
    offset = ((uu * v[1] - u[1] * vv) / det, (u[0] * vv - uu * v[0]) / det)
    return Circle((a[0] + offset[0], a[1] + offset[1]), hypot(*offset), (a, b, c))

def minimum_enclosing_circle(vertices, tol=1e-07):
    pts = convex_hull(vertices, tol)
    if not pts:
        raise ValueError('Circle requires a nonempty bounded region')
    if len(pts) == 1:
        return Circle(pts[0], 0.0, (pts[0],))
    best = None

    def consider(circle):
        nonlocal best
        if circle is None or (best is not None and circle.radius > best.radius + tol):
            return
        radius = max((distance(circle.center, p) for p in pts))
        if radius <= circle.radius + tol:
            circle.radius = max(circle.radius, radius)
            if best is None or circle.radius < best.radius:
                best = circle
    for a, b in combinations(pts, 2):
        consider(Circle(((a[0] + b[0]) / 2, (a[1] + b[1]) / 2), distance(a, b) / 2, (a, b)))
    for a, b, c in combinations(pts, 3):
        consider(_circumcircle(a, b, c))
    if best is None:
        raise ArithmeticError('No enclosing circle found; check numerical scale')
    return best

def region_metrics(region):
    if region.status == 'empty':
        return {'status': 'empty', 'diameter': None, 'radius': None}
    if region.status == 'unbounded':
        return {'status': 'unbounded', 'diameter': inf, 'radius': inf}
    d, endpoints = diameter(region.vertices)
    circle = minimum_enclosing_circle(region.vertices)
    return {'status': region.status, 'vertices': region.vertices, 'diameter': d, 'diameter_endpoints': endpoints, 'center': circle.center, 'radius': circle.radius, 'guaranteed_clearance': circle.radius <= 20.0}

def triangle_counterexample():
    h = 20 * 3 ** 0.5
    vertices = [(0.0, 0.0), (40.0, 0.0), (20.0, h)]
    stations = []
    bearings = []
    for i in range(3):
        a, b = (vertices[i], vertices[(i + 1) % 3])
        unit = ((b[0] - a[0]) / 40, (b[1] - a[1]) / 40)
        stations.append((a[0] - 1000 * unit[0], a[1] - 1000 * unit[1]))
        bearings.append((atan2(unit[1], unit[0]) * 180 / pi + 1) % 360)
    region = localization_region(stations, bearings)
    return {'stations': stations, 'bearings': bearings, **region_metrics(region)}

class ClientError(RuntimeError):
    pass

class ProtocolError(ClientError):
    pass

class ActionRejected(ClientError):

    def __init__(self, status, response):
        self.status, self.response = (status, response)
        super().__init__(f'Simulator rejected action: HTTP {status}, response={response!r}')

class DeadlineExceeded(ClientError):
    pass

class TransportError(ClientError):
    pass

def valid_identifier(value, maximum):
    return isinstance(value, str) and 0 < len(value.encode('utf-8')) <= maximum and (not any((unicodedata.category(c) in ('Cc', 'Cf') for c in value)))

def finite_number(value):
    return isinstance(value, (int, float)) and (not isinstance(value, bool)) and math.isfinite(value)

class HTTPTransport:

    def __init__(self, base_url='http://127.0.0.1:2026', timeout_s=5):
        self.base_url = base_url.rstrip('/')
        self.timeout_s = timeout_s

    def post(self, path, payload):
        body = json.dumps(payload, ensure_ascii=False, allow_nan=False, sort_keys=True, separators=(',', ':')).encode('utf-8')
        request = urllib.request.Request(self.base_url + path, data=body, headers={'Content-Type': 'application/json; charset=utf-8'}, method='POST')
        try:
            with urllib.request.urlopen(request, timeout=self.timeout_s) as response:
                status, raw = (response.status, response.read())
        except urllib.error.HTTPError as error:
            status, raw = (error.code, error.read())
        try:
            result = json.loads(raw.decode('utf-8'))
        except (ValueError, UnicodeError) as error:
            raise ProtocolError('Response is not a valid UTF-8 JSON object') from error
        return (status, result)

class HTTPClient:

    def __init__(self, robot_id, base_url='http://127.0.0.1:2026', log_path=None, transport=None, exit_margin_s=5, timeout_s=5, max_retries=3, retry_delay_s=0.1, clock=time.monotonic, sleeper=time.sleep):
        if not valid_identifier(robot_id, 64):
            raise ValueError('robot_id must be 1..64 UTF-8 bytes without control/format characters')
        if not finite_number(exit_margin_s) or exit_margin_s < 0:
            raise ValueError('exit_margin_s must be finite and nonnegative')
        if not finite_number(timeout_s) or timeout_s <= 0:
            raise ValueError('timeout_s must be positive and finite')
        if isinstance(max_retries, bool) or not isinstance(max_retries, int) or max_retries < 0:
            raise ValueError('max_retries must be a nonnegative integer')
        if not finite_number(retry_delay_s) or retry_delay_s < 0:
            raise ValueError('retry_delay_s must be finite and nonnegative')
        self.robot_id = robot_id
        self.transport = transport if transport is not None else HTTPTransport(base_url, timeout_s)
        self.exit_margin_s, self.timeout_s = (float(exit_margin_s), float(timeout_s))
        self.max_retries, self.retry_delay_s = (max_retries, retry_delay_s)
        self._clock, self._sleep = (clock, sleeper)
        self._lock = threading.Lock()
        self._session = uuid.uuid4().hex
        self._counter = 0
        self._state = 'new'
        self._pending = None
        self.deadline = None
        self.last_virtual_time_s = 0.0
        self.max_virtual_duration_s = 360000.0
        self.log_path = Path(log_path) if log_path is not None else Path('q3_logs') / f'{self._session}.jsonl'
        self.log_path.parent.mkdir(parents=True, exist_ok=True)

    @property
    def remaining_real_s(self):
        return math.inf if self.deadline is None else max(0.0, self.deadline - self._clock())

    @property
    def should_exit(self):
        return self._state == 'active' and (self.remaining_real_s <= self.exit_margin_s or self.last_virtual_time_s >= self.max_virtual_duration_s)

    @property
    def has_pending_action(self):
        return self._pending is not None

    def _log(self, **record):
        record['local_timestamp_ms'] = int(time.time() * 1000)
        with self.log_path.open('a', encoding='utf-8') as stream:
            stream.write(json.dumps(record, ensure_ascii=False, allow_nan=False) + '\n')

    def _check_deadline(self, path):
        if self.deadline is None:
            return
        remaining = self.deadline - self._clock()
        if remaining <= (0 if path == '/exit' else self.exit_margin_s):
            raise DeadlineExceeded('Reality deadline/exit reserve reached; do not issue new search actions')
        if self.last_virtual_time_s >= self.max_virtual_duration_s:
            raise DeadlineExceeded('Virtual time limit reached; simulator may already have closed its interface')

    def _validate_response(self, path, response):
        if not isinstance(response, dict) or type(response.get('accepted')) is not bool:
            raise ProtocolError('Response lacks boolean accepted')
        if not finite_number(response.get('real_timestamp_ms')):
            raise ProtocolError('Response lacks finite real_timestamp_ms')
        if not finite_number(response.get('virtual_time_s')) or response['virtual_time_s'] < 0:
            raise ProtocolError('Response lacks nonnegative virtual_time_s')
        if not response['accepted']:
            return
        if path == '/enter':
            for key in ('remaining_real_duration_s', 'max_real_duration_s', 'max_virtual_duration_s'):
                if not finite_number(response.get(key)) or response[key] < 0:
                    raise ProtocolError(f'Invalid enter field: {key}')
            if response['remaining_real_duration_s'] > 1200:
                raise ProtocolError('remaining_real_duration_s exceeds documented 1200 seconds')
        elif path == '/measure':
            kind = response.get('measure_result')
            if kind not in ('direction', 'near', 'no_signal'):
                raise ProtocolError('Unknown measure_result')
            if kind == 'direction' and (not (finite_number(response.get('svd_deg')) and 0 <= response['svd_deg'] < 360)):
                raise ProtocolError('Invalid svd_deg')
        elif path == '/clear' and response.get('clear_result') not in ('success', 'no_target_in_range'):
            raise ProtocolError('Unknown clear_result')
        elif path == '/exit' and response.get('exit_reason') != 'user_exit':
            raise ProtocolError('Unknown exit_reason')

    def _dispatch_pending(self):
        path, serialized, started = self._pending
        for attempt in range(self.max_retries + 1):
            self._check_deadline(path)
            if isinstance(self.transport, HTTPTransport) and self.deadline is not None:
                reserve = 0 if path == '/exit' else self.exit_margin_s
                self.transport.timeout_s = max(0.001, min(self.timeout_s, self.remaining_real_s - reserve))
            payload = json.loads(serialized)
            self._log(event='request', path=path, payload=payload, attempt=attempt)
            try:
                answer = self.transport.post(path, payload)
                status, response = (200, answer) if isinstance(answer, dict) else answer
                self._log(event='response', path=path, request_id=payload['request_id'], http_status=status, response=response, attempt=attempt)
                self._validate_response(path, response)
                if status != 200 or response['accepted'] is not True:
                    self._pending = None
                    raise ActionRejected(status, response)
            except (OSError, urllib.error.URLError, http.client.HTTPException, ProtocolError) as error:
                self._log(event='transport_error', path=path, request_id=payload['request_id'], attempt=attempt, error=repr(error))
                if attempt == self.max_retries:
                    raise TransportError('Action acknowledgement uncertain; retry_pending() preserves its ID and body') from error
                pause = self.retry_delay_s * 2 ** attempt
                if self.deadline is not None:
                    pause = min(pause, max(0, self.remaining_real_s - (0 if path == '/exit' else self.exit_margin_s)))
                self._sleep(pause)
                continue
            self.last_virtual_time_s = float(response['virtual_time_s'])
            if path == '/enter':
                self.deadline = started + float(response['remaining_real_duration_s'])
                self.max_virtual_duration_s = float(response['max_virtual_duration_s'])
                self._state = 'active'
            elif path == '/exit':
                self._state = 'exited'
            self._pending = None
            return response
        raise AssertionError('unreachable')

    def _action(self, path, point=None, channel=None):
        with self._lock:
            if self._pending is not None:
                raise TransportError('An earlier action is uncertain; call retry_pending() before any new action')
            if path == '/enter' and self._state != 'new':
                raise ClientError('enter is allowed once per client')
            if path != '/enter' and self._state != 'active':
                raise ClientError('Client is not in an active session')
            self._check_deadline(path)
            payload = {'arena_id': 'default', 'robot_id': self.robot_id}
            if point is not None:
                if isinstance(point, dict):
                    if set(point) != {'x', 'y'}:
                        raise ValueError('point mapping must contain exactly x and y')
                    x, y = (point['x'], point['y'])
                else:
                    x, y = point
                if not all((finite_number(v) and abs(v) <= 2000000 for v in (x, y))):
                    raise ValueError('Coordinates must be finite and bounded by 2000000')
                if not finite_number(channel) or not float(channel).is_integer() or (not 1 <= channel <= 20):
                    raise ValueError('channel must be an integer in 1..20')
                payload.update(position={'x': float(x), 'y': float(y)}, channel=int(channel))
            self._counter += 1
            payload['request_id'] = f'{self._session}-{self._counter}'
            serialized = json.dumps(payload, ensure_ascii=False, allow_nan=False, sort_keys=True, separators=(',', ':'))
            self._pending = (path, serialized, self._clock())
            return self._dispatch_pending()

    def retry_pending(self):
        with self._lock:
            if self._pending is None:
                raise ClientError('No pending action')
            return self._dispatch_pending()

    def enter(self):
        return self._action('/enter')

    def measure(self, point, channel):
        return self._action('/measure', point, channel)

    def clear(self, point, channel):
        return self._action('/clear', point, channel)

    def exit(self):
        return self._action('/exit')

class Controller:

    def __init__(self, client, progress=None):
        self.client = client
        self.position = (0.0, 0.0)
        self.channel = 1
        self.observations = defaultdict(list)
        self.no_signals = defaultdict(list)
        self.cleared = set()
        self.scanned = set()
        self.phase = 'search'
        self.progress = progress
        self.virtual_time = 0.0
        self.travel_m = 0.0
        self.measure_count = 0
        self.switch_count = 0
        self.clear_attempts = 0
        self.action_count = 0
        self.phase_time = defaultdict(float)
        self.max_clock_error = 0.0
        self.version = defaultdict(int)
        self.summary_cache = {}
        self.area_cover = None
        self.decisions = []

    def _account(self, point, response, fixed_cost):
        distance = math.dist(self.position, point)
        increment = distance / 5 + fixed_cost
        self.travel_m += distance
        self.virtual_time += increment
        self.position = tuple(point)
        self.action_count += 1
        self.phase_time[self.phase] += increment
        actual = float(response['virtual_time_s'])
        error = abs(actual - self.virtual_time)
        self.max_clock_error = max(self.max_clock_error, error)
        if error > 0.001 + self.action_count * 2e-06:
            raise RuntimeError(f'Virtual clock disagreement: predicted={self.virtual_time}, server={actual}')

    def measure(self, point, channel):
        response = self.client.measure(tuple(point), channel)
        switched = int(channel != self.channel)
        self._account(point, response, 5 + switched)
        self.measure_count += 1
        self.switch_count += switched
        self.channel = channel
        if response['measure_result'] == 'direction':
            reading = float(response['svd_deg'])
            if not any((math.dist(p, point) < 1e-08 for p, _ in self.observations[channel])):
                self.observations[channel].append((tuple(point), reading))
                self.version[channel] += 1
        elif response['measure_result'] == 'no_signal':
            if not any((math.dist(p, point) < 1e-08 for p in self.no_signals[channel])):
                self.no_signals[channel].append(tuple(point))
                self.version[channel] += 1
        return response

    def clear(self, point, channel):
        response = self.client.clear(tuple(point), channel)
        success = response['clear_result'] == 'success'
        self._account(point, response, 3 + 2 * success)
        self.clear_attempts += 1
        if success:
            self.cleared.add(channel)
            if self.progress:
                self.progress(f'cleared channel={channel}; total={len(self.cleared)}; virtual={self.virtual_time:.2f}s')
        return success

    def is_cleared(self, channel):
        return channel in self.cleared

    def pending(self):
        return [c for c in self.observations if c not in self.cleared]

    def channels_to_scan(self):
        return sorted((c for c in range(1, 21) if c not in self.cleared), key=lambda c: (c != self.channel, c))

    def certified(self):
        covered = self.area_cover is not None and self.area_cover.complete
        return len(self.cleared) == 16 or (10 <= len(self.cleared) < 16 and covered and (not self.pending()))

    def result(self):
        return dict(cleared_count=len(self.cleared), cleared_channels=sorted(self.cleared), complete_certificate=self.certified(), scanned_stations=sorted(self.scanned), area_coverage_complete=self.area_cover.complete if self.area_cover is not None else None, census_points=self.area_cover.points if self.area_cover is not None else None, unresolved_coverage_cells=len(self.area_cover.cells) if self.area_cover is not None else None, pending_channels=self.pending(), virtual_time_s=self.virtual_time, average_time_per_cleared_s=self.virtual_time / len(self.cleared) if self.cleared else None, travel_m=self.travel_m, travel_time_s=self.travel_m / 5, measure_count=self.measure_count, measurement_time_s=5 * self.measure_count, switch_count=self.switch_count, switching_time_s=self.switch_count, clear_attempts=self.clear_attempts, clearance_time_s=3 * self.clear_attempts + 2 * len(self.cleared), action_count=self.action_count, phase_time_s=dict(self.phase_time), max_clock_error_s=self.max_clock_error)
ANGLE_ERROR_DEG = 1.005
MAX_INITIAL_RANGE = 1500.0
MIN_RECEPTION_RANGE = 1000.0
CLEAR_RADIUS = 19.9

class LocalizationError(RuntimeError):
    pass

def _validate_observation(position, bearing):
    if len(position) != 2 or not all((isfinite(float(x)) for x in (*position, bearing))):
        raise ValueError('Observation coordinates and bearing must be finite')
    return ((float(position[0]), float(position[1])), float(bearing))

def _append_observation(A, b, position, bearing, distance_bound):
    position, bearing = _validate_observation(position, bearing)
    if not isfinite(distance_bound) or distance_bound <= 0:
        raise ValueError('Distance bound must be positive and finite')
    ai, bi = bearing_halfplanes([position], [bearing], ANGLE_ERROR_DEG)
    theta = radians(bearing)
    unit = (cos(theta), sin(theta))
    A.extend(ai + [unit])
    b.extend(bi + [dot(unit, position) + distance_bound])

def _summary(A, b):
    region = halfplane_intersection(A, b)
    if region.status in ('empty', 'unbounded'):
        raise LocalizationError('Bearing posterior is ' + region.status)
    circle = minimum_enclosing_circle(region.vertices)
    return (circle.center, circle.radius, region.vertices)

def coverage_points(radius=1150.0, rotation_deg=0.0):
    return [(0.0, 0.0)] + [(radius * math.cos(math.radians(rotation_deg) + k * math.pi / 3), radius * math.sin(math.radians(rotation_deg) + k * math.pi / 3)) for k in range(6)]

def covering_radius(radius=1150.0, target_radius=1800.0):
    if radius < 0 or target_radius <= 0:
        raise ValueError('Require radius >= 0 and target_radius > 0')
    if radius >= math.sqrt(3) * target_radius:
        return target_radius
    return max(radius / math.sqrt(3), math.sqrt(target_radius ** 2 + radius ** 2 - math.sqrt(3) * target_radius * radius))

def minimum_feasible_ring_radius(target_radius=1800.0, reception_radius=1000.0):
    if target_radius <= 0 or reception_radius <= 0:
        raise ValueError('Radii must be positive')
    if reception_radius >= target_radius:
        return 0.0
    if reception_radius < target_radius / 2:
        raise ValueError('Seven-site regular arrangement cannot achieve this radius')
    return math.sqrt(3) * target_radius / 2 - math.sqrt(reception_radius ** 2 - target_radius ** 2 / 4)

def optimal_open_route(points, start_index=0):
    n = len(points)
    if not 0 <= start_index < n:
        raise ValueError('Invalid start index')
    if n > 18:
        raise ValueError('Subset DP is intended for at most 18 sites')
    dp = {(1 << start_index, start_index): (0.0, None)}
    for mask in range(1 << n):
        if not mask & 1 << start_index:
            continue
        for last in range(n):
            entry = dp.get((mask, last))
            if entry is None:
                continue
            for nxt in range(n):
                if mask & 1 << nxt:
                    continue
                key = (mask | 1 << nxt, nxt)
                value = entry[0] + math.dist(points[last], points[nxt])
                if key not in dp or value < dp[key][0]:
                    dp[key] = (value, last)
    mask = (1 << n) - 1
    end = min(range(n), key=lambda j: dp.get((mask, j), (math.inf, None))[0])
    length = dp[mask, end][0]
    order = []
    while end is not None:
        order.append(end)
        prev = dp[mask, end][1]
        mask ^= 1 << end
        end = prev
    return (length, list(reversed(order)))

def run():
    radius = 1150.0
    points = coverage_points(radius)
    length, order = optimal_open_route(points)
    assert abs(length - 6 * radius) < 1e-07
    rmin = minimum_feasible_ring_radius()
    assert abs(covering_radius(rmin) - 1000) < 1e-07
    sampled = max((min((math.dist((rho * math.cos(theta), rho * math.sin(theta)), p) for p in points)) for rho in [1800 * j / 120 for j in range(121)] for theta in [2 * math.pi * j / 720 for j in range(720)]))
    assert sampled <= covering_radius(radius) + 1e-07
    print(json.dumps(dict(selected_radius_m=radius, minimum_feasible_radius_m=rmin, points=points, exact_covering_radius_m=covering_radius(radius), reception_margin_m=1000 - covering_radius(radius), independently_sampled_covering_radius_m=sampled, optimal_open_route_length_m=length, optimal_open_route_indices=order, detection_count=140, minimum_channel_switch_count=133, complete_census_time_s=length / 5 + 140 * 5 + 133, caveat='Census time excludes localization and clearing detours; not simulator results.'), indent=2))
DEFAULT_HORIZON = 12

def warm_route(start, centers):
    ids = sorted(centers)
    locations = [start] + [centers[c] for c in ids]
    n = len(ids)
    if not n:
        return []
    d = [[math.dist(a, b) for b in locations] for a in locations]
    best_length, best = (math.inf, None)
    for first in range(1, n + 1):
        todo = set(range(1, n + 1)) - {first}
        route = [0, first]
        while todo:
            nxt = min(todo, key=lambda k: (d[route[-1]][k], k))
            todo.remove(nxt)
            route.append(nxt)
        for _ in range(50):
            improvement, swap = (-1e-07, None)
            for i in range(1, n):
                for j in range(i + 1, n + 1):
                    delta = d[route[i - 1]][route[j]] - d[route[i - 1]][route[i]]
                    if j < n:
                        delta += d[route[i]][route[j + 1]] - d[route[j]][route[j + 1]]
                    if delta < improvement:
                        improvement, swap = (delta, (i, j))
            if swap is None:
                break
            i, j = swap
            route[i:j + 1] = reversed(route[i:j + 1])
        length = sum((d[a][z] for a, z in zip(route, route[1:])))
        if length < best_length:
            best_length, best = (length, route)
    return [ids[k - 1] for k in best[1:]]

def route_distance(start, route, centers):
    return sum((math.dist(a, b) for a, b in zip([start] + [centers[c] for c in route], [centers[c] for c in route])))

def solve_route(start, centers, horizon=DEFAULT_HORIZON):
    if not isinstance(horizon, int) or not 1 <= horizon <= 16:
        raise ValueError('DP horizon must be an integer in 1..16')
    if not all((math.isfinite(v) for p in [start] + list(centers.values()) for v in p)):
        raise ValueError('DP coordinates must be finite')
    warm = warm_route(start, centers)
    active, suffix = (warm[:horizon], warm[horizon:])
    n = len(active)
    warm_cost = route_distance(start, warm, centers)
    if not n:
        return ([], dict(horizon=0, task_count=0, states=0, exact_for_frozen_tasks=True, route_distance_m=0.0, warm_distance_m=0.0))
    points = [centers[c] for c in active]
    distances = [[math.dist(a, b) for b in points] for a in points]
    full = (1 << n) - 1
    dp = array('d', [math.inf]) * ((full + 1) * n)
    parent = array('b', [-1]) * ((full + 1) * n)
    states = 0
    for j, p in enumerate(points):
        dp[(1 << j) * n + j] = math.dist(start, p)
    for mask in range(1, full + 1):
        remaining = full ^ mask
        last_bits = mask
        while last_bits:
            last_bit = last_bits & -last_bits
            j = last_bit.bit_length() - 1
            last_bits ^= last_bit
            value = dp[mask * n + j]
            states += 1
            todo = remaining
            while todo:
                bit = todo & -todo
                k = bit.bit_length() - 1
                todo ^= bit
                index = (mask | bit) * n + k
                candidate = value + distances[j][k]
                if candidate < dp[index]:
                    dp[index] = candidate
                    parent[index] = j
    tail_first = centers[suffix[0]] if suffix else None
    last = min(range(n), key=lambda j: dp[full * n + j] + (math.dist(points[j], tail_first) if tail_first else 0.0))
    order, mask = ([], full)
    while last >= 0:
        order.append(active[last])
        previous = parent[mask * n + last]
        mask ^= 1 << last
        last = previous
    route = list(reversed(order)) + suffix
    cost = route_distance(start, route, centers)
    if cost >= warm_cost - 1e-07:
        route, cost = (warm, warm_cost)
    return (route, dict(horizon=n, task_count=len(centers), states=states, exact_for_frozen_tasks=not suffix, route_distance_m=cost, warm_distance_m=warm_cost))
DISCOVERY_GAIN_M2 = 2400000.0

class AreaCover:

    def __init__(self, max_depth=10):
        self.cells = [(0.0, 0.0, 1800.0, 0)]
        self.points = []
        self.max_depth = max_depth

    @staticmethod
    def outside(x, y, h):
        return math.hypot(max(0.0, abs(x) - h), max(0.0, abs(y) - h)) > 1800.0 + 1e-07

    @staticmethod
    def inside(cell, point):
        x, y, h, _ = cell
        return math.hypot(abs(x - point[0]) + h, abs(y - point[1]) + h) <= 999.9

    def add(self, point):
        if any((math.dist(point, p) < 1e-07 for p in self.points)):
            return
        self.points.append(tuple(point))
        stack = list(self.cells)
        remaining = []
        while stack:
            cell = stack.pop()
            x, y, h, depth = cell
            if self.outside(x, y, h) or any((self.inside(cell, p) for p in self.points)):
                continue
            intersects = any((math.hypot(max(0.0, abs(x - p[0]) - h), max(0.0, abs(y - p[1]) - h)) <= 999.9 for p in self.points))
            if depth < self.max_depth and intersects:
                hh = h / 2
                stack.extend(((x + sx * hh, y + sy * hh, hh, depth + 1) for sx in (-1, 1) for sy in (-1, 1)))
            else:
                remaining.append(cell)
        self.cells = remaining

    @property
    def complete(self):
        return not self.cells

    def gain(self, point):
        return sum((4 * h * h for x, y, h, _ in self.cells if math.dist((x, y), point) <= 999.9))

def clip_linear(poly, normal, limit):
    output = []
    limit += 1e-07
    for a, z in zip(poly, poly[1:] + poly[:1]):
        fa = a[0] * normal[0] + a[1] * normal[1] - limit
        fz = z[0] * normal[0] + z[1] * normal[1] - limit
        if fa <= 0:
            output.append(a)
        if fa < 0 < fz or fz < 0 < fa:
            t = fa / (fa - fz)
            output.append((a[0] + t * (z[0] - a[0]), a[1] + t * (z[1] - a[1])))
    return output

def exclude_disk_hull(poly, point, radius=999.9):
    kept = [p for p in poly if math.dist(p, point) >= radius - 1e-07]
    for a, z in zip(poly, poly[1:] + poly[:1]):
        u = (z[0] - a[0], z[1] - a[1])
        v = (a[0] - point[0], a[1] - point[1])
        aa = u[0] ** 2 + u[1] ** 2
        if aa < 1e-16:
            continue
        bb = 2 * (u[0] * v[0] + u[1] * v[1])
        cc = v[0] ** 2 + v[1] ** 2 - radius ** 2
        disc = bb * bb - 4 * aa * cc
        if disc >= 0:
            for t in ((-bb - math.sqrt(disc)) / (2 * aa), (-bb + math.sqrt(disc)) / (2 * aa)):
                if 0 <= t <= 1:
                    kept.append((a[0] + t * u[0], a[1] + t * u[1]))
    return convex_hull(kept)

def posterior(ctrl, channel, A=None, b=None):
    if A is None:
        A, b = ([], [])
        for point, bearing in ctrl.observations[channel]:
            _append_observation(A, b, point, bearing, 1500.0)
    _, _, poly = _summary(A, b)
    for k in range(96):
        t = k * 2 * math.pi / 96
        poly = clip_linear(poly, (math.cos(t), math.sin(t)), 1800.0)
    for point in ctrl.no_signals[channel]:
        poly = exclude_disk_hull(poly, point)
        if not poly:
            raise LocalizationError('No-signal history contradicted bearing posterior')
    if not poly:
        raise LocalizationError('Target domain contradicted bearing posterior')
    circle = minimum_enclosing_circle(poly)
    return dict(center=circle.center, radius=circle.radius, polygon=poly)

def adaptive_summary(ctrl, channel):
    key = ('adaptive', channel)
    cached = ctrl.summary_cache.get(key)
    if cached is None or cached[0] != ctrl.version[channel]:
        cached = (ctrl.version[channel], posterior(ctrl, channel))
        ctrl.summary_cache[key] = cached
    return cached[1]

def clearance_entry(position, center, vertices):
    distance = math.dist(position, center)
    if distance < 1e-08:
        return center
    unit = ((position[0] - center[0]) / distance, (position[1] - center[1]) / distance)
    step = distance
    for v in vertices:
        w = (center[0] - v[0], center[1] - v[1])
        projection = w[0] * unit[0] + w[1] * unit[1]
        limit = -projection + math.sqrt(max(0.0, projection ** 2 + CLEAR_RADIUS ** 2 - w[0] ** 2 - w[1] ** 2))
        step = min(step, max(0.0, limit - 1e-06))
    return (center[0] + step * unit[0], center[1] + step * unit[1])

def scan_unknown(ctrl, cover, point):
    if cover.complete or len(ctrl.cleared) + len(ctrl.pending()) == 16:
        return
    ctrl.phase = 'adaptive_search'
    unknown = [c for c in range(1, 21) if c not in ctrl.cleared and c not in ctrl.observations]
    for c in sorted(unknown, key=lambda c: (c != ctrl.channel, c)):
        if len(ctrl.cleared) + len(ctrl.pending()) == 16:
            break
        response = ctrl.measure(point, c)
        if response['measure_result'] == 'near':
            if not ctrl.clear(point, c):
                raise LocalizationError('Near observation contradicted by clearance')
    cover.add(point)

def shared_bearings(ctrl, point):
    for c in sorted(ctrl.pending(), key=lambda c: (c != ctrl.channel, c)):
        obs = ctrl.observations[c]
        if len(obs) >= 3 or any((math.dist(point, p) < 40 for p, _ in obs)):
            continue
        info = adaptive_summary(ctrl, c)
        if info['radius'] <= CLEAR_RADIUS or math.dist(point, info['center']) > 1200:
            continue
        old = obs[-1][0]
        center = info['center']
        u = (center[0] - old[0], center[1] - old[1])
        v = (center[0] - point[0], center[1] - point[1])
        sine = abs(u[0] * v[1] - u[1] * v[0]) / max(1.0, math.hypot(*u) * math.hypot(*v))
        if sine < 0.22:
            continue
        ctrl.phase = 'shared_bearings'
        response = ctrl.measure(point, c)
        if response['measure_result'] == 'near' and (not ctrl.clear(point, c)):
            raise LocalizationError('Near shared observation contradicted by clearance')

def fast_service(ctrl, channel, on_stop=None):
    A, b = ([], [])
    for point, bearing in ctrl.observations[channel]:
        _append_observation(A, b, point, bearing, 1500.0)
    q = 1.25 / (2 * math.cos(math.radians(ANGLE_ERROR_DEG)) ** 2)
    previous = None
    for _ in range(12):
        info = posterior(ctrl, channel, A, b)
        center, radius, vertices = (info['center'], info['radius'], info['polygon'])
        if previous is not None and radius > q * previous + 1e-05:
            raise LocalizationError('Adaptive posterior did not contract')
        ctrl.phase = 'adaptive_localization'
        if radius <= CLEAR_RADIUS:
            point = clearance_entry(ctrl.position, center, vertices)
            if not ctrl.clear(point, channel):
                raise LocalizationError('Failed geometrically certified clearance')
            if on_stop:
                on_stop(point)
            return
        a, z = max(((a, z) for a in vertices for z in vertices), key=lambda pair: math.dist(*pair))
        length = math.dist(a, z)
        offset = min(100.0, radius * 0.25)
        normal = (-(z[1] - a[1]) / max(length, 1e-12), (z[0] - a[0]) / max(length, 1e-12))
        options = [(center[0] + s * offset * normal[0], center[1] + s * offset * normal[1]) for s in (-1, 1)]
        point = min(options, key=lambda p: math.dist(ctrl.position, p))
        bound = max((math.dist(point, v) for v in vertices))
        if bound >= 1000:
            raise LocalizationError('Adaptive probe is not reception-safe')
        response = ctrl.measure(point, channel)
        if response['measure_result'] == 'near':
            if not ctrl.clear(point, channel):
                raise LocalizationError('Near adaptive observation contradicted')
            if on_stop:
                on_stop(point)
            return
        if response['measure_result'] != 'direction':
            raise LocalizationError('Missing signal at reception-safe adaptive probe')
        _append_observation(A, b, point, response['svd_deg'], bound)
        previous = radius
    raise LocalizationError('Adaptive localization exceeded analytic round bound')

def partition_uncovered(cover, original, unvisited):
    groups = {i: [] for i in unvisited}
    stack = list(cover.cells)
    while stack:
        cell = stack.pop()
        x, y, h, depth = cell
        if cover.outside(x, y, h):
            continue
        candidates = [i for i in unvisited if cover.inside(cell, original[i])]
        if candidates:
            i = min(candidates, key=lambda i: (math.dist((x, y), original[i]), i))
            groups[i].extend(((x + sx * h, y + sy * h) for sx in (-1, 1) for sy in (-1, 1)))
        elif depth < cover.max_depth:
            hh = h / 2
            stack.extend(((x + sx * hh, y + sy * hh, hh, depth + 1) for sx in (-1, 1) for sy in (-1, 1)))
        else:
            raise LocalizationError('Uncovered cell lost its fallback station')
    return {i: convex_hull(points) for i, points in groups.items() if points}

def feasible_step(start, target, witnesses):
    length = math.dist(start, target)
    if length < 1e-08:
        return start
    u = ((target[0] - start[0]) / length, (target[1] - start[1]) / length)
    step = length
    for w in witnesses:
        v = (start[0] - w[0], start[1] - w[1])
        projection = v[0] * u[0] + v[1] * u[1]
        discriminant = projection ** 2 + 999.9 ** 2 - v[0] ** 2 - v[1] ** 2
        limit = -projection + math.sqrt(max(0.0, discriminant))
        step = min(step, max(0.0, limit - 1e-05))
    p = (start[0] + step * u[0], start[1] + step * u[1])
    return p if all((math.dist(p, w) <= 999.9 for w in witnesses)) else start
ARENA_RADIUS = 1800.0
SCAN_RADIUS = 999.9
PLANNING_RADIUS = 975.0
EDGE_RADIUS = 800.0
AREA_FRACTION = 0.95

def project_arena(point):
    r = math.hypot(*point)
    return tuple(point) if r <= ARENA_RADIUS else tuple((v * ARENA_RADIUS / r for v in point))

class CirclePlanner:

    def __init__(self, cover):
        self.cover = cover
        self.samples = [(float(x), float(y), 10000.0) for x in range(-1650, 1700, 100) for y in range(-1650, 1700, 100) if x * x + y * y <= 1700.0 ** 2]
        weight = math.pi * (1800.0 ** 2 - 1700.0 ** 2) / 180
        self.samples += [(1800 * math.cos(k * math.pi / 90), 1800 * math.sin(k * math.pi / 90), weight) for k in range(180)]
        self.base = [(float(x), float(y)) for x in range(-1800, 1801, 200) for y in range(-1800, 1801, 200) if x * x + y * y <= ARENA_RADIUS ** 2]
        self.cached = None

    def scanned(self, point):
        self.samples = [p for p in self.samples if math.dist(p[:2], point) > SCAN_RADIUS]
        self.cached = None

    def select(self, start):
        if self.cover.complete:
            return (None, dict(method='complete', sampled_new_area_m2=0.0))
        if self.cached is None:
            scored = {}

            def evaluate(point):
                p = tuple((round(v, 5) for v in project_arena(point)))
                if p not in scored:
                    scored[p] = sum((w for x, y, w in self.samples if (x - p[0]) ** 2 + (y - p[1]) ** 2 <= PLANNING_RADIUS ** 2))
            for p in self.base:
                evaluate(p)
            for step in (100.0, 40.0, 10.0):
                seeds = sorted(scored, key=lambda p: (-scored[p], p))[:8]
                for p in seeds:
                    for dx in (-step, 0.0, step):
                        for dy in (-step, 0.0, step):
                            evaluate((p[0] + dx, p[1] + dy))
            best = max(scored.values(), default=0)
            if best:
                near_best = [p for p, g in scored.items() if g >= AREA_FRACTION * best]
                self.cached = (near_best, scored, best)
            else:
                self.cached = ([], scored, 0)
        choices, scored, best = self.cached
        if choices:
            point = min(choices, key=lambda p: (math.dist(start, p), -scored[p], p))
            return (point, dict(method='sampled_area_then_nearest', sampled_new_area_m2=scored[point], best_candidate_area_m2=best, retained_area_fraction=AREA_FRACTION, candidates=len(scored)))
        cell = min(self.cover.cells, key=lambda c: math.dist(start, project_arena(c[:2])))
        point = project_arena(cell[:2])
        return (point, dict(method='certified_hole_repair', sampled_new_area_m2=0.0, unresolved_cell=list(cell)))
MIN_PLAN_SAVING_S = 30.0
SHARED_VIEW_MIN_CHANNELS = 2

def copy_cover(cover):
    clone = AreaCover(cover.max_depth)
    clone.points, clone.cells = (list(cover.points), list(cover.cells))
    return clone

def reference_stations(cover, original):
    probe = copy_cover(cover)
    available, needed = (set(range(1, 7)), [])
    while not probe.complete:
        if not available:
            raise LocalizationError('Original D reference cannot complete coverage')
        i = max(available, key=lambda i: (probe.gain(original[i]), -i))
        available.remove(i)
        needed.append(i)
        probe.add(original[i])
    return needed

def insertion(start, route, centers, point):
    best = (math.inf, 0)
    for i in range(len(route) + 1):
        a = start if i == 0 else centers[route[i - 1]]
        b = centers[route[i]] if i < len(route) else None
        extra = math.dist(a, point) + (math.dist(point, b) - math.dist(a, b) if b else 0.0)
        best = min(best, (extra, i))
    return best

def circle_candidates(start, sources, groups, original, area_point=None):
    points, kinds = ({}, {})

    def add(p, kind, task=None):
        p = project_arena(p)
        if any((math.dist(p, q) < 1e-05 for q in points.values())):
            return
        task = task if task is not None else 100 + len(points)
        points[task], kinds[task] = (p, kind)
    for i in groups:
        add(original[i], 'reference', 20 + i)
    add(start, 'scan_at_current_stop')
    for p in sources.values():
        add(p, 'scan_near_known_source')
    if area_point is not None:
        add(area_point, 'large_new_area')
    for i, witnesses in groups.items():
        targets = [start] + sorted(sources.values(), key=lambda p: math.dist(p, original[i]))[:3]
        for target in targets:
            add(feasible_step(original[i], target, witnesses), 'slide_to_route')
        circle = minimum_enclosing_circle(witnesses)
        if circle.radius < 999.8:
            add(circle.center, 'region_center')
            for target in targets:
                add(feasible_step(circle.center, target, witnesses), 'slide_to_route')
    ids = sorted(groups)
    for j, i in enumerate(ids):
        for k in ids[j + 1:]:
            witnesses = convex_hull(groups[i] + groups[k])
            circle = minimum_enclosing_circle(witnesses)
            if circle.radius < 999.8:
                add(circle.center, 'merge_regions')
                targets = [start] + sorted(sources.values(), key=lambda p: math.dist(p, circle.center))[:2]
                for target in targets:
                    add(feasible_step(circle.center, target, witnesses), 'merge_and_slide')
    masks = {c: sum((1 << j for j, i in enumerate(ids) if all((math.dist(p, w) <= 999.9 for w in groups[i])))) for c, p in points.items()}
    return ({c: p for c, p in points.items() if masks[c]}, {c: m for c, m in masks.items() if m}, kinds)

def covering_routes(start, sources, candidates, masks, full, scan_cost_s, beam_width=3):
    all_centers = dict(sources)
    all_centers.update(candidates)
    initial = warm_route(start, sources)
    initial_cost = route_distance(start, initial, sources) / 5
    beams = {0: [(initial_cost, (), tuple(initial))]}
    for mask in range(full + 1):
        for cost, chosen, route in list(beams.get(mask, [])):
            for c, coverage in masks.items():
                combined = mask | coverage
                if combined == mask:
                    continue
                extra, index = insertion(start, route, all_centers, candidates[c])
                subset = tuple(sorted(chosen + (c,)))
                new_route = route[:index] + (c,) + route[index:]
                entry = (cost + extra / 5 + scan_cost_s, subset, new_route)
                options = beams.setdefault(combined, [])
                previous = next((old for old in options if old[1] == subset), None)
                if previous is not None:
                    if previous[0] <= entry[0]:
                        continue
                    options.remove(previous)
                options.append(entry)
                options.sort()
                del options[beam_width:]
    return ([list(chosen) for _, chosen, _ in beams.get(full, [])], sum((len(b) for b in beams.values())))

def choose_plan_reference(ctrl, cover, original, area_planner):
    infos = {c: adaptive_summary(ctrl, c) for c in ctrl.pending()}
    sources = {c: info['center'] for c, info in infos.items()}
    unknown = [c for c in range(1, 21) if c not in ctrl.cleared and c not in ctrl.observations]
    scan_cost = 6.0 * len(unknown)
    forecast = copy_cover(cover)
    anticipated = []
    if len(ctrl.cleared) + len(sources) < 16:
        for c, p in sources.items():
            if infos[c]['radius'] <= 250 and forecast.gain(p) >= DISCOVERY_GAIN_M2:
                forecast.add(p)
                anticipated.append(c)
    need_search = not forecast.complete and len(ctrl.cleared) + len(sources) < 16
    needed = reference_stations(forecast, original) if need_search else []
    baseline = dict(sources)
    baseline.update({20 + i: original[i] for i in needed})
    if not baseline:
        raise LocalizationError('No tasks remain before a completion certificate')
    horizon = 16 if not needed else 12
    route, stats = solve_route(ctrl.position, baseline, horizon=horizon)
    baseline_cost = stats['route_distance_m'] / 5 + len(needed) * scan_cost
    best = dict(centers=baseline, route=route, dp=stats, estimated_s=baseline_cost, plan_type='reference_d', scan_tasks=[20 + i for i in needed], scan_task_kinds={20 + i: 'reference' for i in needed})
    proposals, beam_states, candidate_count = ([], 0, 0)
    if need_search:
        groups = partition_uncovered(forecast, original, set(range(1, 7)))
        area_point, _ = area_planner.select(ctrl.position)
        candidates, masks, kinds = circle_candidates(ctrl.position, sources, groups, original, area_point)
        candidate_count = len(candidates)
        sets, beam_states = covering_routes(ctrl.position, sources, candidates, masks, (1 << len(groups)) - 1, scan_cost)
        for selected in sets:
            centers = dict(sources)
            centers.update({c: candidates[c] for c in selected})
            proposed, info = solve_route(ctrl.position, centers)
            cost = info['route_distance_m'] / 5 + len(selected) * scan_cost
            proposals.append(dict(estimated_s=cost, scan_count=len(selected), kinds=[kinds[c] for c in selected]))
            if cost < baseline_cost - MIN_PLAN_SAVING_S and cost < best['estimated_s'] - 1e-06:
                best = dict(centers=centers, route=proposed, dp=info, estimated_s=cost, plan_type='joint_circle_route', scan_tasks=selected, scan_task_kinds={c: kinds[c] for c in selected})
    best.update(reference_estimated_s=baseline_cost, alternatives=proposals, beam_states=beam_states, circle_candidate_count=candidate_count, anticipated_scan_sources=anticipated, unknown_channels=unknown, beliefs=infos)
    return best

PLAN_SAMPLE_STEP = 60.0
PLAN_DISK_RADIUS = 999.9
PLAN_SAFETY_MARGIN = 9.0
PLAN_RING_RADII = (1150.0, 1240.0, 1330.0, 1420.0, 1510.0, 1600.0, 1690.0)
PLAN_RING_ANGLES = 18
PLAN_MAX_STATIONS = 12
PLAN_NEW_PLAN_MARGIN_S = 30.0
FUTURE_EXPECTED_TOTAL_SOURCES = 13.5
FUTURE_PRIOR_RADIUS = 1150.0
FUTURE_PENALTY_MIN_RADIUS = 1500.0
FUTURE_PENALTY_RAMP_M = 200.0
FUTURE_SAMPLES = 36
FUTURE_DETOUR_FACTOR = 0.6
PLAN_BOUNDARY_SAMPLES = 720

class StationGrid:
    """Area-uniform sample grid used to score candidate census stations."""

    def __init__(self, step=PLAN_SAMPLE_STEP, arena=ARENA_RADIUS):
        dy = step * math.sqrt(3) / 2
        rows = int(arena / dy) + 1
        points = []
        for j in range(-rows, rows + 1):
            y = j * dy
            if abs(y) > arena:
                continue
            half = int(math.sqrt(max(0.0, arena * arena - y * y)) / step) + 1
            offset = step / 2.0 if j % 2 else 0.0
            for i in range(-half, half + 1):
                x = i * step + offset
                if x * x + y * y <= arena * arena:
                    points.append((x, y))
        for k in range(PLAN_BOUNDARY_SAMPLES):
            angle = 2 * math.pi * k / PLAN_BOUNDARY_SAMPLES
            points.append((arena * math.cos(angle), arena * math.sin(angle)))
        self.points = points
        self.total = len(points)
        self.full = (1 << self.total) - 1
        self.radius = PLAN_DISK_RADIUS - PLAN_SAFETY_MARGIN
        self.cell_area_m2 = step * dy
        self._masks = {}

    @property
    def area_m2(self):
        return self.total * self.cell_area_m2

    def mask(self, point):
        key = (round(point[0], 3), round(point[1], 3))
        cached = self._masks.get(key)
        if cached is None:
            radius2 = self.radius * self.radius
            x0, y0 = point
            value = 0
            bit = 1
            for x, y in self.points:
                if (x - x0) ** 2 + (y - y0) ** 2 <= radius2:
                    value |= bit
                bit <<= 1
            cached = value
            self._masks[key] = cached
        return cached

    @staticmethod
    def uncovered(mask, covered):
        return mask & ~covered

def ring_station_candidates():
    limit = PLAN_DISK_RADIUS - PLAN_SAFETY_MARGIN
    out = []
    for radius in PLAN_RING_RADII:
        if covering_radius(radius) > limit:
            continue
        for k in range(PLAN_RING_ANGLES):
            phase = 2 * math.pi * k / PLAN_RING_ANGLES
            for i in range(6):
                angle = phase + i * math.pi / 3
                out.append((radius * math.cos(angle), radius * math.sin(angle)))
    return out

def cheapest_insertion(route, point):
    best_extra, best_index = (math.inf, 0)
    for i in range(len(route)):
        a = route[i]
        b = route[i + 1] if i + 1 < len(route) else None
        extra = math.dist(a, point) + (math.dist(point, b) - math.dist(a, b) if b else 0.0)
        if extra < best_extra:
            best_extra, best_index = (extra, i + 1)
    return (best_extra, best_index)

def _two_opt(seq):
    """Open-route 2-opt: seq[0] is the fixed start."""
    for _ in range(8):
        changed = False
        for i in range(1, len(seq) - 1):
            a = seq[i - 1]
            for j in range(i + 1, len(seq)):
                b, c = (seq[i], seq[j])
                d = seq[j + 1] if j + 1 < len(seq) else None
                delta = math.dist(a, c) - math.dist(a, b) - (math.dist(c, d) if d else 0.0)
                if d is not None:
                    delta += math.dist(b, d)
                if delta < -1e-09:
                    seq[i:j + 1] = list(reversed(seq[i:j + 1]))
                    changed = True
        if not changed:
            break
    return seq

def _or_opt(seq, rounds=3):
    """Relocate single nodes to their best position (open route, fixed start)."""
    for _ in range(rounds):
        improved = False
        for i in range(1, len(seq)):
            node = seq[i]
            before = seq[i - 1]
            after = seq[i + 1] if i + 1 < len(seq) else None
            saving = math.dist(before, node) + (math.dist(node, after) if after else 0.0)
            if after:
                saving -= math.dist(before, after)
            rest = seq[:i] + seq[i + 1:]
            best_delta, best_index = (0.0, None)
            for j in range(len(rest) + 1):
                a = rest[j - 1] if j else None
                b = rest[j] if j < len(rest) else None
                insertion = (math.dist(a, node) if a else 0.0) + (math.dist(node, b) if b else 0.0)
                if a and b:
                    insertion -= math.dist(a, b)
                if insertion - saving < best_delta - 1e-09:
                    best_delta, best_index = (insertion - saving, j)
            if best_index is not None:
                seq = rest[:best_index] + [node] + rest[best_index:]
                improved = True
        if not improved:
            break
    return seq

def fast_route_length(start, points):
    """Cheap open-route length: nearest neighbour + 2-opt + or-opt."""
    seq = fast_route_order(start, points)
    return sum((math.dist(a, b) for a, b in zip(seq, seq[1:])))

def fast_route_order(start, points):
    """Nearest neighbour + 2-opt + or-opt ordering (start first)."""
    points = list(points)
    count = len(points)
    if not count:
        return [start]
    if count > 60:
        seq = [start] + points
    else:
        remaining = set(range(count))
        order = []
        current = start
        while remaining:
            nxt = min(remaining, key=lambda i: math.dist(current, points[i]))
            order.append(nxt)
            remaining.discard(nxt)
            current = points[nxt]
        seq = [start] + [points[i] for i in order]
    seq = _two_opt(seq)
    seq = _or_opt(seq)
    seq = _two_opt(seq)
    return seq

def ring_phase_candidates(sources):
    """Uniform phases plus phases that park a ring station on a pending source."""
    phases = [2 * math.pi * k / PLAN_RING_ANGLES for k in range(PLAN_RING_ANGLES)]
    for point in sources:
        angle = math.atan2(point[1], point[0]) % (math.pi / 3)
        phases.append(angle)
    unique = []
    for phase in sorted(phases):
        if all((abs(phase - seen) > 0.02 for seen in unique)):
            unique.append(phase)
    while unique and unique[0] + math.pi / 3 - unique[-1] < 0.02:
        unique.pop()
    return unique

def ring_station_options(covered_now, scanned, grid, sources):
    """All (radius, phase) hexagon rings that complete the certificate alone."""
    limit = PLAN_DISK_RADIUS - PLAN_SAFETY_MARGIN
    options = []
    for radius in PLAN_RING_RADII:
        if covering_radius(radius) > limit:
            continue
        for k, phase in enumerate(ring_phase_candidates(sources)):
            ring = []
            union = covered_now
            for i in range(6):
                angle = phase + i * math.pi / 3
                point = (radius * math.cos(angle), radius * math.sin(angle))
                ring.append(point)
                union |= grid.mask(point)
            if union != grid.full:
                continue
            stations = [p for p in ring if all((math.dist(p, q) > 1.0 for q in scanned))]
            options.append(dict(radius=radius, phase=phase, phase_index=k, ring=ring, stations=stations))
    return options

def slide_towards(point, target, fraction):
    return (point[0] + fraction * (target[0] - point[0]), point[1] + fraction * (target[1] - point[1]))

def covers(grid, covered_now, stations):
    union = covered_now
    for point, _ in stations:
        union |= grid.mask(point)
        if union == grid.full:
            return True
    return union == grid.full

def station_plan_cost(start, source_route, stations, scan_cost_s, terms=None):
    ring_points = [point for point, kind in stations if kind != 'at_source']
    route_points = source_route + ring_points
    travel = fast_route_length(start, route_points)
    cost = travel / 5.0 + len(stations) * scan_cost_s
    if terms is not None:
        cost += future_source_penalty(start, route_points, terms)
    return cost

def point_segment_distance(point, a, b):
    dx = b[0] - a[0]
    dy = b[1] - a[1]
    square = dx * dx + dy * dy
    if square < 1e-12:
        return math.dist(point, a)
    t = ((point[0] - a[0]) * dx + (point[1] - a[1]) * dy) / square
    t = max(0.0, min(1.0, t))
    return math.dist(point, (a[0] + t * dx, a[1] + t * dy))

def point_route_distance(point, seq):
    best = min((math.dist(point, node) for node in seq))
    for a, b in zip(seq, seq[1:]):
        value = point_segment_distance(point, a, b)
        if value < best:
            best = value
    return best

def future_source_terms(ctrl, pending_count):
    """Expected number of still-undetected sources and their typical radius."""
    detected = len(ctrl.cleared) + pending_count
    expected = max(0.0, FUTURE_EXPECTED_TOTAL_SOURCES - detected)
    radii = sorted((math.hypot(point[0], point[1]) for point in getattr(ctrl, 'source_positions', {}).values()))
    radius = radii[len(radii) // 2] if radii else None
    if radius is None:
        return (0.0, radius)
    weight = (radius - FUTURE_PENALTY_MIN_RADIUS) / FUTURE_PENALTY_RAMP_M
    weight = max(0.0, min(1.0, weight))
    return (expected * weight, radius)

def future_source_penalty(start, route_points, terms):
    expected, radius = terms
    if expected <= 0.0:
        return 0.0
    radius = FUTURE_PRIOR_RADIUS if radius is None else radius
    seq = fast_route_order(start, route_points)
    total = 0.0
    for k in range(FUTURE_SAMPLES):
        angle = 2 * math.pi * k / FUTURE_SAMPLES
        total += point_route_distance((radius * math.cos(angle), radius * math.sin(angle)), seq)
    return expected * FUTURE_DETOUR_FACTOR * (total / FUTURE_SAMPLES) / 5.0

def refine_station_set(grid, covered_now, stations, sources, start, source_route, scan_cost_s, terms=None, rounds=3):
    """Drop or slide stations while the coverage certificate still holds."""
    best = station_plan_cost(start, source_route, stations, scan_cost_s, terms)
    for _ in range(rounds):
        changed = False
        for i in range(len(stations)):
            point, kind = stations[i]
            if kind == 'at_source':
                alternatives = []
            else:
                alternatives = []
                for target in sorted(sources, key=lambda q: math.dist(point, q))[:3]:
                    if math.dist(point, target) < 5.0:
                        alternatives.append(('at_source', target))
                        continue
                    for fraction in (0.25, 0.5, 0.75, 1.0):
                        alternatives.append((kind, slide_towards(point, target, fraction)))
            for new_kind, alternative in alternatives:
                trial = list(stations)
                trial[i] = (alternative, new_kind)
                if not covers(grid, covered_now, trial):
                    continue
                cost = station_plan_cost(start, source_route, trial, scan_cost_s, terms)
                if cost < best - 1.0:
                    stations, best, changed = (trial, cost, True)
                    break
        if not changed:
            break
    for i in range(len(stations) - 1, -1, -1):
        trial = stations[:i] + stations[i + 1:]
        if not covers(grid, covered_now, trial):
            continue
        cost = station_plan_cost(start, source_route, trial, scan_cost_s, terms)
        if cost < best - 1.0:
            stations, best = (trial, cost)
    return (stations, best)

PLAN_SEARCH_STEPS = 160
PLAN_SEEDS_REFINED = 3
PLAN_OUTER_RADII = (1200.0, 1400.0, 1600.0, 1770.0)

def outer_station_seeds(grid, covered_now, scanned, source_route):
    """Hexagon seeds pushed out towards the arena rim, aligned with sources."""
    phases = [0.0]
    if source_route:
        radii = sorted(math.hypot(p[0], p[1]) for p in source_route)
        if radii[len(radii) // 2] > 900.0:
            for point in source_route:
                phases.append(math.atan2(point[1], point[0]) % (math.pi / 3))
    seeds = []
    for radius in PLAN_OUTER_RADII:
        if covering_radius(radius) > PLAN_DISK_RADIUS - PLAN_SAFETY_MARGIN:
            continue
        for phase in phases[:4]:
            ring = [(radius * math.cos(phase + k * math.pi / 3), radius * math.sin(phase + k * math.pi / 3)) for k in range(6)]
            stations = [(point, 'ring') for point in ring if all((math.dist(point, q) > 1.0 for q in scanned))]
            if covers(grid, covered_now, stations):
                seeds.append(stations)
    return seeds

def optimise_station_set(grid, covered_now, stations, start, source_route, scan_cost_s, terms=None, iterations=PLAN_SEARCH_STEPS):
    """Greedy geometry search: slide stations towards sources while coverage holds."""
    stations = list(stations)
    current = station_plan_cost(start, source_route, stations, scan_cost_s, terms)
    if not stations:
        return (stations, current)
    targets = sorted(source_route, key=lambda p: -math.hypot(p[0], p[1])) if source_route else []
    for step in range(iterations):
        index = step % len(stations)
        point, kind = stations[index]
        radius = math.hypot(point[0], point[1])
        options = []
        if targets:
            near = sorted(targets, key=lambda q: math.dist(point, q))[:2]
            for target in near:
                for fraction in (0.2, 0.45, 0.7):
                    options.append(slide_towards(point, target, fraction))
        for delta in (-260.0, -130.0, 130.0, 260.0):
            scale = (radius + delta) / max(radius, 1e-06)
            if 200.0 < radius + delta < 1800.0:
                options.append((point[0] * scale, point[1] * scale))
        for angle in (0.12, -0.12, 0.26, -0.26):
            if radius > 1.0:
                beta = math.atan2(point[1], point[0]) + angle
                options.append((radius * math.cos(beta), radius * math.sin(beta)))
        for alternative in options:
            trial = list(stations)
            trial[index] = (alternative, kind)
            if not covers(grid, covered_now, trial):
                continue
            cost = station_plan_cost(start, source_route, trial, scan_cost_s, terms)
            if cost < current - 1.0:
                stations, current = (trial, cost)
                break
    return (stations, current)

def repair_station_set(cover, grid, covered, stations, candidates, route):
    """Verify the planned stations against the real quadtree certificate."""
    probe = copy_cover(cover)
    for point, _ in stations:
        probe.add(point)
    pool = list(candidates)
    added = []
    while not probe.complete and len(stations) + len(added) < PLAN_MAX_STATIONS + 4:
        cell = max(probe.cells, key=lambda c: c[2])
        center = project_arena(cell[:2])
        best = None
        for index, (point, kind, mask) in enumerate(pool):
            if math.dist(point, center) > PLAN_DISK_RADIUS:
                continue
            gain = bin(grid.uncovered(mask, covered)).count('1')
            extra, position = cheapest_insertion(route, point)
            if best is None or gain > best[0]:
                best = (gain, index, extra, position)
        if best is None:
            extra, position = cheapest_insertion(route, center)
            point, kind, mask = (center, 'gap_repair', grid.mask(center))
        else:
            _, index, extra, position = best
            point, kind, mask = pool.pop(index)
        stations.append((point, kind))
        added.append(point)
        covered |= mask
        route.insert(position, point)
        probe.add(point)
    return (stations, probe.complete)

def choose_plan_joint(ctrl, cover, sources, unknown):
    """Cheapest station set among ring seeds and source-visit seeds."""
    grid = getattr(ctrl, 'station_grid', None)
    if grid is None:
        grid = ctrl.station_grid = StationGrid()
    scan_cost_s = 6.0 * len(unknown)
    start = ctrl.position
    covered = 0
    for point in cover.points:
        covered |= grid.mask(point)
    if covered == grid.full:
        return None
    source_route = [sources[c] for c in warm_route(start, sources)]
    terms = future_source_terms(ctrl, len(sources))
    pending_key = frozenset(sources)
    sticky = getattr(ctrl, 'planned_stations', None)
    if sticky and getattr(ctrl, 'planned_key', None) == pending_key:
        kept = [(point, kind) for point, kind in sticky if all((math.dist(point, q) > 1.0 for q in cover.points))]
        if kept and covers(grid, covered, kept):
            return build_station_plan(ctrl, cover, grid, sources, source_route, kept, scan_cost_s, covered)
    seeds = []
    options = ring_station_options(covered, cover.points, grid, source_route)
    ranked = []
    for option in options:
        stations = [(point, 'ring') for point in option['stations']]
        cost = station_plan_cost(start, source_route, stations, scan_cost_s, terms)
        ranked.append((cost, option, stations))
    ranked.sort(key=lambda entry: entry[0])
    for _, option, stations in ranked[:1]:
        seeds.append((stations, option))
    for stations in outer_station_seeds(grid, covered, cover.points, source_route):
        seeds.append((stations, None))
    if sticky:
        kept = [(point, kind) for point, kind in sticky if all((math.dist(point, q) > 1.0 for q in cover.points))]
        if kept and covers(grid, covered, kept):
            seeds.insert(0, (kept, None))
    if not seeds:
        return None
    ranked_seeds = sorted(seeds, key=lambda entry: station_plan_cost(start, source_route, entry[0], scan_cost_s, terms))
    best = None
    for stations, option in ranked_seeds[:PLAN_SEEDS_REFINED]:
        if not covers(grid, covered, stations):
            continue
        improved, estimate = optimise_station_set(grid, covered, stations, start, source_route, scan_cost_s, terms)
        improved, estimate = refine_station_set(grid, covered, improved, source_route, start, source_route, scan_cost_s, terms)
        if not covers(grid, covered, improved):
            continue
        if best is None or estimate < best[0]:
            best = (estimate, improved, option)
    if best is None:
        return None
    estimate, stations, option = best
    ctrl.planned_stations = list(stations)
    ctrl.planned_key = pending_key
    return build_station_plan(ctrl, cover, grid, sources, source_route, stations, scan_cost_s, covered, ring_radius=None if option is None else option['radius'])

def build_station_plan(ctrl, cover, grid, sources, source_route, stations, scan_cost_s, covered, ring_radius=None):
    """Turn a chosen station set into an executable, certificate-verified plan."""
    start = ctrl.position
    candidates = [(point, 'ring', grid.mask(point)) for point in ring_station_candidates()]
    candidates.extend(((point, 'at_source', grid.mask(point)) for point in source_route))
    route = [start] + list(source_route) + [p for p, kind in stations if kind != 'at_source']
    stations, verified = repair_station_set(cover, grid, covered, stations, candidates, route)
    tasks = dict(sources)
    kinds = {}
    for index, (point, kind) in enumerate(stations):
        key = 21 + index
        tasks[key] = point
        kinds[key] = kind
    plan_route, stats = solve_route(start, tasks, horizon=16)
    cost = stats['route_distance_m'] / 5 + len(stations) * scan_cost_s
    return dict(centers=tasks, route=plan_route, dp=stats, estimated_s=cost, plan_type='joint_cover_stations', scan_tasks=[21 + i for i in range(len(stations))], scan_task_kinds=kinds, verified=verified, ring_radius=ring_radius, station_kinds=[kind for _, kind in stations], stations_planned=len(stations))

def choose_plan(ctrl, cover, original, area_planner):
    plan = choose_plan_reference(ctrl, cover, original, area_planner)
    sources = {c: info['center'] for c, info in plan['beliefs'].items()}
    try:
        joint = choose_plan_joint(ctrl, cover, sources, plan['unknown_channels'])
    except Exception:
        joint = None
    terms = future_source_terms(ctrl, len(sources))
    reference_estimate = plan['estimated_s'] + future_source_penalty(ctrl.position, [plan['centers'][key] for key in plan['route']], terms)
    plan['reference_estimate_with_future_s'] = reference_estimate
    if joint is not None and joint['estimated_s'] < reference_estimate - PLAN_NEW_PLAN_MARGIN_S:
        joint.update(beliefs=plan['beliefs'], unknown_channels=plan['unknown_channels'], reference_estimated_s=plan['estimated_s'], alternatives=plan['alternatives'], beam_states=plan['beam_states'], circle_candidate_count=plan['circle_candidate_count'], anticipated_scan_sources=plan['anticipated_scan_sources'])
        ctrl.plan_origin = 'joint'
        return joint
    ctrl.plan_origin = 'reference'
    return plan

def scheme_d_enhanced(ctrl):
    cover = AreaCover()
    ctrl.area_cover = cover
    ctrl.source_positions = {}
    ctrl.enhanced_events = []
    ctrl.enhanced_settings = dict(minimum_plan_saving_s=MIN_PLAN_SAVING_S, shared_view_min_channels=SHARED_VIEW_MIN_CHANNELS)
    original = coverage_points()
    area_planner = CirclePlanner(cover)
    synced = 0

    def sync_cover():
        nonlocal synced
        for p in cover.points[synced:]:
            area_planner.scanned(p)
        synced = len(cover.points)

    def scan(point, reason):
        before = ctrl.action_count
        cleared_before = sorted(ctrl.cleared)
        channels = [c for c in range(1, 21) if c not in ctrl.cleared and c not in ctrl.observations]
        scan_unknown(ctrl, cover, point)
        sync_cover()
        if ctrl.action_count > before:
            ctrl.enhanced_events.append(dict(kind='scan', point=list(point), reason=reason, after_action=before, end_action=ctrl.action_count, channels=channels, cleared_before=cleared_before, cleared_after=sorted(ctrl.cleared)))
    scan((0.0, 0.0), 'initial_census')
    if len(ctrl.pending()) >= SHARED_VIEW_MIN_CHANNELS:
        bearings = [math.radians(ctrl.observations[c][-1][1]) for c in ctrl.pending()]
        angle = max((k * math.pi / 36 for k in range(72)), key=lambda t: sum((min(abs(math.sin(t - b)), 0.5) for b in bearings)))
        point = (150 * math.cos(angle), 150 * math.sin(angle))
        ctrl.phase = 'enhanced_shared_probe'
        ctrl.enhanced_events.append(dict(kind='shared_probe', after_action=ctrl.action_count, point=list(point), channels=sorted(ctrl.pending())))
        for c in sorted(ctrl.pending(), key=lambda c: (c != ctrl.channel, c)):
            response = ctrl.measure(point, c)
            if response['measure_result'] == 'near' and (not ctrl.clear(point, c)):
                raise LocalizationError('Near shared probe contradicted by clearance')

    def at_stop(point):
        if cover.gain(point) >= DISCOVERY_GAIN_M2:
            scan(point, 'opportunistic_clearance_stop')
        shared_bearings(ctrl, point)
    for _ in range(160):
        if ctrl.certified():
            return
        plan = choose_plan(ctrl, cover, original, area_planner)
        route, centers, infos = (plan['route'], plan['centers'], plan['beliefs'])
        for channel, info in infos.items():
            ctrl.source_positions[channel] = info['center']
        action = route[0]
        ctrl.decisions.append(dict(after_action=ctrl.action_count, position=list(ctrl.position), cleared=sorted(ctrl.cleared), coverage_cells=len(cover.cells), beliefs={str(c): dict(center=list(i['center']), radius=i['radius'], version=ctrl.version[c], bearing_count=len(ctrl.observations[c]), no_signal_count=len(ctrl.no_signals[c])) for c, i in infos.items()}, candidate_points={str(c): list(p) for c, p in centers.items()}, planned_route=route, chosen_action=action, dp=plan['dp'], plan_type=plan['plan_type'], estimated_remaining_travel_scan_s=plan['estimated_s'], reference_estimated_travel_scan_s=plan['reference_estimated_s'], circle_candidates=plan['circle_candidate_count'], coverage_beam_states=plan['beam_states'], alternatives=plan['alternatives'], anticipated_scan_sources=plan['anticipated_scan_sources'], unknown_channels=plan['unknown_channels'], scan_tasks=plan['scan_tasks']))
        ctrl.decisions[-1]['scan_task_kinds'] = plan['scan_task_kinds']
        if action <= 20:
            fast_service(ctrl, action, at_stop)
        else:
            scan(centers[action], plan['plan_type'])
            ctrl.scanned.add(len(cover.points) - 1)
            shared_bearings(ctrl, centers[action])
    raise LocalizationError('Enhanced planner safety limit reached without completion')

def q4_clip(poly, normal, limit):
    out = []
    limit += 1e-07
    for a, b in zip(poly, poly[1:] + poly[:1]):
        fa = sum((x * y for x, y in zip(a, normal))) - limit
        fb = sum((x * y for x, y in zip(b, normal))) - limit
        if fa <= 0:
            out.append(a)
        if fa < 0 < fb or fb < 0 < fa:
            t = fa / (fa - fb)
            out.append((a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])))
    return out

def q4_posterior(ctrl, channel):
    key = (channel, ctrl.version[channel])
    if key in ctrl.summary_cache:
        return ctrl.summary_cache[key]
    poly = [(-1800.0, -1800.0), (1800.0, -1800.0), (1800.0, 1800.0), (-1800.0, 1800.0)]
    for p, angle in ctrl.observations[channel]:
        aa, bb = bearing_halfplanes([p], [angle], 1.005)
        for a, b in zip(aa, bb):
            poly = q4_clip(poly, a, b)
        for k in range(16):
            t = 2 * math.pi * k / 16
            n = (math.cos(t), math.sin(t))
            poly = q4_clip(poly, n, 1500.0 + sum((x * y for x, y in zip(n, p))))
    for k in range(32):
        t = 2 * math.pi * k / 32
        poly = q4_clip(poly, (math.cos(t), math.sin(t)), 1800.0)
    if not poly:
        raise RuntimeError('Empty positive-bearing intersection: inconsistent observations')
    circle = minimum_enclosing_circle(poly)
    answer = (poly, circle.center, circle.radius)
    ctrl.summary_cache[key] = answer
    return answer

def q4_cross(a, b, c):
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])

def q4_in_convex(p, poly, tolerance=1e-07):
    return len(poly) >= 3 and all((q4_cross(a, b, p) >= -tolerance for a, b in zip(poly, poly[1:] + poly[:1])))

def q4_segment_distance(p, a, b):
    d = (b[0] - a[0], b[1] - a[1])
    t = max(0.0, min(1.0, sum(((p[i] - a[i]) * d[i] for i in (0, 1))) / (d[0] ** 2 + d[1] ** 2)))
    return math.dist(p, (a[0] + t * d[0], a[1] + t * d[1]))

class q4_DirectionalCover:

    def __init__(self, spacing=990.0):
        if not 100 <= spacing < 1000:
            raise ValueError('spacing must be in [100,1000)')
        h = spacing * math.sqrt(3) / 2

        def p(i, j):
            return (spacing * (i + j / 2), h * j)
        triangles = []
        n = math.ceil(3600 / spacing) + 3
        for i in range(-n, n + 1):
            for j in range(-n, n + 1):
                for tri in ([p(i, j), p(i + 1, j), p(i, j + 1)], [p(i + 1, j), p(i + 1, j + 1), p(i, j + 1)]):
                    if q4_in_convex((0.0, 0.0), tri) or min((q4_segment_distance((0.0, 0.0), a, b) for a, b in zip(tri, tri[1:] + tri[:1]))) <= 1800.0 + 1e-07:
                        triangles.append(tri)
        self.triangles = triangles
        self.sites = sorted(set((p for tri in triangles for p in tri)))
        self.cells = list(range(len(triangles)))
        self.points = []
        self.visited = set()
        self.spacing = spacing

    @property
    def complete(self):
        return not self.cells

    def add(self, p):
        p = tuple(p)
        if p in self.visited:
            return
        self.visited.add(p)
        self.points.append(p)
        remaining = []
        for idx in self.cells:
            tri = self.triangles[idx]
            if all((v in self.visited for v in tri)):
                continue
            nearby = [s for s in self.points if max((math.dist(s, v) for v in tri)) <= 999.999]
            hull = convex_hull(nearby)
            if not all((q4_in_convex(v, hull) for v in tri)):
                remaining.append(idx)
        self.cells = remaining

    def candidates(self):
        return sorted({p for i in self.cells for p in self.triangles[i] if p not in self.visited})

def q4_optical_fallback(ctrl, channel):
    poly, _, _ = q4_posterior(ctrl, channel)
    t = math.radians(ctrl.observations[channel][0][1])
    u = (math.cos(t), math.sin(t))
    v = (-u[1], u[0])
    projections = [(sum((a * b for a, b in zip(p, u))), sum((a * b for a, b in zip(p, v)))) for p in poly]
    lo = [min((p[i] for p in projections)) for i in (0, 1)]
    hi = [max((p[i] for p in projections)) for i in (0, 1)]
    nx, ny = [max(1, math.ceil((hi[i] - lo[i]) / 20)) for i in (0, 1)]
    candidates = []
    for i in range(nx):
        for j in range(ny):
            x = lo[0] + (i + 0.5) * (hi[0] - lo[0]) / nx
            y = lo[1] + (j + 0.5) * (hi[1] - lo[1]) / ny
            candidates.append((x * u[0] + y * v[0], x * u[1] + y * v[1]))
    ctrl.fallback_count = getattr(ctrl, 'fallback_count', 0) + 1
    while candidates:
        p = min(candidates, key=lambda p: math.dist(p, ctrl.position))
        candidates.remove(p)
        if ctrl.clear(p, channel):
            return
    raise RuntimeError('Optical cover exhausted: observations or simulator violate assumptions')

def q4_localize(ctrl, channel, step=420.0, lateral=100.0):
    for iteration in range(10):
        poly, center, radius = q4_posterior(ctrl, channel)
        if radius <= 19.9:
            if not ctrl.clear(center, channel):
                raise RuntimeError('Certified 19.9 m enclosing circle failed clearance')
            return
        if len(ctrl.observations[channel]) >= 2:
            if ctrl.clear(center, channel):
                return
            ans = ctrl.measure(center, channel)
            if ans['measure_result'] == 'near':
                if ctrl.clear(center, channel):
                    return
            if ans['measure_result'] == 'direction':
                _, newcenter, newradius = q4_posterior(ctrl, channel)
                if newradius < radius * 0.8:
                    continue
        p, angle = ctrl.observations[channel][-1]
        t = math.radians(angle)
        u = (math.cos(t), math.sin(t))
        v = (-u[1], u[0])
        distance = math.dist(p, center)
        advance = min(step, max(20.0, distance * 0.65))
        baseline = min(lateral, max(12.0, distance * 0.25))
        candidates = [(p[0] + advance * u[0] + sign * baseline * v[0], p[1] + advance * u[1] + sign * baseline * v[1]) for sign in (1, -1)]
        candidates.sort(key=lambda q: math.dist(q, ctrl.position))
        improved = False
        for q in candidates:
            if any((math.dist(q, s) < 1e-05 for s in ctrl.no_signals[channel])):
                continue
            if any((math.dist(q, s) < 1e-05 for s, _ in ctrl.observations[channel])):
                continue
            ans = ctrl.measure(q, channel)
            if ans['measure_result'] == 'near':
                if ctrl.clear(q, channel):
                    return
            if ans['measure_result'] == 'direction':
                improved = True
                break
        if not improved:
            dark = min(candidates, key=lambda q: math.dist(q, ctrl.position))
            for _ in range(6):
                q = ((p[0] + dark[0]) / 2, (p[1] + dark[1]) / 2)
                if math.dist(q, p) < 1.0:
                    break
                ans = ctrl.measure(q, channel)
                if ans['measure_result'] == 'near':
                    if ctrl.clear(q, channel):
                        return
                if ans['measure_result'] == 'direction':
                    improved = True
                    break
                dark = q
        if not improved:
            break
    q4_optical_fallback(ctrl, channel)

def q4_route(start, points):
    pending = list(points)
    out = []
    p = start
    while pending:
        q = min(pending, key=lambda q: math.dist(p, q))
        pending.remove(q)
        out.append(q)
        p = q
    for _ in range(8):
        changed = False
        for i in range(len(out) - 1):
            a = start if i == 0 else out[i - 1]
            for j in range(i + 1, len(out)):
                old = math.dist(a, out[i])
                new = math.dist(a, out[j])
                if j + 1 < len(out):
                    old += math.dist(out[j], out[j + 1])
                    new += math.dist(out[i], out[j + 1])
                if new < old - 1e-06:
                    out[i:j + 1] = reversed(out[i:j + 1])
                    changed = True
        if not changed:
            break
    return out

def q4_policy_joint(ctrl, step=420.0, lateral=100.0, **_):
    cover = q4_RingCover()
    ctrl.area_cover = cover
    action = ('scan', (0.0, 0.0))
    while not ctrl.certified():
        kind, item = action
        if kind == 'scan':
            ctrl.phase = 'search'
            for c in ctrl.channels_to_scan():
                ans = ctrl.measure(item, c)
                if ans['measure_result'] == 'near':
                    ctrl.clear(item, c)
            cover.add(item)
            ctrl.scanned.add(len(cover.points) - 1)
        else:
            ctrl.phase = 'localize'
            q4_localize(ctrl, item, step, lateral)
        if ctrl.certified():
            break
        items = [('scan', p, p) for p in cover.candidates()]
        items += [('clear', c, q4_posterior(ctrl, c)[1]) for c in ctrl.pending()]
        if not items:
            raise RuntimeError('No action without certificate')
        ordered = q4_route(ctrl.position, [x[2] for x in items])
        first = next((x for x in items if x[2] == ordered[0]))
        action = first[:2]

def q4_ring_sites(inner=7, outer=14, r1=950.0, r2=2050.0, phase=0.0):
    return [(0.0, 0.0)] + [(r * math.cos(2 * math.pi * k / n + offset), r * math.sin(2 * math.pi * k / n + offset)) for n, r, offset in ((inner, r1, 0.0), (outer, r2, phase)) for k in range(n)]

class q4_RingCover:

    def __init__(self):
        self.sites = q4_ring_sites(8, 16, 960.0, 1850.0, math.pi / 16)
        self.points = []
        self.visited = set()
        self.cells = list(range(len(self.sites)))

    @property
    def complete(self):
        return not self.cells

    def add(self, p):
        p = tuple(p)
        if p not in self.visited:
            self.points.append(p)
            self.visited.add(p)
        self.cells = [i for i in self.cells if self.sites[i] not in self.visited]

    def candidates(self):
        return [self.sites[i] for i in self.cells]

def execute_q3(client, progress=None):
    ctrl = Controller(client, progress)
    started = time.perf_counter()
    error = None
    exit_response = None
    entered = False
    try:
        client.enter()
        entered = True
        scheme_d_enhanced(ctrl)
    except Exception as exc:
        error = f'{type(exc).__name__}: {exc}'
    finally:
        if entered and (not client.has_pending_action):
            try:
                exit_response = client.exit()
            except Exception as exc:
                error = (error + '; ' if error else '') + f'exit: {type(exc).__name__}: {exc}'
    result = ctrl.result()
    result['planning_method'] = 'joint coverage circle placement and receding-horizon route DP'
    result['planner'] = 'enhanced'
    result['enhanced_events'] = getattr(ctrl, 'enhanced_events', [])
    result['enhanced_settings'] = getattr(ctrl, 'enhanced_settings', {})
    result['decisions'] = ctrl.decisions
    result['dp_decision_count'] = len(ctrl.decisions)
    result['dp_total_states'] = sum((d['dp']['states'] for d in ctrl.decisions))
    result.update(scheme='d', real_runtime_s=time.perf_counter() - started, status='complete' if error is None and ctrl.certified() else 'incomplete', error=error, exit_acknowledged=exit_response is not None, unresolved_request=client.has_pending_action, virtual_target_s=500, within_virtual_target_500s=ctrl.certified() and ctrl.virtual_time <= 500)
    return result

def execute_q4(client, progress=None):
    ctrl = Controller(client, progress)
    started = time.perf_counter()
    error = None
    entered = False
    exited = False
    try:
        client.enter()
        entered = True
        q4_policy_joint(ctrl)
    except Exception as exc:
        error = f'{type(exc).__name__}: {exc}'
    finally:
        if entered and (not client.has_pending_action):
            try:
                client.exit()
                exited = True
            except Exception as exc:
                error = (error + '; ' if error else '') + f'exit: {type(exc).__name__}: {exc}'
    result = ctrl.result()
    result.update(problem=4, planner='joint', options={'mode': 'joint'}, status='complete' if error is None and ctrl.certified() and exited else 'incomplete', error=error, exit_acknowledged=exited, unresolved_request=client.has_pending_action, real_runtime_s=time.perf_counter() - started, fallback_count=getattr(ctrl, 'fallback_count', 0), certificate_method='offline-certified 25-site ring cover', planned_station_count=len(ctrl.area_cover.sites) if ctrl.area_cover else None)
    return result

def main(argv=None):
    parser = argparse.ArgumentParser(prog='robot_dog')
    parser.add_argument('--problem', type=int, choices=(3, 4), required=True)
    parser.add_argument('--robot-id', required=True)
    parser.add_argument('--base-url', default='http://127.0.0.1:2026')
    parser.add_argument('--output-dir', type=Path, default=Path('runs'))
    parser.add_argument('--quiet', action='store_true')
    args = parser.parse_args(argv)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    stamp = time.strftime('%Y%m%d-%H%M%S') + f'-{time.time_ns() % 1000000:06d}'
    prefix = args.output_dir / f'q{args.problem}-{stamp}'
    client = HTTPClient(args.robot_id, base_url=args.base_url, log_path=prefix.with_suffix('.jsonl'), exit_margin_s=10)
    progress = None if args.quiet else lambda s: print(s, flush=True)
    if args.problem == 3:
        result = execute_q3(client, progress)
        decisions = result.pop('decisions')
        decision_path = prefix.with_suffix('.decisions.jsonl')
        decision_path.write_text(''.join((json.dumps(d, ensure_ascii=False) + '\n' for d in decisions)), encoding='utf-8')
        result['decision_log_path'] = str(decision_path)
    else:
        result = execute_q4(client, progress)
    result['environment'] = 'official interface (user-started session)'
    prefix.with_suffix('.json').write_text(json.dumps(result, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0 if result['status'] == 'complete' else 2
if __name__ == '__main__':
    sys.exit(main())
