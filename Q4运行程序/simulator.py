import argparse
import copy
import hashlib
import json
import math
import random
import sys
import threading
import time
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from robot_dog import HTTPClient, execute_q3, execute_q4, finite_number, valid_identifier

class MockTransport:

    def __init__(self, seed=20260910, sources=None, robot_id='local-test', remaining_real_duration_s=1200, clock=time.monotonic, max_virtual_duration_s=360000, error_mode='hash'):
        if not valid_identifier(robot_id, 64):
            raise ValueError('Invalid robot_id')
        if not finite_number(remaining_real_duration_s) or not 0 <= remaining_real_duration_s <= 1200:
            raise ValueError('Invalid remaining real duration')
        self.__seed, self.__robot_id, self.__clock = (seed, robot_id, clock)
        if error_mode not in ('hash', 'positive', 'negative', 'smooth'):
            raise ValueError('Unknown error_mode')
        self.__error_mode = error_mode
        self.__duration = float(remaining_real_duration_s)
        self.__max_virtual = float(max_virtual_duration_s)
        if not math.isfinite(self.__max_virtual) or self.__max_virtual <= 0:
            raise ValueError('Invalid maximum virtual duration')
        rng = random.Random(seed)
        if sources is None:
            sources = []
            for channel in rng.sample(range(1, 21), rng.randint(10, 16)):
                r, a = (1800 * math.sqrt(rng.random()), 2 * math.pi * rng.random())
                sources.append({'position': (r * math.cos(a), r * math.sin(a)), 'channel': channel, 'radius': rng.uniform(1000, 1500)})
        self.__sources = {}
        for source in sources:
            p = source['position']
            x, y = (p['x'], p['y']) if isinstance(p, dict) else p
            c, radius = (source['channel'], source['radius'])
            if not all((finite_number(v) for v in (x, y, c, radius))):
                raise ValueError('Non-numeric source field')
            if math.hypot(x, y) > 1800 + 1e-09 or not 1000 <= radius <= 1500:
                raise ValueError('Synthetic source outside contest spatial/radius limits')
            if not float(c).is_integer() or not 1 <= c <= 20 or c in self.__sources:
                raise ValueError('Source channels must be distinct integers in 1..20')
            self.__sources[int(c)] = {'position': (float(x), float(y)), 'channel': int(c), 'radius': float(radius), 'cleared': False}
        self.__state = 'new'
        self.__position, self.__channel = ((0.0, 0.0), 1)
        self.__virtual_us = 0
        self.__started = None
        self.__ended = None
        self.__reason = None
        self.__cache = {}
        self.__lock = threading.Lock()
        self.__counts = {'measure_count': 0, 'clear_attempts': 0, 'cleared_count': 0, 'switch_count': 0, 'movement_m': 0.0}

    def __base(self, accepted):
        return {'accepted': accepted, 'real_timestamp_ms': int(time.time() * 1000), 'virtual_time_s': self.__virtual_us / 1000000 if accepted else 0}

    def __expire(self):
        if self.__state == 'active':
            if self.__clock() >= self.__started + self.__duration:
                self.__state, self.__reason, self.__ended = ('ended', 'real_timeout', self.__clock())
            elif self.__virtual_us >= self.__max_virtual * 1000000:
                self.__state, self.__reason, self.__ended = ('ended', 'virtual_timeout', self.__clock())

    def __validate(self, path, payload):
        if path not in ('/enter', '/measure', '/clear', '/exit'):
            return 404
        if not isinstance(payload, dict):
            return 400
        fields = {'arena_id', 'robot_id', 'request_id'}
        if path in ('/measure', '/clear'):
            fields |= {'position', 'channel'}
        if not fields.issubset(payload):
            return 400
        if set(payload) - fields:
            return 200
        if not isinstance(payload['arena_id'], str):
            return 400
        if not valid_identifier(payload['robot_id'], 64) or not valid_identifier(payload['request_id'], 128):
            return 400
        if path in ('/measure', '/clear'):
            p, c = (payload['position'], payload['channel'])
            if not isinstance(p, dict) or not {'x', 'y'}.issubset(p):
                return 400
            if set(p) - {'x', 'y'}:
                return 200
            if not all((finite_number(p[k]) and abs(p[k]) <= 2000000 for k in ('x', 'y'))):
                return 400
            if not finite_number(c) or not float(c).is_integer() or (not 1 <= c <= 20):
                return 400
        if payload['arena_id'] != 'default' or payload['robot_id'] != self.__robot_id:
            return 200
        return None

    def __bearing(self, point, source):
        dx, dy = (source['position'][0] - point[0], source['position'][1] - point[1])
        true_deg = math.degrees(math.atan2(dy, dx)) % 360
        key = f"{self.__seed}|{source['channel']}|{float(point[0])!r}|{float(point[1])!r}".encode()
        integer = int.from_bytes(hashlib.sha256(key).digest()[:8], 'big')
        lo, hi = (math.ceil((true_deg - 1) * 100), math.floor((true_deg + 1) * 100))
        if self.__error_mode == 'positive':
            value = hi
        elif self.__error_mode == 'negative':
            value = lo
        elif self.__error_mode == 'smooth':
            value = min(hi, max(lo, round((true_deg + math.sin(point[0] * 0.011 + point[1] * 0.017)) * 100)))
        else:
            value = lo + integer % (hi - lo + 1)
        return value % 36000 / 100

    def post(self, path, payload):
        if not self.__lock.acquire(blocking=False):
            return (409, self.__base(False))
        try:
            error = self.__validate(path, payload)
            if error is not None:
                return (error, self.__base(False))
            fingerprint = path + '\n' + json.dumps(payload, sort_keys=True, separators=(',', ':'), allow_nan=False)
            rid = payload['request_id']
            if rid in self.__cache:
                original, answer = self.__cache[rid]
                return copy.deepcopy(answer) if original == fingerprint else (409, self.__base(False))
            self.__expire()
            if self.__state == 'ended':
                raise ConnectionError('Synthetic test ended; local interface is closed')
            if path == '/enter':
                if self.__state != 'new':
                    return (200, self.__base(False))
                self.__state, self.__started = ('active', self.__clock())
                response = self.__base(True)
                response.update(max_real_duration_s=1200, max_virtual_duration_s=self.__max_virtual, remaining_real_duration_s=self.__duration)
            elif self.__state != 'active':
                return (200, self.__base(False))
            elif path == '/exit':
                self.__state, self.__reason, self.__ended = ('ended', 'user_exit', self.__clock())
                response = self.__base(True)
                response['exit_reason'] = 'user_exit'
            else:
                point = (float(payload['position']['x']), float(payload['position']['y']))
                channel = int(payload['channel'])
                distance = math.dist(self.__position, point)
                self.__virtual_us += round(distance / 5 * 1000000)
                self.__counts['movement_m'] += distance
                self.__position = point
                source = self.__sources.get(channel)
                live = source is not None and (not source['cleared'])
                source_distance = math.dist(point, source['position']) if live else math.inf
                if path == '/measure':
                    switched = int(self.__channel != channel)
                    self.__counts['switch_count'] += switched
                    self.__counts['measure_count'] += 1
                    self.__channel = channel
                    self.__virtual_us += (5 + switched) * 1000000
                    response = self.__base(True)
                    if source_distance > (source['radius'] if live else -1):
                        response['measure_result'] = 'no_signal'
                    elif source_distance <= 5:
                        response['measure_result'] = 'near'
                    else:
                        response.update(measure_result='direction', svd_deg=self.__bearing(point, source))
                else:
                    self.__counts['clear_attempts'] += 1
                    success = source_distance <= 20
                    if success:
                        source['cleared'] = True
                        self.__counts['cleared_count'] += 1
                    self.__virtual_us += (5 if success else 3) * 1000000
                    response = self.__base(True)
                    response['clear_result'] = 'success' if success else 'no_target_in_range'
            answer = (200, response)
            self.__cache[rid] = (fingerprint, copy.deepcopy(answer))
            return answer
        finally:
            self.__lock.release()

    def report(self):
        with self.__lock:
            self.__expire()
            if self.__state != 'ended':
                raise RuntimeError('Synthetic truth is only available AFTER the run ends')
            total = len(self.__sources)
            cleared = self.__counts['cleared_count']
            virtual = self.__virtual_us / 1000000
            return {'environment': 'self-built synthetic simulator; not official', 'seed': self.__seed, 'source_count': total, 'cleared_count': cleared, 'clearance_ratio': cleared / total if total else None, 'virtual_time_s': virtual, 'average_clear_time_s': virtual / cleared if cleared else None, 'program_runtime_s': self.__ended - self.__started if self.__started is not None else 0, 'exit_reason': self.__reason, **copy.deepcopy(self.__counts), 'sources': copy.deepcopy(list(self.__sources.values()))}

def make_server(transport, port=2027):

    class Handler(BaseHTTPRequestHandler):

        def log_message(self, *args):
            pass

        def send_json(self, status, response):
            body = json.dumps(response, allow_nan=False).encode('utf-8')
            self.send_response(status)
            self.send_header('Content-Type', 'application/json; charset=utf-8')
            self.send_header('Content-Length', str(len(body)))
            self.end_headers()
            self.wfile.write(body)

        def error(self, status):
            self.send_json(status, {'accepted': False, 'real_timestamp_ms': int(time.time() * 1000), 'virtual_time_s': 0})

        def do_GET(self):
            self.error(405 if self.path in ('/enter', '/exit', '/measure', '/clear') else 404)

        def do_POST(self):
            content_type = self.headers.get('Content-Type', '').lower().replace(' ', '')
            if content_type not in ('application/json', 'application/json;charset=utf-8') or self.headers.get('Content-Encoding', 'identity') != 'identity':
                self.error(415)
                return
            try:
                size = int(self.headers.get('Content-Length', '-1'))
                if size > 65536:
                    self.error(413)
                    return
                if size < 0:
                    raise ValueError('missing content length')

                def unique(pairs):
                    obj = {}
                    for k, v in pairs:
                        if k in obj:
                            raise ValueError('duplicate key')
                        obj[k] = v
                    return obj
                payload = json.loads(self.rfile.read(size).decode('utf-8'), object_pairs_hook=unique, parse_constant=lambda value: (_ for _ in ()).throw(ValueError(value)))

                def depth(value):
                    return 1 + max((depth(v) for v in (value.values() if isinstance(value, dict) else value)), default=0) if isinstance(value, (list, dict)) else 0
                if depth(payload) > 16:
                    raise ValueError('nested JSON too deep')
            except (ValueError, UnicodeError, RecursionError):
                self.error(400)
                return
            try:
                status, response = transport.post(self.path, payload)
            except ConnectionError:
                self.close_connection = True
                return
            self.send_json(status, response)
    return ThreadingHTTPServer(('127.0.0.1', port), Handler)

def sources_for(seed, kind='uniform', n=None):
    rng = random.Random(seed)
    n = n or rng.randint(10, 16)
    channels = rng.sample(range(1, 21), n)
    sources = []
    for i, c in enumerate(channels):
        a = rng.random() * 2 * math.pi
        if kind == 'uniform':
            r = 1800 * math.sqrt(rng.random())
        elif kind == 'boundary':
            r = 1800
        elif kind == 'cluster':
            r = 160 * math.sqrt(rng.random())
        elif kind == 'origin_near':
            r = 3 * rng.random()
        else:
            raise ValueError(kind)
        sources.append(dict(channel=c, position=(r * math.cos(a), r * math.sin(a)), radius=1000 if kind == 'boundary' else rng.uniform(1000, 1500)))
    if kind == 'boundary':
        sources[0]['position'] = (1800.0, 0.0)
        sources[1]['position'] = (-1800.0, 0.0)
    return sources

class MockTransport4:

    def __init__(self, seed=20260910, sources=None, robot_id='local-test', remaining_real_duration_s=1200, clock=time.monotonic, max_virtual_duration_s=360000, error_mode='hash'):
        if not valid_identifier(robot_id, 64):
            raise ValueError('Invalid robot_id')
        if not finite_number(remaining_real_duration_s) or not 0 <= remaining_real_duration_s <= 1200:
            raise ValueError('Invalid remaining real duration')
        self.__seed, self.__robot_id, self.__clock = (seed, robot_id, clock)
        if error_mode not in ('hash', 'positive', 'negative', 'smooth'):
            raise ValueError('Unknown error_mode')
        self.__error_mode = error_mode
        self.__duration = float(remaining_real_duration_s)
        self.__max_virtual = float(max_virtual_duration_s)
        if not math.isfinite(self.__max_virtual) or self.__max_virtual <= 0:
            raise ValueError('Invalid maximum virtual duration')
        rng = random.Random(seed)
        if sources is None:
            sources = []
            for channel in rng.sample(range(1, 21), rng.randint(10, 16)):
                r, a = (1800 * math.sqrt(rng.random()), 2 * math.pi * rng.random())
                sources.append({'position': (r * math.cos(a), r * math.sin(a)), 'channel': channel, 'radius': rng.uniform(1000, 1500), 'heading': None if len(sources) % 2 == 0 else rng.uniform(0, 360)})
        self.__sources = {}
        for source in sources:
            p = source['position']
            x, y = (p['x'], p['y']) if isinstance(p, dict) else p
            c, radius = (source['channel'], source['radius'])
            if not all((finite_number(v) for v in (x, y, c, radius))):
                raise ValueError('Non-numeric source field')
            if math.hypot(x, y) > 1800 + 1e-09 or not 1000 <= radius <= 1500:
                raise ValueError('Synthetic source outside contest spatial/radius limits')
            if not float(c).is_integer() or not 1 <= c <= 20 or c in self.__sources:
                raise ValueError('Source channels must be distinct integers in 1..20')
            self.__sources[int(c)] = {'position': (float(x), float(y)), 'channel': int(c), 'radius': float(radius), 'heading': source.get('heading'), 'cleared': False}
            heading = source.get('heading')
            if heading is not None and (not finite_number(heading)):
                raise ValueError('Heading must be finite or None')
        self.__state = 'new'
        self.__position, self.__channel = ((0.0, 0.0), 1)
        self.__virtual_us = 0
        self.__started = None
        self.__ended = None
        self.__reason = None
        self.__cache = {}
        self.__lock = threading.Lock()
        self.__counts = {'measure_count': 0, 'clear_attempts': 0, 'cleared_count': 0, 'switch_count': 0, 'movement_m': 0.0}

    def __base(self, accepted):
        return {'accepted': accepted, 'real_timestamp_ms': int(time.time() * 1000), 'virtual_time_s': self.__virtual_us / 1000000 if accepted else 0}

    def __expire(self):
        if self.__state == 'active':
            if self.__clock() >= self.__started + self.__duration:
                self.__state, self.__reason, self.__ended = ('ended', 'real_timeout', self.__clock())
            elif self.__virtual_us >= self.__max_virtual * 1000000:
                self.__state, self.__reason, self.__ended = ('ended', 'virtual_timeout', self.__clock())

    def __validate(self, path, payload):
        if path not in ('/enter', '/measure', '/clear', '/exit'):
            return 404
        if not isinstance(payload, dict):
            return 400
        fields = {'arena_id', 'robot_id', 'request_id'}
        if path in ('/measure', '/clear'):
            fields |= {'position', 'channel'}
        if not fields.issubset(payload):
            return 400
        if set(payload) - fields:
            return 200
        if not isinstance(payload['arena_id'], str):
            return 400
        if not valid_identifier(payload['robot_id'], 64) or not valid_identifier(payload['request_id'], 128):
            return 400
        if path in ('/measure', '/clear'):
            p, c = (payload['position'], payload['channel'])
            if not isinstance(p, dict) or not {'x', 'y'}.issubset(p):
                return 400
            if set(p) - {'x', 'y'}:
                return 200
            if not all((finite_number(p[k]) and abs(p[k]) <= 2000000 for k in ('x', 'y'))):
                return 400
            if not finite_number(c) or not float(c).is_integer() or (not 1 <= c <= 20):
                return 400
        if payload['arena_id'] != 'default' or payload['robot_id'] != self.__robot_id:
            return 200
        return None

    def __bearing(self, point, source):
        dx, dy = (source['position'][0] - point[0], source['position'][1] - point[1])
        true_deg = math.degrees(math.atan2(dy, dx)) % 360
        key = f"{self.__seed}|{source['channel']}|{float(point[0])!r}|{float(point[1])!r}".encode()
        integer = int.from_bytes(hashlib.sha256(key).digest()[:8], 'big')
        lo, hi = (math.ceil((true_deg - 1) * 100), math.floor((true_deg + 1) * 100))
        if self.__error_mode == 'positive':
            value = hi
        elif self.__error_mode == 'negative':
            value = lo
        elif self.__error_mode == 'smooth':
            value = min(hi, max(lo, round((true_deg + math.sin(point[0] * 0.011 + point[1] * 0.017)) * 100)))
        else:
            value = lo + integer % (hi - lo + 1)
        return value % 36000 / 100

    def post(self, path, payload):
        if not self.__lock.acquire(blocking=False):
            return (409, self.__base(False))
        try:
            error = self.__validate(path, payload)
            if error is not None:
                return (error, self.__base(False))
            fingerprint = path + '\n' + json.dumps(payload, sort_keys=True, separators=(',', ':'), allow_nan=False)
            rid = payload['request_id']
            if rid in self.__cache:
                original, answer = self.__cache[rid]
                return copy.deepcopy(answer) if original == fingerprint else (409, self.__base(False))
            self.__expire()
            if self.__state == 'ended':
                raise ConnectionError('Synthetic test ended; local interface is closed')
            if path == '/enter':
                if self.__state != 'new':
                    return (200, self.__base(False))
                self.__state, self.__started = ('active', self.__clock())
                response = self.__base(True)
                response.update(max_real_duration_s=1200, max_virtual_duration_s=self.__max_virtual, remaining_real_duration_s=self.__duration)
            elif self.__state != 'active':
                return (200, self.__base(False))
            elif path == '/exit':
                self.__state, self.__reason, self.__ended = ('ended', 'user_exit', self.__clock())
                response = self.__base(True)
                response['exit_reason'] = 'user_exit'
            else:
                point = (float(payload['position']['x']), float(payload['position']['y']))
                channel = int(payload['channel'])
                distance = math.dist(self.__position, point)
                self.__virtual_us += round(distance / 5 * 1000000)
                self.__counts['movement_m'] += distance
                self.__position = point
                source = self.__sources.get(channel)
                live = source is not None and (not source['cleared'])
                source_distance = math.dist(point, source['position']) if live else math.inf
                if path == '/measure':
                    switched = int(self.__channel != channel)
                    self.__counts['switch_count'] += switched
                    self.__counts['measure_count'] += 1
                    self.__channel = channel
                    self.__virtual_us += (5 + switched) * 1000000
                    response = self.__base(True)
                    visible = live and (source['heading'] is None or (point[0] - source['position'][0]) * math.cos(math.radians(source['heading'])) + (point[1] - source['position'][1]) * math.sin(math.radians(source['heading'])) >= -1e-10)
                    if not visible or source_distance > source['radius']:
                        response['measure_result'] = 'no_signal'
                    elif source_distance <= 5:
                        response['measure_result'] = 'near'
                    else:
                        response.update(measure_result='direction', svd_deg=self.__bearing(point, source))
                else:
                    self.__counts['clear_attempts'] += 1
                    success = source_distance <= 20
                    if success:
                        source['cleared'] = True
                        self.__counts['cleared_count'] += 1
                    self.__virtual_us += (5 if success else 3) * 1000000
                    response = self.__base(True)
                    response['clear_result'] = 'success' if success else 'no_target_in_range'
            answer = (200, response)
            self.__cache[rid] = (fingerprint, copy.deepcopy(answer))
            return answer
        finally:
            self.__lock.release()

    def report(self):
        with self.__lock:
            self.__expire()
            if self.__state != 'ended':
                raise RuntimeError('Synthetic truth is only available AFTER the run ends')
            total = len(self.__sources)
            cleared = self.__counts['cleared_count']
            virtual = self.__virtual_us / 1000000
            return {'environment': 'self-built synthetic simulator; not official', 'seed': self.__seed, 'source_count': total, 'cleared_count': cleared, 'clearance_ratio': cleared / total if total else None, 'virtual_time_s': virtual, 'average_clear_time_s': virtual / cleared if cleared else None, 'program_runtime_s': self.__ended - self.__started if self.__started is not None else 0, 'exit_reason': self.__reason, **copy.deepcopy(self.__counts), 'sources': copy.deepcopy(list(self.__sources.values()))}

def sources_for_q4(seed, scenario='mixed', count=None):
    rng = random.Random(seed)
    n = count or rng.randint(10, 16)
    channels = rng.sample(range(1, 21), n)
    out = []
    for i, c in enumerate(channels):
        a = rng.uniform(0, 2 * math.pi)
        r = 1800 * math.sqrt(rng.random())
        if scenario == 'outward':
            r = rng.uniform(1770, 1800)
        if scenario == 'cluster':
            r = 300 * math.sqrt(rng.random())
        if scenario == 'tangent':
            r = 1800
        direction = None if scenario == 'omni' or (scenario in ('mixed', 'cluster') and i % 2 == 0) else rng.uniform(0, 360)
        if scenario == 'outward':
            direction = math.degrees(a)
        if scenario == 'tangent':
            direction = math.degrees(a) + 90
        out.append(dict(position=(r * math.cos(a), r * math.sin(a)), channel=c, radius=1000 if scenario in ('outward', 'tangent', 'all_directional') else rng.uniform(1000, 1500), heading=direction))
    return out
Q3_SCENARIOS = ('uniform', 'boundary', 'cluster', 'origin_near')
Q4_SCENARIOS = ('mixed', 'omni', 'all_directional', 'outward', 'tangent', 'cluster')

def run_q3_case(scenario, n, seed):
    sources = sources_for(seed, scenario, n)
    transport = MockTransport(seed=seed, sources=sources, robot_id='local-test')
    log_path = Path('runs') / f'sim-q3-{scenario}-{seed}.jsonl'
    client = HTTPClient('local-test', transport=transport, log_path=log_path)
    result = execute_q3(client)
    try:
        report = transport.report()
    except Exception:
        report = {'source_count': len(sources), 'cleared_count': result['cleared_count'], 'virtual_time_s': result['virtual_time_s']}
    ok = result['status'] == 'complete' and report['cleared_count'] == report['source_count']
    return dict(problem=3, scenario=scenario, seed=seed, source_count=report['source_count'], cleared_count=report['cleared_count'], virtual_time_s=report['virtual_time_s'], per_source_s=report['virtual_time_s'] / report['cleared_count'] if report['cleared_count'] else None, status=result['status'], error=result['error'], ok=ok)

def run_q4_case(scenario, n, seed):
    sources = sources_for_q4(seed, scenario, n)
    transport = MockTransport4(seed=seed, sources=sources, robot_id='local-test')
    log_path = Path('runs') / f'sim-q4-{scenario}-{seed}.jsonl'
    client = HTTPClient('local-test', transport=transport, log_path=log_path)
    result = execute_q4(client)
    try:
        report = transport.report()
    except Exception:
        report = {'source_count': len(sources), 'cleared_count': result['cleared_count'], 'virtual_time_s': result['virtual_time_s']}
    ok = result['status'] == 'complete' and report['cleared_count'] == report['source_count']
    return dict(problem=4, scenario=scenario, seed=seed, source_count=report['source_count'], cleared_count=report['cleared_count'], virtual_time_s=report['virtual_time_s'], per_source_s=report['virtual_time_s'] / report['cleared_count'] if report['cleared_count'] else None, status=result['status'], error=result['error'], fallback_count=result['fallback_count'], ok=ok)

def summarize(name, rows):
    good = [r for r in rows if r['ok']]
    per = [r['virtual_time_s'] / r['cleared_count'] for r in good]
    return dict(group=name, cases=len(rows), success=len(good), success_rate=len(good) / len(rows), avg_virtual_time_s=sum((r['virtual_time_s'] for r in good)) / len(good) if good else None, avg_per_source_s=sum(per) / len(per) if per else None, rows=rows)

def bench_q3(cases, seed_base):
    groups = []
    for si, scenario in enumerate(('uniform', 'boundary', 'cluster')):
        rows = [run_q3_case(scenario, (10, 13, 16)[i % 3], seed_base + si * 1000 + i) for i in range(cases)]
        groups.append(summarize(f'q3 {scenario}', rows))
    return {'groups': groups, 'ok': all((g['success'] == g['cases'] for g in groups))}

def bench_q4(cases, seed_base):
    groups = []
    for si, scenario in enumerate(Q4_SCENARIOS):
        rows = [run_q4_case(scenario, None, seed_base + si * 1000 + i) for i in range(cases)]
        groups.append(summarize(f'q4 {scenario}', rows))
    return {'groups': groups, 'ok': all((g['success'] == g['cases'] for g in groups))}

def main(argv=None):
    parser = argparse.ArgumentParser(prog='simulator')
    parser.add_argument('--problem', type=int, choices=(3, 4), required=True)
    parser.add_argument('--scenario', default=None, choices=('uniform', 'boundary', 'cluster', 'origin_near', 'mixed', 'omni', 'all_directional', 'outward', 'tangent'))
    parser.add_argument('--sources', type=int, default=None)
    parser.add_argument('--seed', type=int, default=20260910)
    parser.add_argument('--bench', action='store_true')
    parser.add_argument('--cases', type=int, default=6)
    parser.add_argument('--seed-base', type=int, default=2026091900)
    args = parser.parse_args(argv)
    if args.problem == 3:
        scenario = args.scenario or 'uniform'
        if scenario not in Q3_SCENARIOS:
            parser.error(f'Q3 scenarios: {Q3_SCENARIOS}')
    else:
        scenario = args.scenario or 'mixed'
        if scenario not in Q4_SCENARIOS:
            parser.error(f'Q4 scenarios: {Q4_SCENARIOS}')
    if args.bench:
        summary = bench_q3(args.cases, args.seed_base) if args.problem == 3 else bench_q4(args.cases, args.seed_base)
        print(json.dumps(summary, ensure_ascii=False, indent=2))
        return 0 if summary['ok'] else 2
    if args.problem == 3:
        row = run_q3_case(scenario, args.sources, args.seed)
    else:
        row = run_q4_case(scenario, args.sources, args.seed)
    print(json.dumps(row, ensure_ascii=False, indent=2))
    return 0 if row['ok'] else 2
if __name__ == '__main__':
    sys.exit(main())
