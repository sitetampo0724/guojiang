#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
离线蒙特卡洛自检 —— 不需要模拟器即可验证策略正确性与耗时。

用法:
    python selftest.py --problem 3 --trials 200
    python selftest.py --problem 4 --trials 200
    python selftest.py --coverage 20000     # 仅验证普查几何覆盖(问题3/4)

自检内容:
    1. 用与模拟器一致的物理规则(示向度误差、有效半径 1000-1500m、
       近距离阈值、清除半径、定向覆盖等)模拟干扰源;
    2. 跑完整策略(普查 -> 最近邻 -> 逼近定位 -> 清除);
    3. 统计: 完全清除率、平均/最差虚拟总时间、平均每源时间、普查覆盖检查。
"""

import argparse
import math
import random
import statistics
import sys

import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import robot as R


class SimJammer:
    def __init__(self, pos, channel, radius, directional=False, heading=0.0):
        self.pos = pos
        self.channel = channel
        self.radius = radius
        self.directional = directional
        self.heading = heading % 360.0
        self.cleared = False


class FakeClient:
    """与真实 Client 接口一致(enter/exit/measure/clear, pos/virtual_time)。"""

    def __init__(self, jammers, rng):
        self.jammers = jammers
        self.rng = rng
        self.pos = (0.0, 0.0)
        self.channel = 1
        self.virtual_time = 0.0
        self.remaining = 1200.0
        self.real_start = 0.0
        self._err_cache = {}        # 同一地点误差固定(与题目一致)
        self.n_measure = 0
        self.n_clear = 0

    # -- 物理 --
    def _error(self, pos):
        key = (round(pos[0], 1), round(pos[1], 1))
        if key not in self._err_cache:
            self._err_cache[key] = self.rng.uniform(-1.0, 1.0)
        return self._err_cache[key]

    def _find(self, channel):
        for j in self.jammers:
            if j.channel == channel and not j.cleared:
                return j
        return None

    def _detectable(self, pos, j):
        d = R.dist(pos, j.pos)
        if d > j.radius:
            return False
        if j.directional and abs(R.ang_diff(R.bearing(j.pos, pos), j.heading)) > 90.0:
            return False
        return True

    # -- 接口 --
    def enter(self):
        return {"accepted": True}

    def exit(self):
        return {"accepted": True}

    def real_time_left(self):
        return 1e9

    def measure(self, pos, channel):
        move = R.dist(self.pos, pos)
        self.virtual_time += move / R.MOVE_SPEED
        if channel != self.channel:
            self.virtual_time += R.SWITCH_SEC
        self.virtual_time += R.MEASURE_SEC
        self.pos = (float(pos[0]), float(pos[1]))
        self.channel = int(channel)
        self.n_measure += 1
        j = self._find(channel)
        if j is None or not self._detectable(self.pos, j):
            return "no_signal", None
        d = R.dist(self.pos, j.pos)
        if d <= R.NEAR_M:
            return "near", None
        true_b = R.bearing(self.pos, j.pos)
        return "direction", R.ang_norm(true_b + self._error(self.pos))

    def clear(self, pos, channel):
        move = R.dist(self.pos, pos)
        self.virtual_time += move / R.MOVE_SPEED
        self.pos = (float(pos[0]), float(pos[1]))
        self.n_clear += 1
        j = self._find(channel)
        if j is not None and R.dist(self.pos, j.pos) <= R.CLEAR_RADIUS:
            j.cleared = True
            self.virtual_time += R.CLEAR_OK_SEC
            return True
        self.virtual_time += R.CLEAR_FAIL_SEC
        return False


def random_point_in_disk(rng, radius):
    r = radius * math.sqrt(rng.random())
    a = rng.uniform(0.0, 360.0)
    return (r * math.cos(math.radians(a)), r * math.sin(math.radians(a)))


def make_case(rng, problem):
    n = rng.randint(10, 16)
    channels = rng.sample(range(R.CH_MIN, R.CH_MAX + 1), n)
    jammers = []
    for ch in channels:
        pos = random_point_in_disk(rng, 1800.0)
        radius = rng.uniform(1000.0, 1500.0)
        directional = (problem == 4 and rng.random() < 0.5)
        heading = rng.uniform(0.0, 360.0)
        jammers.append(SimJammer(pos, ch, radius, directional, heading))
    return jammers


def run_trial(problem, rng, verbose_fail=False):
    jammers = make_case(rng, problem)
    client = FakeClient(jammers, rng)
    try:
        R.run_strategy(client, problem=problem, verbose=False)
    except Exception as e:
        if verbose_fail:
            print("  异常: %r" % e)
        return None
    all_ok = all(j.cleared for j in jammers)
    return {
        "ok": all_ok,
        "n": len(jammers),
        "vtime": client.virtual_time,
        "per": client.virtual_time / len(jammers),
        "measures": client.n_measure,
        "clears": client.n_clear,
    }


def check_coverage(problem, samples, seed=1):
    """普查几何覆盖: 对随机(干扰源位置, 定向方向), 是否必有普查点
    在 1000m 内(最小接收半径)且位于有效覆盖角度范围内。"""
    rng = random.Random(seed)
    pts = R.census_points(problem)
    bad = 0
    for _ in range(samples):
        g = random_point_in_disk(rng, 1800.0)
        heading = rng.uniform(0.0, 360.0)
        found = False
        for p in pts:
            if R.dist(p, g) > 1000.0:
                continue
            # 问题3全为全向源; 问题4按最坏情况(定向源)检查覆盖角度
            if problem == 3 or abs(R.ang_diff(R.bearing(g, p), heading)) <= 90.0:
                found = True
                break
        if not found:
            bad += 1
    return bad, samples


def max_nearest_distance(problem, samples=200000, seed=7):
    """数值验证: 圆域内任一点到最近普查点的最大距离(应 < 1000m)。"""
    rng = random.Random(seed)
    pts = R.census_points(problem)
    worst = 0.0
    for _ in range(samples):
        g = random_point_in_disk(rng, 1800.0)
        d = min(R.dist(p, g) for p in pts)
        worst = max(worst, d)
    return worst


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--problem", type=int, choices=[3, 4], default=3)
    ap.add_argument("--trials", type=int, default=200)
    ap.add_argument("--coverage", type=int, default=20000,
                    help="普查覆盖检查的采样数(设为0可跳过)")
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    if args.coverage:
        for prob in (3, 4):
            bad, n = check_coverage(prob, args.coverage)
            print("[coverage Q%d] 未覆盖 %d / %d 采样 (%.4f%%)"
                  % (prob, bad, n, 100.0 * bad / n))
        print("[worst-dist Q3] 最近普查点最大距离 %.1f m (需<1000)"
              % max_nearest_distance(3))

    rng = random.Random(args.seed)
    results = []
    fails = 0
    for t in range(args.trials):
        r = run_trial(args.problem, rng)
        if r is None:
            fails += 1
            continue
        if not r["ok"]:
            fails += 1
        results.append(r)
        if (t + 1) % 50 == 0:
            print("  ... %d/%d" % (t + 1, args.trials))

    if not results:
        print("无有效试验")
        return
    ok = sum(1 for r in results if r["ok"])
    vt = [r["vtime"] for r in results]
    per = [r["per"] for r in results]
    print("\n===== 问题%d 蒙特卡洛 (%d 次) =====" % (args.problem, len(results)))
    print("完全清除率: %d/%d (%.1f%%), 异常/失败 %d 次"
          % (ok, len(results), 100.0 * ok / len(results), fails))
    print("虚拟总时间: 均值 %.0f s, 中位 %.0f s, 最差 %.0f s"
          % (statistics.mean(vt), statistics.median(vt), max(vt)))
    print("平均每源:   均值 %.0f s, 中位 %.0f s, 最差 %.0f s"
          % (statistics.mean(per), statistics.median(per), max(per)))
    print("平均检测次数 %.0f, 平均清除次数 %.1f"
          % (statistics.mean([r["measures"] for r in results]),
             statistics.mean([r["clears"] for r in results])))


if __name__ == "__main__":
    main()
