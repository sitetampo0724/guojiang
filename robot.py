#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CUMCM 2026 B 题 —— 机器狗自动定位清除程序

用法:
    python robot.py --robot <参赛队号> --problem 3|4 [--url http://127.0.0.1:2026]

前置条件:
    1. 模拟器已启动、已联网登录;
    2. 在模拟器中选定 问题3/问题4 的演练或正式测试, 并点击开始;
    3. 倒计时结束、接口开放后运行本程序(程序会自动重试 /enter)。

策略概述(详见论文):
    普查(几何保证覆盖的布局, 自适应停测) -> 扫全部20频道并记录示向度
    -> 按 2-opt 全局顺序逐源: 示向条形域裁剪求定位区域, 沿线推进/交点跳测/
       侧向跳收敛, 区域直径 <=60m 时探针网格清除(几何保证必中) -> 全部清除。
    定向源越界时退回最近有信号点爬行; 爬行失效时沿示向线强制步进(凸性保证)。
"""

import argparse
import json
import math
import os
import sys
import time
from datetime import datetime
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

# ---------------------------------------------------------------------------
# 基本几何
# ---------------------------------------------------------------------------

def ang_norm(a):
    return a % 360.0


def unit(a_deg):
    r = math.radians(a_deg)
    return (math.cos(r), math.sin(r))


def bearing(p, q):
    """p 指向 q 的方位角(度, [0,360)), x轴正向逆时针为正。"""
    return math.degrees(math.atan2(q[1] - p[1], q[0] - p[0])) % 360.0


def ang_diff(a, b):
    """有向角 a-b, 取值 (-180, 180]。"""
    d = (a - b) % 360.0
    return d - 360.0 if d > 180.0 else d


def dist(p, q):
    return math.hypot(q[0] - p[0], q[1] - p[1])


def vadd(p, v):
    return (p[0] + v[0], p[1] + v[1])


def vsub(p, q):
    return (p[0] - q[0], p[1] - q[1])


def vscale(v, s):
    return (v[0] * s, v[1] * s)


def vcross(a, b):
    return a[0] * b[1] - a[1] * b[0]


# ---------------------------------------------------------------------------
# 凸多边形工具(半平面裁剪)
# ---------------------------------------------------------------------------

def clip_halfplane(poly, p, a_deg, keep_positive=True, eps=1e-9):
    """用半平面 cross(unit(a), x-p) >= 0 (或 <= 0) 裁剪凸多边形。"""
    u = unit(a_deg)

    def h(x):
        c = vcross(u, vsub(x, p))
        return c >= -eps if keep_positive else c <= eps

    out = []
    n = len(poly)
    for i in range(n):
        cur = poly[i]
        nxt = poly[(i + 1) % n]
        hc, hn = h(cur), h(nxt)
        if hc:
            out.append(cur)
        if hc != hn:
            vc = vcross(u, vsub(cur, p))
            vn = vcross(u, vsub(nxt, p))
            t = vc / (vc - vn)
            out.append((cur[0] + (nxt[0] - cur[0]) * t,
                        cur[1] + (nxt[1] - cur[1]) * t))
    return out


def clip_strip(poly, p, b_deg, err_deg=1.0):
    """检测点 p 处示向度 b±err 的条形不确定域裁剪。"""
    poly = clip_halfplane(poly, p, b_deg - err_deg, True)
    poly = clip_halfplane(poly, p, b_deg + err_deg, False)
    return poly


def regular_gon(center, radius, n, start_deg=0.0):
    return [vadd(center, vscale(unit(start_deg + 360.0 * i / n), radius))
            for i in range(n)]


def polygon_diameter(poly):
    """多边形直径: 任意两顶点距离的最大值(凸多边形直径必在顶点对上取得)。"""
    n = len(poly)
    return max((dist(poly[i], poly[j])
                for i in range(n) for j in range(i + 1, n)), default=0.0)


def polygon_centroid(poly):
    n = len(poly)
    return (sum(p[0] for p in poly) / n, sum(p[1] for p in poly) / n)


# 目标区域(半径1800)略微外扩的36边形, 用作定位区域的初始裁剪域
ARENA_GON = regular_gon((0.0, 0.0), 1900.0, 36)

CH_MIN, CH_MAX = 1, 20
MOVE_SPEED = 5.0          # m/s
MEASURE_SEC = 5.0         # 每次检测耗时
SWITCH_SEC = 1.0          # 切换频道耗时
CLEAR_FAIL_SEC = 3.0      # 清除未发现耗时
CLEAR_OK_SEC = 5.0        # 清除成功耗时
NEAR_M = 5.0              # 近距离阈值
CLEAR_RADIUS = 20.0       # 清除半径
BEARING_ERR = 1.0         # 示向度误差界(度)
CENSUS_REGION_D = 60.0    # 普查中某频道定位区域直径 <=60m 即跳过重复测量
                          # (与探针覆盖保证衔接: d<=60 时探针网格必覆盖质心 60m 盘)


# ---------------------------------------------------------------------------
# 普查点布局
# ---------------------------------------------------------------------------

def census_points(problem):
    """普查路径点序列(已按行走顺序排列)。

    问题3: 中心 + 半径 1150 的六边形环。
        覆盖证明: 设环半径 R, 最坏点为圆域边缘且位于两相邻环点角平分线上,
        距离 f(R) = sqrt(1800^2 + R^2 - 2*1800*R*cos30°) <= 1000,
        解得 R >= 1123。取 R=1150(余量 11m, 最坏 988.7m < 1000m 最小接收半径),
        故所有全向干扰源必被检出; 且巡回路程 7R 随 R 减小而缩短, 1150 已近最优。
        (R=1559 是"边缘最坏点与内部最坏点同时最优"的解, 但对覆盖约束并非必要。)
    问题4: 中心 + 内环(100,6@38°) + 中环(975,6@38°) + 六环(1050,6@0°)
        + 外环(1908,12@25°), 共 31 点。定向覆盖判据: 对圆域内任意源位置 g,
        1000 m 内普查点对 g 的示向角集合的最大间隙 < 180°, 即任意定向
        方向的源正前方半平面内必有普查点。布局由环形族模式搜索求得
        (8 参数: 四环半径+起始角, 细网格约束最坏间隙 <= 170°,
        目标 = 2-opt 开放巡回长); 六种网格(径向步长 10m/5m/2m,
        含半格偏移)复核最坏间隙 169.99°~172.95° < 180°,
        100 万随机采样零漏检。巡回路程 22.0 km -> 18.9 km。
    """
    pts = [(0.0, 0.0)]
    if problem == 4:
        # 内环: 覆盖近圆心源(中心点单点无法覆盖定向方向)
        pts += regular_gon((0.0, 0.0), 100.0, 6, start_deg=38.0)
        # 中环(975, 6@38°): 与 1050 六环错开半径+角度,
        # 使 r~900-1300 源前方半平面必有近距离普查点
        pts += regular_gon((0.0, 0.0), 975.0, 6, start_deg=38.0)
    if problem == 3:
        # 环半径 1150: 巡回最短且 f(R)<=1000 仍成立(见 docstring 证明)
        pts += regular_gon((0.0, 0.0), 1150.0, 6)
    else:
        # 六环 1050: 全向覆盖条件 f(R)<=1000 只需 R>=1123(见问题3证明),
        # 与 975 中环协同满足定向覆盖(模式搜索在约束边界内取下限)
        pts += regular_gon((0.0, 0.0), 1050.0, 6)
    if problem == 4:
        # 外环 1908(12 点, start 25°): 边缘源 1000m 内仍有外环点
        pts += regular_gon((0.0, 0.0), 1908.0, 12, start_deg=25.0)
        # 巡回顺序: 模式搜索内置 2-opt 解(18.9 km)。先内环, 再 975/1050
        # 双环交错编织, 最后 1908 外环
        order = [0, 1, 6, 5, 4, 3, 2, 8, 15, 9, 16, 10, 17, 11, 18, 12,
                 13, 7, 14, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 19]
        pts = [pts[i] for i in order]
    return pts


def channel_region(dets):
    """由一组 (检测点, 示向度) 求定位区域多边形。"""
    region = list(ARENA_GON)
    for (p, b) in dets:
        region = clip_strip(region, p, b, BEARING_ERR)
    return region


# ---------------------------------------------------------------------------
# HTTP 客户端
# ---------------------------------------------------------------------------

class FatalError(Exception):
    """不可重试的错误(请求格式、冲突等), 必须终止本次测试。"""


class Client:
    """模拟器 HTTP 客户端: 串行发送、幂等重试、JSONL 日志、虚拟时间跟踪。"""

    def __init__(self, base_url, robot_id, log_path):
        self.base_url = base_url.rstrip("/")
        self.robot_id = robot_id
        self.log_path = log_path
        self.pos = (0.0, 0.0)
        self.channel = 1            # 测向机当前频道
        self.virtual_time = 0.0     # 最近 accepted 响应的虚拟时刻
        self.remaining = None       # /enter 返回的剩余现实秒数
        self.real_start = None
        self._seq = 0
        self._log_fp = open(log_path, "a", encoding="utf-8")

    # ---- 底层通信 ----

    def _next_id(self, prefix):
        self._seq += 1
        return "%s-%d" % (prefix, self._seq)

    def _write_log(self, record):
        record["wall_time"] = time.time()
        self._log_fp.write(json.dumps(record, ensure_ascii=False) + "\n")
        self._log_fp.flush()

    def _post(self, path, payload, prefix, retries=10):
        """发送一个动作。网络故障重试时复用同一 request_id 与请求体(幂等)。"""
        payload = dict(payload)
        payload["request_id"] = payload.get("request_id") or self._next_id(prefix)
        body = json.dumps(payload).encode("utf-8")
        attempt = 0
        while True:
            attempt += 1
            rec = {"path": path, "request": payload, "attempt": attempt}
            try:
                req = Request(self.base_url + path, data=body,
                              headers={"Content-Type": "application/json"},
                              method="POST")
                with urlopen(req, timeout=10) as resp:
                    data = json.loads(resp.read().decode("utf-8"))
                rec["response"] = data
                self._write_log(rec)
                return data
            except HTTPError as e:
                if e.code in (400, 409, 413, 415):
                    rec["http_error"] = e.code
                    try:
                        rec["error_body"] = e.read().decode("utf-8", "replace")
                    except Exception:
                        pass
                    self._write_log(rec)
                    raise FatalError("HTTP %d on %s: %s" % (e.code, path, payload))
                # 429/5xx: 退避后同请求重试
                rec["http_error"] = e.code
                self._write_log(rec)
                if attempt > retries:
                    raise
                time.sleep(min(2.0 * attempt, 10.0))
            except (URLError, TimeoutError, ConnectionError, OSError) as e:
                rec["net_error"] = repr(e)
                self._write_log(rec)
                if attempt > retries:
                    raise
                time.sleep(min(0.5 * attempt, 5.0))

    def _action(self, path, x=None, y=None, channel=None, prefix="act"):
        payload = {"arena_id": "default", "robot_id": self.robot_id}
        if x is not None:
            payload["position"] = {"x": float(x), "y": float(y)}
        if channel is not None:
            payload["channel"] = int(channel)
        return self._post(path, payload, prefix)

    # ---- 高层动作 ----

    def enter(self, wait_sec=90):
        """进入目标区域。倒计时期间接口未开放, 自动重试。"""
        t0 = time.time()
        while True:
            try:
                data = self._action("/enter", prefix="enter")
            except (URLError, OSError):
                if time.time() - t0 > wait_sec:
                    raise
                time.sleep(1.0)
                continue
            if data.get("accepted") is True:
                self.remaining = data.get("remaining_real_duration_s", 1200)
                self.real_start = time.time()
                self.virtual_time = data.get("virtual_time_s", 0.0)
                print("[enter] 剩余现实时间 %s s" % self.remaining, flush=True)
                return data
            if time.time() - t0 > wait_sec:
                raise FatalError("/enter 未被接受: %s" % data)
            time.sleep(1.0)

    def exit(self):
        try:
            data = self._action("/exit", prefix="exit")
            if data.get("accepted"):
                print("[exit] 虚拟时刻 %.1f s, 原因 %s"
                      % (self.virtual_time, data.get("exit_reason")), flush=True)
            return data
        except Exception as e:
            print("[exit] 调用失败(可能已结束): %r" % e, flush=True)
            return None

    def measure(self, pos, channel):
        """移动+检测。返回 (measure_result, svd_deg 或 None)。"""
        data = self._action("/measure", pos[0], pos[1], channel, prefix="measure")
        if data.get("accepted") is not True:
            raise FatalError("/measure 未被接受: %s" % data)
        self.pos = (float(pos[0]), float(pos[1]))
        self.channel = int(channel)
        self.virtual_time = data.get("virtual_time_s", self.virtual_time)
        return data.get("measure_result"), data.get("svd_deg")

    def clear(self, pos, channel):
        """移动+清除。返回是否成功。"""
        data = self._action("/clear", pos[0], pos[1], channel, prefix="clear")
        if data.get("accepted") is not True:
            raise FatalError("/clear 未被接受: %s" % data)
        self.pos = (float(pos[0]), float(pos[1]))
        self.virtual_time = data.get("virtual_time_s", self.virtual_time)
        return data.get("clear_result") == "success"

    def real_time_left(self):
        if self.remaining is None or self.real_start is None:
            return None
        return self.remaining - (time.time() - self.real_start)


# ---------------------------------------------------------------------------
# 策略: 定位并清除单个干扰源
# ---------------------------------------------------------------------------

def _probe_points(c, d):
    """探针清除布点(质心 + 单/双环), 覆盖圆盘(c, d) 内任一点到某探针 <= 20m。

    布局与证明(环半径 r, 相邻探针夹角半角 a):
    - d<=25:  中心 + 8 环@0.7d,  a=22.5°:  f(rho) 端点最大, f(d)=0.326d^2, f(20)<... <=400;
    - d<=50:  中心 + 12 环@0.7d, a=15°:   f(d)=0.1374d^2<=400, f(20)=400+0.49d^2-27.05d<=400 (d<=52);
    - d<=70:  中心 + 12 环@0.45d + 12 环@0.85d:
        外缘 f(d)=(1+0.7225-1.7cos15°)d^2=0.0803d^2<=400 (d<=70);
        内环带 f(20)=400+0.2025d^2-17.4d<=400; 环间最坏 0.0776d^2 (d<=65 宽松)。
    干扰源必在定位区域内, 故探针清除必然成功, 与定向源朝向无关(光学只认距离)。"""
    pts = [c]
    if d <= 25.0:
        for k in range(8):
            pts.append(vadd(c, vscale(unit(45.0 * k), 0.7 * d)))
    elif d <= 50.0:
        for k in range(12):
            pts.append(vadd(c, vscale(unit(30.0 * k), 0.7 * d)))
    else:
        for k in range(12):
            pts.append(vadd(c, vscale(unit(30.0 * k), 0.45 * d)))
        for k in range(12):
            pts.append(vadd(c, vscale(unit(30.0 * k + 15.0), 0.85 * d)))
    return pts


def strip_intersection_jump(dets):
    """从已有检测中选一对示向中心线求交点, 用于长轴未界定时直接跳测,
    省掉沿示向线走超再折返的时间。选夹角最大且交点在场地内的一对;
    无合适对(夹角<12°或交点出界)返回 None, 调用方回退沿线推进。"""
    best_ang = math.radians(12.0)
    best_pt = None
    for i in range(len(dets)):
        for j in range(i + 1, len(dets)):
            (p, b1), (q, b2) = dets[i], dets[j]
            a1, a2 = b1 % 180.0, b2 % 180.0
            ang = abs(a1 - a2)
            if ang > 90.0:
                ang = 180.0 - ang
            if ang < math.degrees(best_ang) - 1e-9:
                continue
            u1, u2 = unit(a1), unit(a2)
            den = u1[0] * u2[1] - u2[0] * u1[1]      # cross(u1, u2)
            if abs(den) < 1e-9:
                continue
            w = vsub(q, p)
            t = (w[0] * u2[1] - u2[0] * w[1]) / den
            pt = vadd(p, vscale(u1, t))
            if dist((0.0, 0.0), pt) > 1900.0:
                continue
            if ang >= math.degrees(best_ang) - 1e-9:
                best_ang = math.radians(ang)
                best_pt = pt
    return best_pt


def locate_and_clear(client, ch, dets, cleared, verbose=True):
    """对频道 ch 的干扰源执行 逼近-定位-清除。

    dets: [(检测点, 示向度), ...], 会被本函数追加新检测。
    返回 True 表示已清除。
    定向干扰源处理: 记录最近有信号点(必在覆盖扇形内), 一旦 no_signal
    即退回该点并转入小步爬行(朝区域质心步进, 扇形为凸集, 直线段保持在覆盖内),
    爬行步长逐次减半直至定位成功。
    """
    region = channel_region(dets)
    if not region:                       # 数值兜底(理论上不会发生)
        region = list(ARENA_GON)
    step_phase = 0                       # 0=远距离大步, 1=已掉头/近距离小步
    walk = None                          # [anchor, bearing, 已沿线前进距离]
    creep = False                        # 爬行模式(定向源边界保守移动)
    creep_step = 60.0
    creep_iters = 0                      # 爬行模式已用迭代数
    last_good = dets[-1][0]              # 最近有信号的检测点
    side = [1]                           # 侧向跳左右交替
    jump_once = []                       # 一次性交点跳测目标(非空时优先于沿线推进)
    force_walk = [False]                 # 质心测不到信号(定向源扇区外)时回退示向线推进

    def long_axis():
        n = len(region)
        pair = max(((region[i], region[j]) for i in range(n)
                    for j in range(i + 1, n)), key=lambda pq: dist(*pq),
                   default=((0.0, 0.0), (1.0, 0.0)))
        return bearing(pair[0], pair[1])

    def lateral_from(p, hop):
        """沿区域长轴垂直方向侧向跳(左右交替), 用于快速压扁细长定位区域。
        注意: 长轴须模 180° 定向, 否则顶点对顺序翻转会使方向恒定, 导致两点振荡。"""
        s = side[0]
        side[0] = -side[0]
        return vadd(p, vscale(unit(long_axis() % 180.0 + 90.0 * s), hop))

    def update_region(pos, svd):
        nonlocal region
        dets.append((pos, svd))
        new = clip_strip(region, pos, svd, BEARING_ERR)
        region = new if new else channel_region(dets)
        if not region:
            region = list(ARENA_GON)

    def retreat():
        """退回最近有信号点并转入/加速爬行模式。返回是否成功退回。"""
        nonlocal creep, creep_step, walk, region
        creep = True
        creep_step = max(12.0, creep_step / 2.0)
        walk = None
        res, svd = client.measure(last_good, ch)
        if res == "direction":
            update_region(last_good, svd)
            return True
        return res == "near" and client.clear(last_good, ch)

    def desperate_probe(center):
        """爬行失效后的兜底: 围绕 center 直接探针清除(与定向朝向无关)。
        探针也失败(典型: 定向源质心在扇区外)则回退示向线推进
        (扇形为凸集, 检测点到源的线段全程有信号, 必收敛)。返回是否已清除。"""
        nonlocal creep, creep_step, creep_iters, walk
        for p in _probe_points(center, 45.0):
            if client.clear(p, ch):
                cleared.add(ch)
                if verbose:
                    print("[ch%02d] 绝望探针清除成功 虚拟时刻 %.0f s"
                          % (ch, client.virtual_time), flush=True)
                return True
        creep = False
        creep_iters = 0
        walk = None
        force_walk[0] = True
        return False

    for it in range(80):
        d = polygon_diameter(region)
        c = polygon_centroid(region)
        cur = client.pos
        dc = dist(cur, c)

        # ---- 区域直径<=60m: 探针清除(网格必覆盖质心 60m 盘, 与定向朝向无关) ----
        # 普查阶段的跳过规则(区域<=60m 才停测)保证进入本函数时若区域已收敛
        # 必在此清除成功; 万一未中(数值边界兜底), 继续迭代细化而非放弃。
        if d <= 60.0:
            probes = _probe_points(c, d)
            # 各环从距当前位置最近的点开始, 省一次径向跳转
            for lo, hi in ((1, 13), (13, len(probes))):
                if hi - lo <= 2:
                    continue
                ring_pts = probes[lo:hi]
                ref = cur if lo == 1 else probes[lo - 1]
                k = min(range(len(ring_pts)), key=lambda i: dist(ring_pts[i], ref))
                probes[lo:hi] = ring_pts[k:] + ring_pts[:k]
            for p in probes:
                if client.clear(p, ch):
                    cleared.add(ch)
                    if verbose:
                        print("[ch%02d] 探针清除成功 @(%.0f,%.0f) 虚拟时刻 %.0f s"
                              % (ch, p[0], p[1], client.virtual_time), flush=True)
                    return True
            continue

        # ---- 爬行模式: 小步朝质心移动(全向源几乎不会进入) ----
        if creep:
            creep_iters += 1
            if creep_iters > 22:
                # 爬行超时: 机器人很可能已距源很近(爬行曾获近距离示向),
                # 围绕当前位置直接探针清除(与定向朝向无关), 不再恋战
                if desperate_probe(cur):
                    return True
                continue
            if dc > 1e-6 and dc > creep_step:
                step = min(creep_step, dc)
                nxt = vadd(cur, vscale(vsub(c, cur), step / dc))
            else:
                # 接近/到达质心但区域未收敛: 沿条带垂直方向侧移
                nxt = lateral_from(cur, creep_step)
            res, svd = client.measure(nxt, ch)
            if res == "direction":
                update_region(nxt, svd)
                last_good = nxt
            elif res == "near":
                if client.clear(nxt, ch):
                    cleared.add(ch)
                    return True
            else:
                if not retreat() and creep_step <= 12.0:
                    # retreat 失败且步长已减到最小: 绝望探针兜底
                    if desperate_probe(cur):
                        return True
            continue

        # ---- 区域长轴未界定(d>350 或质心不可达): 沿最近检测点的示向线推进 ----
        # 首次且已有>=2条示向线时, 直接跳测两示向中心线交点(估计误差
        # ~L*2°/sin夹角 << 1000m 接收半径), 省掉沿线走超再折返; 跳测失败
        # (定向源扇区外)则回退到原有沿线推进。
        if d > 350.0 or force_walk[0]:
            did_jump = False
            if walk is None:
                p0, b0 = min(dets, key=lambda pb: dist(pb[0], cur))
                walk = [p0, b0, 0.0]
                if len(dets) >= 2:
                    jp = strip_intersection_jump(dets)
                    if jp is not None:
                        jump_once.append(jp)
            step = 400.0 if step_phase == 0 else 200.0
            if jump_once:
                nxt = jump_once.pop()
                did_jump = True
            else:
                walk[2] += step
                nxt = vadd(walk[0], vscale(unit(walk[1]), walk[2]))
            res, svd = client.measure(nxt, ch)
            if res == "direction":
                force_walk[0] = False        # 已获得信号, 回到正常收敛
                if not did_jump and abs(ang_diff(svd, walk[1])) > 90.0:
                    step_phase = 1           # 示向翻转: 已走过干扰源
                    walk = [nxt, svd, 0.0]   # 掉头沿新示向线回走
                update_region(nxt, svd)
                last_good = nxt
            elif res == "near":
                if client.clear(nxt, ch):
                    cleared.add(ch)
                    return True
            elif did_jump:
                # 交点跳测无信号(定向源扇区外/估计出界): 回退沿线推进
                p0, b0 = min(dets, key=lambda pb: dist(pb[0], nxt))
                walk = [p0, b0, 0.0]
                step_phase = 0
            else:  # no_signal: 走出定向覆盖边界, 退回中点再退回安全点
                back = vadd(nxt, vscale(unit(walk[1]), -step / 2.0))
                res2, svd2 = client.measure(back, ch)
                if res2 == "direction":
                    update_region(back, svd2)
                    last_good = back
                    walk = [back, svd2, 0.0]
                    step_phase = 1
                elif res2 == "near":
                    if client.clear(back, ch):
                        cleared.add(ch)
                        return True
                elif not retreat():
                    return False
            continue

        # ---- 区域已界定(d<=350): 快速收敛 ----
        walk = None
        if dc > 300.0:
            nxt = c                     # 远离区域: 一次转运直达质心(顺带检测)
        else:
            # 垂直于区域长轴横向跳测(交会角≈90°, 一步把区域压到 ~2*D*sin1°)
            hop = min(300.0, max(60.0, d))
            nxt = lateral_from(cur, hop)
        res, svd = client.measure(nxt, ch)
        if res == "direction":
            update_region(nxt, svd)
            last_good = nxt
        elif res == "near":
            if client.clear(nxt, ch):
                cleared.add(ch)
                return True
        else:  # no_signal: 多半跳出定向覆盖边界, 退回爬行
            if not retreat():
                return False

    return False


# ---------------------------------------------------------------------------
# 策略主流程
# ---------------------------------------------------------------------------

def tsp_order(channels, detections, start):
    """第二阶段频道清除顺序: 最近邻构造 + 2-opt 改良(开放路径, 首端固定)。
    以各频道定位区域质心为近似位置, 清除后机器人实际位于源附近,
    与质心近似误差(<=区域半径)相对频道间距离可忽略。"""
    pts = {c: polygon_centroid(channel_region(detections[c]))
           for c in channels}
    order = []
    cur = start
    pool = list(channels)
    while pool:
        nxt = min(pool, key=lambda c: dist(cur, pts[c]))
        order.append(nxt)
        cur = pts[nxt]
        pool.remove(nxt)
    improved = True
    while improved:
        improved = False
        for i in range(1, len(order) - 1):
            for j in range(i + 1, len(order)):
                a, b, cc = pts[order[i - 1]], pts[order[i]], pts[order[j]]
                e = pts[order[j + 1]] if j + 1 < len(order) else None
                old = dist(a, b) + (dist(cc, e) if e else 0.0)
                new = dist(a, cc) + (dist(b, e) if e else 0.0)
                if new < old - 1e-9:
                    order[i:j + 1] = order[i:j + 1][::-1]
                    improved = True
    return order


def run_strategy(client, problem=3, verbose=True, min_real_left=45.0):
    """完整执行一次测试。返回 (清除数, 检测到的频道集合)。"""
    client.enter()
    detections = {ch: [] for ch in range(CH_MIN, CH_MAX + 1)}
    cleared = set()

    def log(msg):
        if verbose:
            print(msg, flush=True)

    # ---- 第 1 步: 普查(扫描 + 初定位合一) ----
    census = census_points(problem)
    reordered = False
    for i, p in enumerate(census):
        # 自适应重排剩余普查点(只做一次, 过半后): 此时已发现源群的质心
        # 已知, 让巡回收束在源群附近, 缩短第二阶段首段转运。覆盖不受顺序
        # 影响(每个普查点都必须走到), 仅多点贪心开销 ~0.35 系数折中。
        if not reordered and i + 3 >= len(census) // 2 and len(census) - i - 1 >= 3:
            found = [ch for ch in detections if detections[ch]]
            if len(found) >= 3:
                tgt = (sum(polygon_centroid(channel_region(detections[c]))[0] for c in found) / len(found),
                       sum(polygon_centroid(channel_region(detections[c]))[1] for c in found) / len(found))
                rest = census[i + 1:]
                ordered = []
                cur = p
                while rest:
                    nxt = min(rest, key=lambda q: dist(cur, q) + 0.35 * dist(q, tgt))
                    ordered.append(nxt)
                    cur = nxt
                    rest.remove(nxt)
                census[i + 1:] = ordered
                reordered = True
        for ch in range(CH_MIN, CH_MAX + 1):
            if ch in cleared:
                continue
            # 已有两次及以上检测且定位区域直径已收敛(<=60m, 探针可保证
            # 覆盖): 该频道交给第 2 步处理, 跳过重复测量(省 6 s/次)
            if len(detections[ch]) >= 2 and \
                    polygon_diameter(channel_region(detections[ch])) \
                    <= CENSUS_REGION_D:
                continue
            res, svd = client.measure(p, ch)
            if res == "direction":
                detections[ch].append((p, svd))
            elif res == "near":
                # 距源 <=5m: 直接精确定位并清除
                if client.clear(p, ch):
                    cleared.add(ch)
                    log("[census] ch%02d near->清除成功 @(%.0f,%.0f)"
                        % (ch, p[0], p[1]))
        log("[census] %d/%d 点 (%.0f,%.0f) 虚拟时刻 %.0f s 已发现 %d 频道"
            % (i + 1, len(census), p[0], p[1],
               client.virtual_time,
               sum(1 for v in detections.values() if v)))

    # ---- 第 2 步: 逐个定位清除(2-opt 全局顺序) ----
    planned = []                         # 待清除频道顺序(tsp_order 计算)
    while True:
        todo = [ch for ch in detections if detections[ch] and ch not in cleared]
        if not todo:
            break
        left = client.real_time_left()
        if left is not None and left < min_real_left:
            log("[warn] 现实时间不足(%.0f s), 提前退出" % left)
            break
        # 静态计划一次到底: 每步全量重排(2-opt)反而更差——质心近似误差
        # 会误导后续链; 一次性 2-opt + 逐段执行最优
        while planned and planned[0] not in todo:
            planned.pop(0)
        if not planned:
            planned = tsp_order(todo, detections, client.pos)
        ch = planned.pop(0)
        before = len(cleared)
        ok = locate_and_clear(client, ch, detections[ch], cleared, verbose)
        log("[ch%02d] %s (已清除 %d, 虚拟时刻 %.0f s)"
            % (ch, "OK" if ok else "FAILED", len(cleared), client.virtual_time))
        if len(cleared) == before and not ok:
            log("[warn] 频道 %d 定位失败, 跳过" % ch)
            detections[ch] = []     # 放弃该频道, 避免死循环(理论上不会发生)

    client.exit()
    log("[done] 清除 %d 个干扰源, 退出时虚拟时刻 %.0f s"
        % (len(cleared), client.virtual_time))
    return len(cleared), {c for c, v in detections.items() if v}


# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description="CUMCM2026 B 题机器狗程序")
    ap.add_argument("--robot", required=True, help="参赛队号(即 robot_id)")
    ap.add_argument("--problem", type=int, choices=[3, 4], default=3)
    ap.add_argument("--url", default="http://127.0.0.1:2026")
    args = ap.parse_args()

    log_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "logs")
    os.makedirs(log_dir, exist_ok=True)
    log_path = os.path.join(log_dir, "run-%s.jsonl" % datetime.now().strftime("%Y%m%d-%H%M%S"))

    client = Client(args.url, args.robot, log_path)
    print("日志文件: %s" % log_path, flush=True)
    try:
        run_strategy(client, problem=args.problem)
    except FatalError as e:
        print("致命错误: %s" % e, flush=True)
        client.exit()
        sys.exit(1)


if __name__ == "__main__":
    main()
