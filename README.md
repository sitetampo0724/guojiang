# CUMCM 2026 B 题 机器狗自动定位清除程序

## 文件

- `robot.py` — 正式比赛程序：HTTP 客户端（幂等重试、JSONL 日志）+ 几何库 + 定位清除策略
- `selftest.py` — 离线蒙特卡洛自检（不需要模拟器）

## 正式运行

前置：模拟器已启动、已联网登录，并在模拟器中选定测试类型（问题3/问题4 演练或正式）点击开始。

倒计时结束、界面显示接口就绪后：

```bash
python robot.py --robot <参赛队号> --problem 3
python robot.py --robot <参赛队号> --problem 4
```

- `--robot`：参赛队号（即 robot_id，需与登录队号一致）
- `--problem`：3 或 4
- `--url`：模拟器地址，默认 `http://127.0.0.1:2026`

程序在倒计时期间会自动重试 `/enter`；每次运行生成日志 `logs/run-<时间戳>.jsonl`。

## 离线自检（无需模拟器）

```bash
python selftest.py --problem 3 --trials 300     # 蒙特卡洛 300 局
python selftest.py --problem 4 --trials 300
python selftest.py --coverage 300000            # 普查几何覆盖验证
```

当前测试结果（300 局 × 2 种子，均 **100% 完全清除**）：
Q3 平均虚拟总时间 ≈ 4220 s（~330 s/源），Q4 ≈ 8350 s（~655 s/源）；
普查覆盖 30 万采样零漏检。

## 策略说明

策略原理、几何覆盖证明与优化细节见论文正文；代码中关键设计
（普查布局参数与证明、探针覆盖证明、凸性兜底论证）在各函数 docstring 内。
