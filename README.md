# CUMCM 2026 B 题 机器狗自动定位清除干扰源程序

## 程序文件

- `robot.py` — 支撑官方模拟器的算法文件
- `selftest.py` — 离线自检
- `logs` — 存放运行日志

## 正式运行

前置条件：确保模拟器已正常启动，在模拟器中选定测试类型。

倒计时结束、界面显示接口就绪后：

```powershell
python robot.py --robot <参赛队号> --problem 3
python robot.py --robot <参赛队号> --problem 4
```

- `--robot`：参赛队号
- `--problem`：3 或 4
- `--url`：模拟器地址，默认 `http://127.0.0.1:2026`

程序在倒计时期间会自动重试 `/enter`；每次运行生成日志 `logs/run-<时间戳>.jsonl`。

## 离线自检

```bash
python selftest.py --problem 3 --trials 300     
python selftest.py --problem 4 --trials 300
python selftest.py --coverage 300000            
```


## 策略说明

策略原理、几何覆盖证明与优化细节请见论文正文。
