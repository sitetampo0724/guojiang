"""End-to-end check of the delivered command line against the local mock server.

Runs `robot_dog.main()` exactly as the official run would (HTTP interface,
output files, exit code) against the self-built simulator's HTTP server.

Usage:
    python3 lab/cli_check.py --problem 3 --scenario boundary --sources 13 --seed 2026092904
"""

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import robot_dog  # noqa: E402
import simulator  # noqa: E402


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('--problem', type=int, default=3)
    parser.add_argument('--scenario', default='boundary')
    parser.add_argument('--sources', type=int, default=13)
    parser.add_argument('--seed', type=int, default=2026092904)
    parser.add_argument('--port', type=int, default=2137)
    parser.add_argument('--output-dir', default='/tmp/q3_ultimate_cli')
    args = parser.parse_args(argv)

    if args.problem == 3:
        sources = simulator.sources_for(args.seed, args.scenario, args.sources)
        transport = simulator.MockTransport(seed=args.seed, sources=sources, robot_id='local-test')
    else:
        sources = simulator.sources_for_q4(args.seed, args.scenario)
        transport = simulator.MockTransport4(seed=args.seed, sources=sources, robot_id='local-test')

    class DirectTransport(robot_dog.HTTPTransport):
        """Same interface as HTTPTransport, wired straight to the mock server.

        The sandbox blocks TCP binding, so this exercises the whole CLI path
        (argument parsing, output files, exit code) without a real socket.
        """

        def __init__(self, base_url='http://127.0.0.1:2026', timeout_s=5):
            self.base_url = base_url
            self.timeout_s = timeout_s

        def post(self, path, payload):
            return transport.post(path, payload)

    robot_dog.HTTPTransport = DirectTransport
    code = robot_dog.main(['--problem', str(args.problem), '--robot-id', 'local-test',
                           '--base-url', f'http://127.0.0.1:{args.port}',
                           '--output-dir', args.output_dir, '--quiet'])
    report = transport.report()
    out_dir = Path(args.output_dir)
    produced = sorted(p.name for p in out_dir.glob('*'))
    summary = json.loads((out_dir / produced[[i for i, n in enumerate(produced) if n.endswith('.json')][0]]).read_text())
    print(json.dumps(dict(exit_code=code, status=summary['status'],
                          cleared=summary['cleared_count'], total=report['source_count'],
                          virtual_s=round(summary['virtual_time_s'], 1),
                          files=produced,
                          decisions=summary.get('dp_decision_count'),
                          certificate=summary['complete_certificate']), ensure_ascii=False, indent=2))
    return 0 if code == 0 else 1


if __name__ == '__main__':
    sys.exit(main())
