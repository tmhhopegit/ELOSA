import os, subprocess, sys, tempfile
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from elosa_io import write_input, read_groups, brute

EXE = os.environ.get('ELOSA_EXE', os.path.join(os.path.dirname(__file__), '..', 'elosa'))
rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
trials = int(sys.argv[2]) if len(sys.argv) > 2 else 300
tmp = tempfile.mkdtemp()
fails = 0
for trial in range(trials):
    n = int(rng.integers(1, 12)); V = int(rng.integers(1, 80))
    L = rng.random((n, V)) < rng.uniform(0.05, 0.8)
    use = rng.random(n) < 0.7 if rng.random() < 0.4 else None
    if use is not None and not use.any():
        use[0] = True
    minv = int(rng.choice([1, 1, 2, 4])); maxp = int(rng.choice([0, 0, 2, 3]))
    order = str(rng.choice(['ascending', 'descending', 'input']))
    threads = int(rng.choice([1, 2, 4]))
    fin = os.path.join(tmp, 'in.bin'); fout = os.path.join(tmp, 'g')
    write_input(fin, L, use)
    r = subprocess.run([EXE, fin, '--out', fout, '--covering', '--threads', str(threads), '--min-voxels', str(minv),
                        '--max-patients', str(maxp), '--order', order, '--quiet'], capture_output=True, text=True)
    if r.returncode != 0:
        if use is not None or L.any():
            print('error', r.stderr); fails += 1
        continue
    got = read_groups(fout)
    gs = {g: (s, o, c) for g, s, o, c in got}
    exp = brute(L, use, minv, maxp)
    ok = len(gs) == len(got) and gs == exp and int(r.stdout) == len(exp)
    if not ok:
        fails += 1
        if fails < 4:
            print('MISMATCH trial', trial, n, V, minv, maxp, order, threads, 'got', len(got), len(gs), 'exp', len(exp))
            print(' missing', list(set(exp) - set(gs))[:5], 'extra', list(set(gs) - set(exp))[:5])
            print(' diffs', [(g, gs[g], exp[g]) for g in gs if g in exp and gs[g] != exp[g]][:5])
print('trials', trials, 'failures', fails)
