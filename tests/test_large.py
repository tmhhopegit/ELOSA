"""Multi-word bitsets (n > 64), non-grouped patients and thread hand-offs vs brute force."""
import os, subprocess, sys, tempfile
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from elosa_io import write_input, read_groups, brute
from synth import synth_lesions

EXE = os.environ.get('ELOSA_EXE', os.path.join(os.path.dirname(__file__), '..', 'elosa'))
rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
tmp = tempfile.mkdtemp(); fails = 0; trials = int(sys.argv[2]) if len(sys.argv) > 2 else 12
for trial in range(trials):
    n = int(rng.integers(65, 200))
    L = synth_lesions(n, shape=(14, 16, 14), seed=int(rng.integers(1e9)), mean_radius=1.6)
    use = rng.random(n) < 0.6 if trial % 2 else None
    fin = os.path.join(tmp, 'in.bin'); fout = os.path.join(tmp, 'g')
    write_input(fin, L, use)
    exp = brute(L, use)
    for th in (1, 3):
        r = subprocess.run([EXE, fin, '--out', fout, '--covering', '--threads', str(th), '--quiet'], capture_output=True, text=True)
        got = read_groups(fout)
        gs = {g: (s, o, c) for g, s, o, c in got}
        ok = r.returncode == 0 and len(gs) == len(got) and gs == exp
        print(f'n={n} use={"subset" if use is not None else "all"} threads={th} groups={len(exp)} {"ok" if ok else "MISMATCH"}', flush=True)
        fails += not ok
print('failures', fails)
