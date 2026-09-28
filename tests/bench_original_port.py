"""Time the NumPy port of the original search: python3 bench_original_port.py 50,100"""
import os, sys, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from synth import synth_lesions
from original_port import original, atoms_of
for n in [int(x) for x in sys.argv[1].split(',')]:
    L = synth_lesions(n, seed=n)
    t = time.time()
    A, c = atoms_of(L)
    t1 = time.time()
    out = original(A)
    t2 = time.time()
    print(f'n={n} atoms={len(A)} groups_found={len(set(out))} unique_rows_time={t1-t:.2f}s search={t2-t1:.2f}s per_group={1e6*(t2-t1)/max(1,len(out)):.0f}us', flush=True)
