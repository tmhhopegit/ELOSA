import os, subprocess, sys, time
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from synth import synth_lesions
from elosa_io import write_input
EXE = os.path.join(os.path.dirname(__file__), '..', 'elosa')
for n in [int(x) for x in sys.argv[1].split(',')]:
    L = synth_lesions(n, seed=n)
    write_input('/tmp/b.bin', L)
    for th in [1, int(sys.argv[2]) if len(sys.argv) > 2 else 4]:
        t = time.time()
        r = subprocess.run([EXE, '/tmp/b.bin', '--threads', str(th), '--quiet'], capture_output=True, text=True)
        print(f'n={n} voxels={L.shape[1]} lesioned={L.any(0).sum()} threads={th} groups={r.stdout.strip()} time={time.time()-t:.2f}s', flush=True)
