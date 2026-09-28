"""Read/write ELOSA binary files (mirrors the MATLAB functions; used by the tests)."""
import numpy as np


def write_input(path, lesions, use=None):
    L = np.asarray(lesions, dtype=bool)
    n, V = L.shape
    use = np.ones(n, bool) if use is None else np.asarray(use, bool)
    with open(path, 'wb') as f:
        f.write(b'ELOSAIN1')
        f.write(np.array([n, V], '<u4').tobytes())
        f.write(use.astype(np.uint8).tobytes())
        f.write(np.ascontiguousarray(L.T).astype(np.uint8).tobytes())


def read_groups(prefix):
    """Returns a list of (group, region_voxels, n_other, covering_others or None)."""
    raw = open(prefix + '.groups', 'rb').read()
    assert raw[:8] == b'ELOSAGR2'
    n, flags = np.frombuffer(raw[8:16], '<u4')
    h = np.frombuffer(raw[16:], '<u4').reshape(-1, 4).astype(np.int64)
    mem = np.fromfile(prefix + '.members', '<u2' if flags & 2 else '<u4').astype(np.int64) - 1
    out, i = [], 0
    for k, other, lo, hi in h:
        g = frozenset(int(x) for x in mem[i:i + k]); i += k
        cov = None
        if flags & 1:
            cov = frozenset(int(x) for x in mem[i:i + other]); i += other
        out.append((g, int(lo) | (int(hi) << 32), int(other), cov))
    assert i == len(mem)
    return out


def brute(L, use=None, min_voxels=1, max_patients=0):
    """All groups by definition: closures of every patient set with a non-empty region."""
    L = np.asarray(L, bool)
    n, V = L.shape
    use = np.ones(n, bool) if use is None else np.asarray(use, bool)
    cols = {}
    for v in range(V):
        c = L[:, v]
        if (c & use).any():
            key = tuple(np.flatnonzero(c))
            cols[key] = cols.get(key, 0) + 1
    atoms = [(frozenset(k), w) for k, w in cols.items()]
    U = frozenset(np.flatnonzero(use))
    base = {a & U for a, _ in atoms}
    fam = set(base)
    frontier = set(base)
    while frontier:
        new = set()
        for s in frontier:
            for b in base:
                t = s & b
                if t and t not in fam:
                    new.add(t)
        fam |= new
        frontier = new
    res = {}
    for g in fam:
        if not g:
            continue
        inside = [(a, w) for a, w in atoms if g <= a]
        sup = sum(w for _, w in inside)
        cover = frozenset.intersection(*[a for a, _ in inside])
        assert cover & U == g
        if sup >= min_voxels and (max_patients == 0 or len(g) <= max_patients):
            res[g] = (sup, len(cover - U), frozenset(int(x) for x in cover - U))
    return res
