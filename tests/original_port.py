"""NumPy port of the original ExhaustiveRegionFinder (CalcGroupsByVoxel + FindRecursive +
GrowGroup_v3), used to check what it finds and as a rough speed reference.
Run: python3 original_port.py [seed]  - compares it with the definition on random data."""
import numpy as np, itertools, sys
sys.setrecursionlimit(100000)

def atoms_of(L):
    # L: n x V bool; unique columns -> rows of atoms (A x n), with counts
    cols, counts = np.unique(L.T, axis=0, return_counts=True)
    return cols.astype(bool), counts

def brute_closed(atoms):
    # intersection closure of non-empty atoms, keep non-empty sets
    fam = set()
    base = [frozenset(np.flatnonzero(a)) for a in atoms if a.any()]
    frontier = set(base)
    fam |= frontier
    while frontier:
        new = set()
        for s in frontier:
            for b in base:
                t = s & b
                if t and t not in fam:
                    new.add(t)
        fam |= new
        frontier = new
    return fam

def original(atoms):
    # port of FindRecursive + GrowGroup_v3 (0-based indices)
    D = atoms  # rows incl zero row
    n = D.shape[1]
    out = []
    def grow(G, last, KG, R):
        out.append(frozenset(np.flatnonzero(G)))
        if len(R) == 0:
            return
        uKG, y = np.unique(KG.T, axis=0, return_inverse=True)
        y = np.asarray(y).ravel()
        after = R > last
        pos = np.arange(len(R))
        before_pats = set(y[~after])
        cand = [(y[i], i) for i in pos[after] if y[i] not in before_pats]
        seen = set(); kc = []
        for pat, i in cand:
            if pat not in seen:
                seen.add(pat); kc.append((pat, i))
        for pat, i in kc:
            rows = uKG[pat].astype(bool)
            KeyGp = KG[rows]
            AllIn = KeyGp.all(axis=0)
            if not AllIn[:i].any():
                sums = KeyGp.sum(axis=0)
                T = G.copy(); T[R] = AllIn
                remove = AllIn | (sums <= 1)
                grow(T, R[i], KeyGp[:, ~remove], R[~remove])
    for p in range(n):
        KeyGroups = D[D[:, p]]
        if KeyGroups.shape[0] > 1:
            Gp = KeyGroups.all(axis=0)
            if np.flatnonzero(Gp)[0] >= p:
                KG = KeyGroups[:, ~Gp]; KR = np.flatnonzero(~Gp)
                grow(Gp, p, KG, KR)
    return out

if __name__ == '__main__':
    rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
    stats = []
    for trial in range(300):
        n = rng.integers(3, 9); V = rng.integers(5, 60)
        L = rng.random((n, V)) < rng.uniform(0.2, 0.7)
        A, c = atoms_of(L)
        truth = brute_closed(A)
        o = original(A)
        so = set(o)
        missing = truth - so; extra = so - truth
        stats.append((len(truth), len(so), len(missing), len(extra), len(o) - len(so)))
    s = np.array(stats)
    print('trials', len(s), 'with missing', (s[:,2]>0).sum(), 'with extra', (s[:,3]>0).sum(),
          'total truth', s[:,0].sum(), 'total missing', s[:,2].sum(), 'dups', s[:,4].sum())
