"""Synthetic stroke-like lesions for benchmarks: 1-3 ellipsoid blobs per patient,
mostly in one hemisphere, clustered around a few 'vascular territory' centres."""
import numpy as np


def synth_lesions(n, shape=(46, 55, 46), seed=0, mean_radius=4.0):
    rng = np.random.default_rng(seed)
    X, Y, Z = np.meshgrid(*[np.arange(s) for s in shape], indexing='ij')
    c = np.array(shape) / 2
    brain = ((X - c[0]) / (shape[0] / 2.2)) ** 2 + ((Y - c[1]) / (shape[1] / 2.2)) ** 2 + ((Z - c[2]) / (shape[2] / 2.2)) ** 2 <= 1
    idx = np.flatnonzero(brain)
    pts = np.stack([X.ravel()[idx], Y.ravel()[idx], Z.ravel()[idx]], 1).astype(float)
    centres = np.array([[shape[0] * 0.3, shape[1] * 0.5, shape[2] * 0.55],
                        [shape[0] * 0.35, shape[1] * 0.35, shape[2] * 0.45],
                        [shape[0] * 0.3, shape[1] * 0.65, shape[2] * 0.6],
                        [shape[0] * 0.4, shape[1] * 0.5, shape[2] * 0.4]])
    L = np.zeros((n, len(idx)), bool)
    for p in range(n):
        for _ in range(rng.integers(1, 4)):
            ctr = centres[rng.integers(len(centres))] + rng.normal(0, 4, 3)
            rad = np.abs(rng.normal(mean_radius, mean_radius / 2, 3)) + 1
            L[p] |= (((pts - ctr) / rad) ** 2).sum(1) <= 1
    return L
