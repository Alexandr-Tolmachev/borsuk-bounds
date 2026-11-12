"""Floating-point ball model for proposing tangent halfspaces."""
import numpy as np

class Ball:
    def __init__(self, c, R): self.c = np.asarray(c, float); self.R = float(R)
    def inside(self, X, tol): return np.linalg.norm(X - self.c, axis=1) <= self.R + tol
