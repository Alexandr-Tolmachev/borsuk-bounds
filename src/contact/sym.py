"""Sign symmetries of the two-sided frame coordinates."""
import numpy as np
GROUPS = {'flip3': [(), (3,)], 'flip34': [(), (3,), (4,), (3, 4)]}

def signs(name, n):
    out = []
    # only the two-sided frame coordinates 3..n-2 may be flipped (x_{n-1} <= 1/2 is one-sided)
    assert all(3 <= j <= n - 2 for flip in GROUPS[name] for j in flip), 'flip of a one-sided coordinate'
    for flip in GROUPS[name]:
        s = np.ones(n); s[list(flip)] = -1; out.append(s)
    return out
