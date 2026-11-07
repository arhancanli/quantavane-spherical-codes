"""The polish's linear programme at a converged code, and its dual multipliers compared with the equilibrium stress.

  python3 -I certify/lp_dual.py CODE.txt D N [--eps 1e-8]

Builds the linear programme of the sequential LP polish exactly as engine/spherical/sphopt.py writes it (variables:
tangent moves D in [-1, 1]^(N D) and t; for every pair with <x_i, x_j> > c - 4 eps the row
<x_i, D_j> + <x_j, D_i> - t <= (c - <x_i, x_j>) / eps; tangency <x_i, D_i> = 0; minimise t), solves it with HiGHS's
interior-point method as the polish does, and compares the multipliers of the pair rows (y = -marginals) with the
equilibrium stress of the contacts computed separately (least squares on the contacts within 1e-10, normalised to sum
1). Prints the optimal t, the sum of y, the largest |y - stress| relative to the largest stress, and the weight of y
outside the contacts. Not part of the certificate; it illustrates Proposition 4.1 of the paper.
"""
import os, sys
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import sphcert as S


def lp_as_in_sphopt(X, eps):
    n, d = X.shape
    G = X @ X.T
    np.fill_diagonal(G, -2.0)
    c = G.max()
    I, J = np.nonzero(np.triu(G > c - 4 * eps, 1))
    k = len(I)
    nv = n * d + 1
    rows = np.repeat(np.arange(k), 2 * d + 1)
    cols = np.concatenate([np.stack([J[:, None] * d + np.arange(d), I[:, None] * d + np.arange(d)], 1).reshape(k, 2 * d),
                           np.full((k, 1), n * d)], 1).ravel()
    vals = np.concatenate([X[I], X[J], -np.ones((k, 1))], 1).ravel()
    A_ub = coo_matrix((vals, (rows, cols)), shape=(k, nv)).tocsr()
    b_ub = (c - G[I, J]) / eps
    A_eq = coo_matrix((X.ravel(), (np.repeat(np.arange(n), d), np.arange(n * d))), shape=(n, nv)).tocsr()
    cost = np.zeros(nv)
    cost[-1] = 1.0
    bounds = [(-1.0, 1.0)] * (n * d) + [(None, None)]
    res = linprog(cost, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=np.zeros(n), bounds=bounds, method="highs-ipm")
    return res, list(zip(I.tolist(), J.tolist())), c


def stress_lstsq(X, E):
    """Least-squares stress on the pairs E: tangential forces sum to zero at every point, total weight 1."""
    n, d = X.shape
    A = np.zeros((n * d + 1, len(E)))
    for e, (i, j) in enumerate(E):
        A[i * d:(i + 1) * d, e] = X[j] - (X[j] @ X[i]) * X[i]
        A[j * d:(j + 1) * d, e] = X[i] - (X[i] @ X[j]) * X[j]
    A[-1] = 1.0
    b = np.zeros(n * d + 1)
    b[-1] = 1.0
    w, *_ = np.linalg.lstsq(A, b, rcond=None)
    return w


def main():
    path, d, n = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    eps = float(sys.argv[sys.argv.index("--eps") + 1]) if "--eps" in sys.argv else 1e-8
    X = S.read_code(path, d, n)
    X = X / np.linalg.norm(X, axis=1, keepdims=True)
    res, pairs, c = lp_as_in_sphopt(X, eps)
    y = -res.ineqlin.marginals
    G = X @ X.T
    contacts = [(i, j) for i, j in pairs if G[i, j] >= c - 1e-10]
    w = stress_lstsq(X, contacts)
    idx = {p: k for k, p in enumerate(pairs)}
    yc = np.array([y[idx[p]] for p in contacts])
    off = sum(y[k] for k, p in enumerate(pairs) if G[p] < c - 1e-10)
    print(f"{os.path.basename(path)}: LP rows {len(pairs)}, contacts {len(contacts)}, optimal t {res.x[-1]:.2e}, "
          f"sum y {y.sum():.12f}, min y on contacts {yc.min():.3e}, weight of y off the contacts {off:.2e}, "
          f"max |y - stress| / max stress {np.abs(yc - w).max() / np.abs(w).max():.2e}, min stress {w.min():.3e}")


if __name__ == "__main__":
    main()
