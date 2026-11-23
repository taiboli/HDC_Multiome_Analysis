"""
Drop-in replacement for spire.omp.omp (same inputs, same output) that is faster for many TFs.

spire.omp.omp updates the residual correlations of every peak as alpha0 - G[:, S] beta, which gathers a
(TFs x peaks x |S|) block every step. Here the residual is formed in pseudocell space, R = Y - D_S beta
(gather of n x peaks x |S|, n = pseudocells << TFs), and the correlations come from one BLAS product D^T R.
Selection, least-squares refit (normal equations on the Gram sub-matrix) and stopping rules are unchanged, so the
selected supports and coefficients are identical up to floating-point rounding (checked in __main__).
"""
import os
import numpy as np
import scipy.sparse as sp


def omp_fast(D, Y, L=10, eps=0.0, chunk_size=4096, ridge=1e-10):
    D = np.asarray(D, dtype=np.float64)
    Y = np.asarray(Y, dtype=np.float64)
    n, T = D.shape
    if Y.ndim == 1:
        Y = Y[:, None]
    P = Y.shape[1]
    L = int(min(L, T, n))
    G = D.T @ D
    rows, cols, vals = [], [], []
    for start in range(0, P, chunk_size):
        stop = min(start + chunk_size, P)
        Yc = Y[:, start:stop]
        m = stop - start
        alpha0 = np.ascontiguousarray((D.T @ Yc).T)  # peaks x TFs (row-major: contiguous per-peak gathers)
        ynorm2 = np.einsum("ij,ij->j", Yc, Yc)
        corr = alpha0.copy()
        support = np.full((m, L), -1, dtype=np.int64)
        coef = np.zeros((m, L))
        active = ynorm2 > eps
        k_used = np.zeros(m, dtype=np.int64)
        for k in range(L):
            if not active.any():
                break
            idx = np.where(active)[0]
            full = len(idx) == m
            c = np.abs(corr) if full else np.abs(corr[idx])
            if k > 0:
                c[np.arange(len(idx))[:, None], support[idx, :k]] = -1.0
            new = np.argmax(c, axis=1)
            ok = c[np.arange(len(idx)), new] > 1e-12
            idx, new = idx[ok], new[ok]
            if idx.size == 0:
                break
            support[idx, k] = new
            S = support[idx, : k + 1]
            Gs = G[S[:, :, None], S[:, None, :]] + ridge * np.eye(k + 1)[None]
            b = alpha0[idx[:, None], S]
            beta = np.linalg.solve(Gs, b[..., None])[..., 0]
            coef[idx, : k + 1] = beta
            k_used[idx] = k + 1
            R = Yc[:, idx] - np.einsum("nak,ak->na", D[:, S], beta)  # residuals in pseudocell space
            if len(idx) == m:
                np.matmul(R.T, D, out=corr)
            else:
                corr[idx] = R.T @ D
            res2 = ynorm2[idx] - np.einsum("ak,ak->a", beta, b)
            still = np.zeros(m, dtype=bool)
            still[idx] = res2 > eps
            active = active & still
        kk = k_used
        if kk.any():
            j = np.repeat(np.arange(m), kk)
            pos = np.concatenate([np.arange(x) for x in kk])
            rows.append(support[j, pos])
            cols.append(start + j)
            vals.append(coef[j, pos])
    if rows:
        r, c, v = np.concatenate(rows), np.concatenate(cols), np.concatenate(vals)
    else:
        r = c = np.array([], dtype=np.int64)
        v = np.array([], dtype=np.float64)
    return sp.csc_matrix((v, (r, c)), shape=(T, P))


if __name__ == "__main__":
    import sys
    import time

    sys.path.insert(0, os.environ["SPIRE_SRC"])
    from spire.omp import normalize_dictionary, omp

    rng = np.random.default_rng(0)
    for n, T, P in ((300, 1300, 8192), (60, 40, 500)):
        X = rng.normal(size=(n, T)) + rng.normal(size=(n, 1))  # shared component: collinear TFs
        D, _ = normalize_dictionary(X)
        B = sp.random(T, P, density=5 / T, random_state=1).toarray()
        Y = D @ B + 0.1 * rng.normal(size=(n, P))
        Y -= Y.mean(0)
        t = time.time(); A = omp(D, Y, L=10); t1 = time.time() - t
        t = time.time(); F = omp_fast(D, Y, L=10); t2 = time.time() - t
        same = (A != 0).astype(int) - (F != 0).astype(int)
        print(f"n={n} T={T} P={P}: spire {t1:.2f}s fast {t2:.2f}s | support differences {abs(same).sum()} | "
              f"max |coef diff| {abs(A - F).max():.2e}")
