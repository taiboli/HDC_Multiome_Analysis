"""
NumPy port of scMultiMap (Su et al. 2025 Nat Commun; github.com/ChangSuBiostats/scMultiMap v1.0.1).

Same estimators as R/scMultiMap_IRLS.R (step 1: moment estimates + IRLS for per-sample means, variance sigma^2
and a common NB overdispersion per biological sample) and R/scMultiMap_WLS*.R (step 2: WLS covariance estimate,
z-test, correlation = covar / sqrt(sigma2_gene * sigma2_peak)). Differences are implementation only:
  - step 1 runs over features in chunks (the per-sample median theta couples all features, so every IRLS
    iteration is one pass over all chunks), so the full dense cell x feature matrix is never held in memory;
  - step 2 evaluates blocks of genes against all peaks of the block's window with matrix products, using
    g = w_gene * W_peak  =>  sum_c x_c X_cp / g_cp * s_c = (x s / w_gene) . (X_p / W_p), etc.
Validated against the R functions in 23v_validate_scmultimap.R.
"""
import numpy as np
import scipy.sparse as sp
from scipy.stats import norm


def _dense(Xcsc, j0, j1):
    return Xcsc[:, j0:j1].toarray().astype(np.float64)


def irls(Xcsc, sd, bcode, n_b, irls=True, chunk=2000, max_iter=10, verbose=False):
    """Xcsc: cells x features (CSC); sd: seq depth per cell; bcode: int sample code per cell (0..n_b-1).
    Returns dict(mu=F x B, sigma_sq=F, theta=B (per-sample median theta))."""
    n, F = Xcsc.shape
    B = np.zeros((n, n_b))
    B[np.arange(n), bcode] = 1
    sd = sd.astype(np.float64)
    sd2 = sd ** 2
    seq_4 = np.sum(sd ** 4)
    chunks = [(j, min(j + chunk, F)) for j in range(0, F, chunk)]

    # moment estimates
    mu = np.empty((F, n_b))
    sigma2 = np.empty(F)
    den0 = B.T @ sd2  # B
    for j0, j1 in chunks:
        X = _dense(Xcsc, j0, j1)
        m = (B.T @ (X * sd[:, None])).T / den0[None, :]
        M = sd[:, None] * m[:, bcode].T
        sigma2[j0:j1] = (((X - M) ** 2 - M) * sd2[:, None]).sum(0) / seq_4
        mu[j0:j1] = m
    if not irls:
        return dict(mu=mu, sigma_sq=sigma2, theta=None)

    with np.errstate(divide="ignore", invalid="ignore"):
        theta = mu ** 2 / sigma2[:, None]
    j, delta = 0, np.inf

    def med(th):
        out = np.empty(n_b)
        for b in range(n_b):
            t = th[:, b]
            t = t[(t > 0) & ~np.isnan(t)]
            out[b] = np.median(t)
        return out

    while delta > 0.05 and j <= max_iter:
        theta_prev = theta.copy()
        tmed = med(theta)
        tc = tmed[bcode]  # per cell
        mu_new = np.empty_like(mu)
        for j0, j1 in chunks:
            X = _dense(Xcsc, j0, j1)
            M = sd[:, None] * mu[j0:j1, bcode].T
            w = M + M ** 2 / tc[:, None]
            w[np.isnan(w) | (w <= 0)] = 1
            num = B.T @ (X / w * sd[:, None])                   # B x f
            den = (sd2[:, None] / w).T @ B                       # f x B
            mu_new[j0:j1] = num.T / den
            M = sd[:, None] * mu_new[j0:j1, bcode].T
            h = (M ** 2 / tc[:, None] + M) ** 2
            h[h <= 0] = 1
            sigma2[j0:j1] = (((X - M) ** 2 - M) / h * sd2[:, None]).sum(0) / ((sd2 ** 2)[:, None] / h).sum(0)
        mu = mu_new
        with np.errstate(divide="ignore", invalid="ignore"):
            theta = mu ** 2 / sigma2[:, None]
            ok = (theta > 0) & (theta_prev > 0)
            r = np.abs(np.log(theta[ok] / theta_prev[ok]))
        j += 1
        delta = np.nanmax(r) if np.any(~np.isnan(r)) else 0
        if verbose:
            print(f"  IRLS iter {j}: delta {delta:.4g}", flush=True)
    return dict(mu=mu, sigma_sq=sigma2, theta=med(theta), iters=j, delta=delta)


def _side(Xblock, sd, mu, sigma2, theta_b, bcode):
    """centered counts, NB weights, variances for a block of features (cells x f)"""
    M = sd[:, None] * mu[:, bcode].T
    tc = theta_b[bcode][:, None]
    W = M + M ** 2 / tc
    W[np.isnan(W) | (W <= 0)] = 1
    V = M + (sd ** 2)[:, None] * sigma2[None, :]
    return Xblock - M, W, V


def wls_block(xg, sd_g, pg, xp, sd_p, pp, bcode):
    """xg: cells x G gene counts, xp: cells x P peak counts (dense); pg/pp: IRLS dicts subset to these features.
    Returns covar, test_stat, pval as G x P arrays (all gene-peak combinations of the block)."""
    Xg, Wg, Vg = _side(xg, sd_g, pg["mu"], pg["sigma_sq"], pg["theta"], bcode)
    Xp, Wp, Vp = _side(xp, sd_p, pp["mu"], pp["sigma_sq"], pp["theta"], bcode)
    s = sd_g * sd_p
    s2 = s ** 2
    nume = (Xg * (s / 1)[:, None] / Wg).T @ (Xp / Wp)
    est = nume / ((s2[:, None] / Wg).T @ (1 / Wp))
    deno = ((Vg / Wg ** 2) * s2[:, None]).T @ (Vp / Wp ** 2)
    ts = nume / np.sqrt(deno)
    return est, ts, 2 * norm.sf(np.abs(ts))
