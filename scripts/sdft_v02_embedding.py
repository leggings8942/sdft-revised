"""
sdft_v02_embedding.py
SDFT revised v0.2 — Takens 遅延座標埋め込みと Layer 1 構造場変数

時系列 x(t) を Takens の埋め込み定理に基づいて d 次元空間に持ち上げ、
軌道が描く多様体 M の構造的特徴量を抽出する。

【Layer 1：構造場変数】
  S_g  = mean log Vol(k-NN ball)   軌道体積エントロピー   Derived
  D    = Grassberger-Procaccia     相関次元               Derived
  λ₁   = Rosenstein                最大リアプノフ指数     Derived

【Takens 埋め込み】
  τ  = 相互情報量の最初の極小値   Derived (auto) / Assumed (manual)
  d  = False Nearest Neighbors     Derived (auto) / Assumed (manual)

設計原則：
  1. 確率分布を使わない（密度推定は v0.3 以降）
  2. 計算は埋め込み後の点群で完結する
  3. 全変数に証拠分類ラベルを付与する
"""

from __future__ import annotations

import math
from typing import List, Tuple, Optional

import numpy as np
from scipy.spatial import KDTree

EPS = 1e-12


# ============================================================
# 証拠分類
# ============================================================

class Evidence:
    OBSERVED   = "Observed"
    DERIVED    = "Derived"
    ASSUMED    = "Assumed"
    HYPOTHESIS = "Hypothesis"


# ============================================================
# Takens 遅延座標埋め込み
# ============================================================

def estimate_delay(
    x: List[float],
    max_lag: int = 50,
    bins: int = 16,
) -> int:
    """
    相互情報量の最初の極小値から最適遅延 τ を推定する。

    References: Fraser & Swinney (1986)
    証拠分類: Derived
    """
    arr = np.asarray(x, dtype=float)
    n = len(arr)
    max_lag = min(max_lag, n // 4)

    if max_lag < 2:
        return 1

    mi_values = []
    for lag in range(1, max_lag):
        x1 = arr[:n - lag]
        x2 = arr[lag:]
        H2d, _, _ = np.histogram2d(x1, x2, bins=bins)
        pxy = H2d / (H2d.sum() + EPS)
        px = pxy.sum(axis=1)
        py = pxy.sum(axis=0)
        hx  = -np.sum(px[px > 0] * np.log(px[px > 0] + EPS))
        hy  = -np.sum(py[py > 0] * np.log(py[py > 0] + EPS))
        hxy = -np.sum(pxy[pxy > 0] * np.log(pxy[pxy > 0] + EPS))
        mi_values.append(hx + hy - hxy)

    # 最初の局所極小を探す
    for i in range(1, len(mi_values) - 1):
        if mi_values[i] < mi_values[i - 1] and mi_values[i] <= mi_values[i + 1]:
            return i + 1

    return 1


def estimate_embedding_dim(
    x: List[float],
    tau: int,
    max_dim: int = 10,
    r_threshold: float = 10.0,
    fnn_threshold: float = 0.02,
) -> int:
    """
    False Nearest Neighbors (FNN) 法で最適埋め込み次元 d を推定する。

    References: Kennel, Brown & Abarbanel (1992)
    証拠分類: Derived
    """
    arr = np.asarray(x, dtype=float)
    n = len(arr)

    for d in range(1, max_dim + 1):
        N_d  = n - d * tau
        N_d1 = n - (d + 1) * tau  # d+1 次元用
        if N_d1 < 10:
            return d

        emb_d  = np.column_stack([arr[i * tau:i * tau + N_d1] for i in range(d)])
        emb_d1 = np.column_stack([arr[i * tau:i * tau + N_d1] for i in range(d + 1)])

        tree = KDTree(emb_d)
        dists, indices = tree.query(emb_d, k=2)
        nn_dist = dists[:, 1]
        nn_idx  = indices[:, 1]

        fnn_count = 0
        total     = 0
        for i in range(N_d1):
            j = nn_idx[i]
            if nn_dist[i] < EPS:
                continue
            extra_dist = abs(emb_d1[i, -1] - emb_d1[j, -1])
            if extra_dist / (nn_dist[i] + EPS) > r_threshold:
                fnn_count += 1
            total += 1

        fnn_frac = fnn_count / (total + 1)
        if fnn_frac < fnn_threshold:
            return d

    return max_dim


def takens_embed(
    x: List[float],
    d: int,
    tau: int,
) -> np.ndarray:
    """
    時系列を d 次元遅延座標空間に埋め込む。

    y(t) = (x(t), x(t+τ), x(t+2τ), ..., x(t+(d-1)τ))

    Returns: (N, d) ndarray,  N = len(x) - (d-1)*τ
    証拠分類: Derived
    """
    arr = np.asarray(x, dtype=float)
    n = len(arr)
    N = n - (d - 1) * tau
    if N < 1:
        raise ValueError(f"系列長 {n} が短すぎます (d={d}, τ={tau})")
    return np.column_stack([arr[i * tau:i * tau + N] for i in range(d)])


def auto_embed(
    x: List[float],
    max_dim: int = 10,
    max_lag: int = 50,
) -> Tuple[np.ndarray, int, int]:
    """
    遅延 τ と次元 d を自動推定して埋め込みを返す。

    Returns: (embedded, d, tau)
    """
    tau = estimate_delay(x, max_lag=max_lag)
    d   = estimate_embedding_dim(x, tau, max_dim=max_dim)
    emb = takens_embed(x, d, tau)
    return emb, d, tau


# ============================================================
# Layer 1: S_g（軌道体積エントロピー）
# ============================================================

def calc_volume_entropy(
    embedded: np.ndarray,
    k: int = 5,
) -> float:
    """
    軌道体積エントロピーを計算する。

    S_g = mean_i [ log Vol(B_k(y_i)) ]
        ≈ mean_i [ d · log r_k(y_i) + const ]

    r_k(y_i) は点 y_i の k 番目最近傍距離。
    d-球の体積 V_d(r) ∝ r^d なので log V ∝ d·log(r)。

    S_g 大：軌道が空間に広がっている（散逸）
    S_g 小：軌道がコンパクト（秩序）

    証拠分類: Derived
    """
    N, d = embedded.shape
    if N < k + 1:
        return 0.0
    tree = KDTree(embedded)
    dists, _ = tree.query(embedded, k=k + 1)
    r_k = dists[:, -1]
    r_k = np.maximum(r_k, EPS)
    log_vols = d * np.log(r_k)
    return float(np.mean(log_vols))


# ============================================================
# Layer 1: D（相関次元・Grassberger-Procaccia）
# ============================================================

def calc_correlation_dim(
    embedded: np.ndarray,
    n_r: int = 15,
    max_samples: int = 2000,
) -> float:
    """
    Grassberger-Procaccia 法により相関次元 D を計算する。

    C(r) = (2 / N(N-1)) Σ_{i<j} Θ(r - ‖y_i - y_j‖)

    D = lim_{r→0} d log C(r) / d log r

    References: Grassberger & Procaccia (1983)
    証拠分類: Derived
    """
    N, d = embedded.shape
    if N < 30:
        return 0.0

    # 大きなデータはサブサンプル
    if N > max_samples:
        idx = np.random.choice(N, max_samples, replace=False)
        pts = embedded[idx]
    else:
        pts = embedded

    M = len(pts)
    tree = KDTree(pts)

    # ペア距離からスケール範囲を決定
    sample_idx = np.random.choice(M, min(M, 200), replace=False)
    sample_dists = []
    for i in sample_idx:
        dd, _ = tree.query(pts[i], k=min(M, 10))
        sample_dists.extend(dd[1:].tolist())
    sample_dists = np.array(sample_dists)
    sample_dists = sample_dists[sample_dists > EPS]

    if len(sample_dists) < 10:
        return 0.0

    r_min = float(np.percentile(sample_dists, 5))
    r_max = float(np.percentile(sample_dists, 60))
    if r_min <= 0 or r_max <= r_min * 2:
        return 0.0

    r_values = np.logspace(np.log10(r_min), np.log10(r_max), n_r)

    log_r = []
    log_C = []
    total_pairs = M * (M - 1) / 2.0

    for r in r_values:
        count = 0
        for i in range(M):
            neighbors = tree.query_ball_point(pts[i], r)
            count += len(neighbors) - 1  # 自分を除く
        C_r = count / (2.0 * total_pairs + EPS)
        if C_r > EPS:
            log_r.append(np.log(r))
            log_C.append(np.log(C_r))

    if len(log_r) < 4:
        return 0.0

    coeffs = np.polyfit(log_r, log_C, 1)
    D = float(coeffs[0])
    return max(0.0, D)


# ============================================================
# Layer 1: λ₁（最大リアプノフ指数・Rosenstein 法）
# ============================================================

def calc_lyapunov(
    embedded: np.ndarray,
    dt: float = 1.0,
    min_temporal_sep: Optional[int] = None,
    max_iter: Optional[int] = None,
) -> float:
    """
    Rosenstein 法により最大リアプノフ指数 λ₁ を推定する。

    各点の最近傍（時間的に離れた）を見つけ、時間発展に沿った
    距離の発散率を線形近似する。

    λ₁ < 0：安定（軌道が収束）→ 予測可能・持続的
    λ₁ ≈ 0：中立（限周期軌道）
    λ₁ > 0：カオス（軌道が発散）→ 予測不能

    References: Rosenstein, Collins & De Luca (1993)
    証拠分類: Derived
    """
    N, d = embedded.shape
    if N < 20:
        return 0.0

    if min_temporal_sep is None:
        min_temporal_sep = d + 1
    if max_iter is None:
        max_iter = min(N // 4, 200)

    tree = KDTree(embedded)
    divergence = np.zeros(max_iter)
    counts     = np.zeros(max_iter)

    # 各点から最近傍（時間的に離れた）を見つける
    k_search = min(N, min_temporal_sep * 3 + 5)
    dists_all, indices_all = tree.query(embedded, k=k_search)

    for i in range(N - max_iter):
        nn_idx = -1
        for ki in range(1, k_search):
            j = indices_all[i, ki]
            if abs(i - j) >= min_temporal_sep and j + max_iter <= N:
                nn_idx = j
                break
        if nn_idx < 0:
            continue

        d0 = np.linalg.norm(embedded[i] - embedded[nn_idx])
        if d0 < EPS:
            continue

        for step in range(max_iter):
            ii = i + step
            jj = nn_idx + step
            if ii >= N or jj >= N:
                break
            dk = np.linalg.norm(embedded[ii] - embedded[jj])
            if dk > EPS:
                divergence[step] += np.log(dk)
                counts[step] += 1

    valid = counts > 0
    if not np.any(valid):
        return 0.0

    avg_div = np.full_like(divergence, np.nan)
    avg_div[valid] = divergence[valid] / counts[valid]

    # 初期線形領域をフィット
    t_vals = np.where(valid)[0].astype(float) * dt
    y_vals = avg_div[valid]

    if len(t_vals) < 3:
        return 0.0

    n_fit = max(3, len(t_vals) // 3)
    coeffs = np.polyfit(t_vals[:n_fit], y_vals[:n_fit], 1)
    return float(coeffs[0])


# ============================================================
# 統合解析関数
# ============================================================

def analyze_layer1(
    x: List[float],
    d: Optional[int] = None,
    tau: Optional[int] = None,
    k_entropy: int = 5,
) -> dict:
    """
    時系列 x に対して Layer 1 の全変数を計算する。

    Parameters
    ----------
    x   : 時系列データ
    d   : 埋め込み次元（None で自動推定）
    tau : 遅延（None で自動推定）

    Returns
    -------
    dict: embedded, d, tau, S_g, D, lambda_1 と証拠分類
    """
    arr = np.asarray(x, dtype=float)

    # 埋め込みパラメータ
    if tau is None:
        tau = estimate_delay(arr)
        tau_ev = Evidence.DERIVED
    else:
        tau_ev = Evidence.ASSUMED

    if d is None:
        d = estimate_embedding_dim(arr, tau)
        d_ev = Evidence.DERIVED
    else:
        d_ev = Evidence.ASSUMED

    embedded = takens_embed(arr, d, tau)

    S_g      = calc_volume_entropy(embedded, k=k_entropy)
    D        = calc_correlation_dim(embedded)
    lambda_1 = calc_lyapunov(embedded)

    return {
        'embedded':   embedded,
        'd':          d,
        'tau':        tau,
        'd_evidence': d_ev,
        'tau_evidence': tau_ev,
        'S_g':        S_g,
        'D':          D,
        'lambda_1':   lambda_1,
        'S_g_evidence':      Evidence.DERIVED,
        'D_evidence':        Evidence.DERIVED,
        'lambda_1_evidence': Evidence.DERIVED,
    }


# ============================================================
# 動作確認
# ============================================================

if __name__ == "__main__":
    np.random.seed(42)
    n = 1000

    # ① Lorenz アトラクタ風（カオス系、λ₁ > 0）
    from scipy.integrate import solve_ivp

    def lorenz(t, state, sigma=10, rho=28, beta=8/3):
        x, y, z = state
        return [sigma*(y-x), x*(rho-z)-y, x*y-beta*z]

    sol = solve_ivp(lorenz, [0, 50], [1, 1, 1], max_step=0.05,
                    t_eval=np.linspace(0, 50, n))
    x_lorenz = sol.y[0]

    print("=" * 60)
    print("SDFT v0.2 Embedding — 動作確認")
    print("=" * 60)

    # Lorenz x 成分の解析
    print("\n--- Lorenz x成分 ---")
    tau_l = estimate_delay(x_lorenz)
    d_l   = estimate_embedding_dim(x_lorenz, tau_l)
    print(f"  τ = {tau_l},  d = {d_l}")

    emb_l = takens_embed(x_lorenz, d_l, tau_l)
    print(f"  embedded shape: {emb_l.shape}")

    S_g_l = calc_volume_entropy(emb_l)
    D_l   = calc_correlation_dim(emb_l)
    lam_l = calc_lyapunov(emb_l)
    print(f"  S_g = {S_g_l:.4f}")
    print(f"  D   = {D_l:.4f}  (Lorenz 期待値 ≈ 2.05)")
    print(f"  λ₁  = {lam_l:.4f}  (期待: > 0, カオス)")

    # ② 正弦波（規則的、λ₁ ≈ 0 or < 0）
    t = np.linspace(0, 20 * np.pi, n)
    x_sine = np.sin(t) + 0.05 * np.random.randn(n)

    print("\n--- 正弦波 + ノイズ ---")
    result_sine = analyze_layer1(x_sine)
    print(f"  τ = {result_sine['tau']},  d = {result_sine['d']}")
    print(f"  S_g = {result_sine['S_g']:.4f}")
    print(f"  D   = {result_sine['D']:.4f}  (期待: ≈ 1.0)")
    print(f"  λ₁  = {result_sine['lambda_1']:.4f}  (期待: ≤ 0)")

    # ③ ホワイトノイズ（高次元充填、λ₁ → ∞）
    x_noise = np.random.randn(n)
    print("\n--- ホワイトノイズ ---")
    result_noise = analyze_layer1(x_noise)
    print(f"  τ = {result_noise['tau']},  d = {result_noise['d']}")
    print(f"  S_g = {result_noise['S_g']:.4f}")
    print(f"  D   = {result_noise['D']:.4f}  (期待: → d)")
    print(f"  λ₁  = {result_noise['lambda_1']:.4f}  (期待: > 0)")

    print("\n" + "=" * 60)
    print("動作確認完了")
