"""
sdft_v02_tension.py
SDFT revised v0.2 — 𝓣 レジーム間テンション（多様体版）

𝓣 = (T₁, T₂, T₃, T₄, T₅) ∈ ℝ⁵

設計原則（v0.0.6 から継続）：
  1. ハイパーパラメータ（重み係数）を持たない
  2. スカラーに集約しない
  3. 各成分は独立した不整合の次元
  4. 加算・重み付け合成を行わない

【多様体版の各成分】
  T₁ : 体積スケール差       |log Vol(M_A) - log Vol(M_B)|
  T₂ : 測地距離             d_g(p_A, p_B) on k-NN graph
  T₃ : 平行移動の不整合      局所 PCA フレーム間の回転角（簡易版）
  T₄ : ホロノミー            曲率差 |κ_A - κ_B| / max(|κ_A|, |κ_B|)（簡易版）
  T₅ : 埋め込み歪み          1 - σ_min/σ_max of local linear map A→B
"""

from __future__ import annotations

import math
from typing import List, Tuple, NamedTuple

import numpy as np
from scipy.spatial import KDTree
from scipy.linalg import svd

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
# データ構造
# ============================================================

class RegimeTensionVector(NamedTuple):
    T1: float; T2: float; T3: float; T4: float; T5: float
    T1_evidence: str = Evidence.DERIVED
    T2_evidence: str = Evidence.DERIVED
    T3_evidence: str = Evidence.ASSUMED   # 簡易版
    T4_evidence: str = Evidence.ASSUMED   # 簡易版
    T5_evidence: str = Evidence.DERIVED

    def as_array(self) -> np.ndarray:
        return np.array([self.T1, self.T2, self.T3, self.T4, self.T5])

    def summary(self) -> str:
        return (
            f"𝓣 = (\n"
            f"  T₁={self.T1:.6f} [{self.T1_evidence}]  体積スケール差\n"
            f"  T₂={self.T2:.6f} [{self.T2_evidence}]  測地距離\n"
            f"  T₃={self.T3:.6f} [{self.T3_evidence}]  平行移動不整合\n"
            f"  T₄={self.T4:.6f} [{self.T4_evidence}]  ホロノミー（曲率差）\n"
            f"  T₅={self.T5:.6f} [{self.T5_evidence}]  埋め込み歪み\n"
            f")"
        )


# ============================================================
# T₁：体積スケール差
# ============================================================

def calc_T1(
    embedded_A: np.ndarray,
    embedded_B: np.ndarray,
    k: int = 5,
) -> Tuple[float, str]:
    """
    T₁ = |log Vol(M_A) - log Vol(M_B)|

    Vol(M) を k-NN ball 体積の総和で近似する。

    T₁ 小：2レジームのスケールが近い
    T₁ 大：スケール差が大きい（一方が広い空間を占め、他方がコンパクト）

    証拠分類: Derived
    """
    def _log_vol(pts: np.ndarray) -> float:
        N, d = pts.shape
        if N < k + 1:
            return 0.0
        tree = KDTree(pts)
        dists, _ = tree.query(pts, k=k + 1)
        r_k = dists[:, -1]
        r_k = np.maximum(r_k, EPS)
        return float(np.mean(d * np.log(r_k)))

    vol_A = _log_vol(embedded_A)
    vol_B = _log_vol(embedded_B)
    T1 = abs(vol_A - vol_B)

    n_min = min(len(embedded_A), len(embedded_B))
    ev = Evidence.DERIVED if n_min >= 30 else Evidence.ASSUMED
    return float(T1), ev


# ============================================================
# T₂：測地距離
# ============================================================

def calc_T2(
    embedded_A: np.ndarray,
    embedded_B: np.ndarray,
    k: int = 10,
) -> Tuple[float, str]:
    """
    T₂ = d_g(p_A, p_B)

    k-NN グラフ上の Dijkstra 最短経路で測地距離を計算する。

    化学ポテンシャルとの関係：μ = T₂（測地距離がそのまま移動コスト）

    証拠分類: Derived
    """
    from sdft_v02_geometry import calc_geodesic_distance
    T2 = calc_geodesic_distance(embedded_A, embedded_B, k=k)

    n_min = min(len(embedded_A), len(embedded_B))
    ev = Evidence.DERIVED if n_min >= 30 else Evidence.ASSUMED
    return float(T2), ev


# ============================================================
# T₃：平行移動の不整合（簡易版）
# ============================================================

def calc_T3(
    embedded_A: np.ndarray,
    embedded_B: np.ndarray,
    k: int = 10,
) -> Tuple[float, str]:
    """
    T₃ = 局所 PCA フレーム間の回転角（簡易版）

    各レジームの主成分方向を局所的な「フレーム」とみなし、
    そのフレーム間の回転角度を平行移動の不整合として測る。

    ‖R - I‖_F / √d  ∈ [0, 1] に正規化

    R = V_A^T · V_B（局所 PCA の主成分行列間の回転行列）

    T₃ 小：フレームが整合（双方向に自然）
    T₃ 大：フレームが捻れている（非対称）

    ⚠️ v0.2 簡易版。v0.2.1 で Singer-Wu vector diffusion に置換予定。
    証拠分類: Assumed（簡易近似）
    """
    d = min(embedded_A.shape[1], embedded_B.shape[1])
    if d < 2 or len(embedded_A) < d + 1 or len(embedded_B) < d + 1:
        return 0.0, Evidence.ASSUMED

    # PCA フレーム
    cov_A = np.cov(embedded_A.T)
    cov_B = np.cov(embedded_B.T)
    _, V_A = np.linalg.eigh(cov_A)
    _, V_B = np.linalg.eigh(cov_B)

    # 回転行列 R = V_A^T · V_B
    R = V_A.T @ V_B

    # フロベニウスノルムでの差
    diff = np.linalg.norm(R - np.eye(d), 'fro') / math.sqrt(d)
    T3 = float(min(diff, 2.0))  # 理論的最大 ≈ √2 だが安全マージン

    return T3, Evidence.ASSUMED


# ============================================================
# T₄：ホロノミー（簡易版・曲率差）
# ============================================================

def calc_T4(
    embedded_A: np.ndarray,
    embedded_B: np.ndarray,
    k: int = 10,
    n_edges: int = 200,
) -> Tuple[float, str]:
    """
    T₄ = |κ_A - κ_B| / max(|κ_A|, |κ_B|)

    各レジームの Ollivier-Ricci 曲率の平均を計算し、
    その相対差をホロノミーの代理指標とする。

    T₄ 小：曲率が近い（多様体が似た形状）
    T₄ 大：曲率が違う（閉ループで情報が歪む）

    ⚠️ v0.2 簡易版。v0.2.1 で離散ホロノミーに置換予定。
    証拠分類: Assumed（近似推定）
    """
    from sdft_v02_geometry import calc_ollivier_ricci
    kappa_A = calc_ollivier_ricci(embedded_A, k=k, n_edges=n_edges)
    kappa_B = calc_ollivier_ricci(embedded_B, k=k, n_edges=n_edges)

    denom = max(abs(kappa_A), abs(kappa_B), EPS)
    T4 = float(abs(kappa_A - kappa_B) / denom)
    T4 = min(1.0, max(0.0, T4))

    return T4, Evidence.ASSUMED


# ============================================================
# T₅：埋め込み歪み
# ============================================================

def calc_T5(
    embedded_A: np.ndarray,
    embedded_B: np.ndarray,
    k: int = 5,
) -> Tuple[float, str]:
    """
    T₅ = 1 - σ_min(J) / σ_max(J)

    レジーム A → B への局所線形写像 f の Jacobian J の
    特異値比で「情報の潰れ具合」を測る。

    f は A の各点を B の k 最近傍で線形回帰して推定する。

    T₅ = 0：等長写像（情報損失なし）
    T₅ = 1：ランク落ち（情報が完全に潰れる）

    証拠分類: Derived（同次元・十分なサンプル時）
    """
    N_A, d_A = embedded_A.shape
    N_B, d_B = embedded_B.shape

    if d_A != d_B or N_A < d_A + 1 or N_B < d_A + 1:
        return 1.0, Evidence.ASSUMED

    d = d_A

    # A の各点について B での最近傍を見つけ、局所線形写像を推定
    combined = np.vstack([embedded_A, embedded_B])
    tree_B = KDTree(embedded_B)

    # A のサンプル点を使って局所 Jacobian を推定
    n_sample = min(N_A, 100)
    sample_idx = np.random.choice(N_A, n_sample, replace=False)

    singular_ratios = []
    for i in sample_idx:
        # A での局所近傍
        tree_A = KDTree(embedded_A)
        k_local = min(k * 2, N_A - 1)
        _, nbr_A_idx = tree_A.query(embedded_A[i], k=k_local + 1)
        nbr_A_idx = nbr_A_idx[1:]  # 自分を除く

        # 対応する B での最近傍
        _, nbr_B_idx = tree_B.query(embedded_A[nbr_A_idx], k=1)
        nbr_B_idx = np.atleast_1d(nbr_B_idx).flatten()

        if len(nbr_A_idx) < d + 1:
            continue

        # 局所線形回帰 B = J · A + b
        dA = embedded_A[nbr_A_idx] - embedded_A[i]
        dB = embedded_B[nbr_B_idx] - embedded_B[nbr_B_idx[0]]

        try:
            J, _, _, _ = np.linalg.lstsq(dA, dB, rcond=None)
            _, s, _ = svd(J)
            if len(s) > 0 and s[0] > EPS:
                ratio = float(s[-1] / s[0])
                singular_ratios.append(ratio)
        except Exception:
            continue

    if not singular_ratios:
        return 1.0, Evidence.ASSUMED

    # 平均特異値比
    mean_ratio = float(np.mean(singular_ratios))
    T5 = float(max(0.0, min(1.0, 1.0 - mean_ratio)))

    ev = Evidence.DERIVED if n_sample >= 30 else Evidence.ASSUMED
    return T5, ev


# ============================================================
# 統合計算
# ============================================================

def calc_regime_tension_vector(
    embedded_A: np.ndarray,
    embedded_B: np.ndarray,
    k: int = 10,
) -> RegimeTensionVector:
    """
    𝓣 = (T₁, T₂, T₃, T₄, T₅) を計算する。

    設計原則：
      ・ハイパーパラメータなし
      ・スカラー集約しない
      ・各成分は独立した次元
    """
    T1, ev1 = calc_T1(embedded_A, embedded_B, k=min(k, 5))
    T2, ev2 = calc_T2(embedded_A, embedded_B, k=k)
    T3, ev3 = calc_T3(embedded_A, embedded_B, k=k)
    T4, ev4 = calc_T4(embedded_A, embedded_B, k=k)
    T5, ev5 = calc_T5(embedded_A, embedded_B, k=min(k, 5))
    return RegimeTensionVector(
        T1=T1, T2=T2, T3=T3, T4=T4, T5=T5,
        T1_evidence=ev1, T2_evidence=ev2,
        T3_evidence=ev3, T4_evidence=ev4,
        T5_evidence=ev5,
    )


# ============================================================
# 動作確認
# ============================================================

if __name__ == "__main__":
    np.random.seed(42)

    print("=" * 60)
    print("SDFT v0.2 Tension — 動作確認")
    print("=" * 60)

    # 同じ分布（𝓣 ≈ 0 のはず）
    A = np.random.randn(200, 3) * 0.5
    B = np.random.randn(200, 3) * 0.5

    print("\n--- 同一分布 ---")
    t_same = calc_regime_tension_vector(A, B, k=10)
    print(t_same.summary())

    # 異なるスケール・位置
    C = np.random.randn(200, 3) * 0.5 + [5, 0, 0]
    D = np.random.randn(200, 3) * 2.0 + [-3, 0, 0]

    print("\n--- 異なるスケール・位置 ---")
    t_diff = calc_regime_tension_vector(C, D, k=10)
    print(t_diff.summary())

    print("\n--- 比較 ---")
    print(f"  T₁: 同一={t_same.T1:.4f} vs 異={t_diff.T1:.4f}  (異なる方が大)")
    print(f"  T₂: 同一={t_same.T2:.4f} vs 異={t_diff.T2:.4f}  (異なる方が大)")

    print("\n" + "=" * 60)
    print("動作確認完了")
