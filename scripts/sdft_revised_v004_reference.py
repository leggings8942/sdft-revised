"""
sdft_revised_v004_reference.py
SDFT revised v0.0.4 参照実装

設計思想：
  - 後方互換性なし（v0.0.3 以前の実装と互換性を持たない）
  - 全変数を情報幾何学・熱力学的アナロジーに基づいて定義
  - グランドポテンシャル Φ による介入最適化の設計方針を実装
  - 観測可能性を最優先（全変数に証拠分類ラベルを付与）

コア変数：
  S  = k_B·ln2·H_bit          [J/K]   Derived
  D  = higuchi_fd(x)           [1-2]   Derived
  H  = hurst_exponent(x)       [0-1]   Derived
  G  = Fisher情報行列           [-]     Derived
  T  = tr(G⁻¹)/(d·k_B)        [K]     Derived
  U  = tr(G)/d                 [J相当] Derived
  F  = U - T·S                 [J]     Derived
  Φ  = F - μ·N                 [J]     Derived/Observed
  𝓣  = (T₁,T₂,T₃,T₄,T₅)      [R⁵]    各成分による
  μ  = 2√T₂                   [J/人]  Derived
  N  = 移動する人数             [人]    Observed
  P  = sign(dΦ)·‖dΦ‖_G        [J]     Derived
  V  = ∫√detG dθ               [m³相当] Derived
"""

from __future__ import annotations

import math
import statistics
from dataclasses import dataclass, field
from typing import List, Tuple, Optional, NamedTuple
import numpy as np

# ============================================================
# 物理定数
# ============================================================

K_B: float = 1.380649e-23   # ボルツマン定数 [J/K]
LN2: float = math.log(2)    # ln(2)
EPS: float = 1e-12


# ============================================================
# 証拠分類ラベル
# ============================================================

class Evidence:
    OBSERVED   = "Observed"    # 実データから直接計測
    DERIVED    = "Derived"     # 観測値から確定的に算出
    ASSUMED    = "Assumed"     # 分析者が設定した仮定値
    HYPOTHESIS = "Hypothesis"  # 解釈的・推論的推定


# ============================================================
# データ構造
# ============================================================

@dataclass
class SDFTState:
    """
    単一時点での SDFT 状態変数をまとめたデータクラス。
    全変数に証拠分類ラベルを付与する。
    """
    # 構造場変数（時系列から計算）
    H_bit: float = 0.0   # シャノンエントロピー [bit]
    S:     float = 0.0   # 熱力学的エントロピー [J/K]
    D:     float = 1.0   # フラクタル次元 [1.0-2.0]
    H:     float = 0.5   # Hurst 指数 [0.0-1.0]

    # Fisher 行列由来（NumPy 配列として保持）
    G:     Optional[np.ndarray] = None  # Fisher情報行列 [d×d]
    G_inv: Optional[np.ndarray] = None  # G の逆行列

    # Fisher 行列から導出されるスカラー量
    T:     float = 0.0   # 情報的温度 [K]
    U:     float = 0.0   # 内部エネルギー [J相当]
    F:     float = 0.0   # 自由エネルギー [J]

    # グランドポテンシャル関連
    N:     float = 0.0   # 移動する人数 [人]
    mu:    float = 0.0   # 化学ポテンシャル = 2√T₂ [J/人]
    Phi:   float = 0.0   # グランドポテンシャル [J]
    V:     float = 0.0   # 統計多様体の体積 [m³相当]

    # 証拠分類
    evidence_S:   str = Evidence.DERIVED
    evidence_D:   str = Evidence.DERIVED
    evidence_H:   str = Evidence.DERIVED
    evidence_T:   str = Evidence.DERIVED
    evidence_U:   str = Evidence.DERIVED
    evidence_F:   str = Evidence.DERIVED
    evidence_N:   str = Evidence.OBSERVED
    evidence_mu:  str = Evidence.DERIVED
    evidence_Phi: str = Evidence.DERIVED
    evidence_V:   str = Evidence.DERIVED


class RegimeTensionVector(NamedTuple):
    """
    レジーム間テンション 𝓣 の5次元ベクトル。
    集約・スカラー化を行わない。各成分を独立した次元として保持する。
    """
    T1: float  # 対数分配関数差（スケール差）
    T2: float  # Bhattacharyya距離（状態点の幾何学的ずれ）
    T3: float  # α-ダイバージェンス（非対称な意味的差・α≠0）
    T4: float  # 統計的ホロノミー（多様体の曲率差）
    T5: float  # 通信路損失率（旧社会層）

    T1_evidence: str = Evidence.DERIVED
    T2_evidence: str = Evidence.DERIVED
    T3_evidence: str = Evidence.DERIVED
    T4_evidence: str = Evidence.ASSUMED
    T5_evidence: str = Evidence.DERIVED

    T3_alpha: float = 1.0  # α ≠ 0 を必ず明示

    def as_array(self) -> np.ndarray:
        return np.array([self.T1, self.T2, self.T3, self.T4, self.T5])

    def summary(self) -> str:
        return (
            f"𝓣 = (\n"
            f"  T₁={self.T1:.6f} [{self.T1_evidence}]  スケール差\n"
            f"  T₂={self.T2:.6f} [{self.T2_evidence}]  幾何学的距離\n"
            f"  T₃={self.T3:.6f} [{self.T3_evidence}] α={self.T3_alpha}  非対称差\n"
            f"  T₄={self.T4:.6f} [{self.T4_evidence}]  曲率差\n"
            f"  T₅={self.T5:.6f} [{self.T5_evidence}]  通信路損失\n"
            f")"
        )


@dataclass
class InterventionResult:
    """介入最適化の結果"""
    phi_current: float       # 現在の Φ
    phi_predicted: float     # τ 後の予測 Φ
    phi_threshold: float     # 閾値 Φ₀
    needs_intervention: bool # 介入が必要か
    delta_phi_need: float    # 必要な Φ 回復量

    # 各操作変数の感度
    sensitivity_N:  float = 0.0   # ∂Φ/∂N
    sensitivity_T2: float = 0.0   # ∂Φ/∂T₂
    sensitivity_theta: float = 0.0 # ∂Φ/∂θ（近似）

    # 各操作変数の介入効率（コストで割った感度）
    efficiency_N:  float = 0.0
    efficiency_T2: float = 0.0
    efficiency_theta: float = 0.0

    # 優先操作変数と必要操作量
    best_variable: str = ""
    delta_u_optimal: float = 0.0

    # priority_score
    priority_N:  float = 0.0
    priority_T2: float = 0.0
    priority_theta: float = 0.0


# ============================================================
# Layer 1: S（シャノンエントロピー → 熱力学的エントロピー）
# ============================================================

def calc_H_bit(x: List[float], bins: int = 32) -> float:
    """
    シャノンエントロピー H_bit を計算する。

    H_bit = -Σ p_i log₂(p_i)  [bit]

    証拠分類：Derived（時系列データが存在する場合）

    Parameters
    ----------
    x    : 時系列データ（数値リスト）
    bins : ヒストグラムのビン数

    Returns
    -------
    H_bit : float  [bit]
    """
    if len(x) < 2:
        return 0.0
    min_v, max_v = min(x), max(x)
    if max_v - min_v < EPS:
        return 0.0
    bin_width = (max_v - min_v) / bins
    counts = [0] * bins
    for v in x:
        idx = min(int((v - min_v) / bin_width), bins - 1)
        counts[idx] += 1
    n = len(x)
    h_bit = 0.0
    for c in counts:
        if c > 0:
            p = c / n
            h_bit -= p * math.log2(p + EPS)
    return float(h_bit)


def calc_S(x: List[float], bins: int = 32) -> float:
    """
    熱力学的エントロピー S を計算する。

    S = k_B · ln(2) · H_bit  [J/K]

    シャノンエントロピー H_bit（bit 単位）を
    ボルツマン定数 k_B と ln(2) を介して
    熱力学的エントロピー S [J/K] に変換する。

    証拠分類：Derived

    Returns
    -------
    S : float  [J/K]
    """
    H_bit = calc_H_bit(x, bins)
    return float(K_B * LN2 * H_bit)


# ============================================================
# Layer 1: D（フラクタル次元）
# ============================================================

def calc_D(x: List[float], k_max: int = 10) -> float:
    """
    ヒグチ法によるフラクタル次元 D を計算する。

    D ∈ [1.0, 2.0]
      D → 1.0：単純な曲線（低複雑性）
      D → 2.0：完全ランダム（高複雑性）

    証拠分類：Derived

    Parameters
    ----------
    x     : 時系列データ
    k_max : 計算するスケールの最大数

    Returns
    -------
    D : float  [1.0, 2.0]
    """
    n = len(x)
    if n < k_max + 2:
        return 1.0
    lk = []
    for k in range(1, k_max + 1):
        lmk = []
        for m in range(1, k + 1):
            idxs = list(range(m - 1, n, k))
            if len(idxs) < 2:
                continue
            length = sum(
                abs(x[idxs[i]] - x[idxs[i - 1]])
                for i in range(1, len(idxs))
            )
            norm = (n - 1) / (len(idxs) - 1) / k
            lmk.append(length * norm)
        if lmk:
            lk.append((math.log(k + EPS), math.log(statistics.mean(lmk) + EPS)))
    if len(lk) < 2:
        return 1.0
    xs = [p[0] for p in lk]
    ys = [p[1] for p in lk]
    x_mean = statistics.mean(xs)
    y_mean = statistics.mean(ys)
    num = sum((xi - x_mean) * (yi - y_mean) for xi, yi in zip(xs, ys))
    den = sum((xi - x_mean) ** 2 for xi in xs)
    fd = abs(num / (den + EPS))
    return float(max(1.0, min(2.0, fd)))


# ============================================================
# Layer 1: H（Hurst 指数）
# ============================================================

def calc_H(x: List[float], max_lag: int = 20) -> float:
    """
    Hurst 指数 H を計算する。

    H ∈ [0.0, 1.0]
      H > 0.5：持続性（トレンドが続く傾向）
      H = 0.5：ランダムウォーク
      H < 0.5：反持続性（反転しやすい）

    証拠分類：Derived

    Parameters
    ----------
    x       : 時系列データ
    max_lag : ラグの最大値

    Returns
    -------
    H : float  [0.0, 1.0]
    """
    n = len(x)
    if n < max_lag + 2:
        return 0.5
    lags = range(2, min(max_lag, n // 2) + 1)
    log_lags, log_stds = [], []
    for lag in lags:
        diffs = [x[i] - x[i - lag] for i in range(lag, n)]
        if len(diffs) < 2:
            continue
        std = statistics.stdev(diffs)
        if std > EPS:
            log_lags.append(math.log(lag))
            log_stds.append(math.log(std))
    if len(log_lags) < 2:
        return 0.5
    x_mean = statistics.mean(log_lags)
    y_mean = statistics.mean(log_stds)
    num = sum((xi - x_mean) * (yi - y_mean) for xi, yi in zip(log_lags, log_stds))
    den = sum((xi - x_mean) ** 2 for xi in log_lags)
    slope = num / (den + EPS)
    return float(max(0.0, min(1.0, slope / 2.0)))


# ============================================================
# Layer 2: Fisher 情報行列
# ============================================================

def calc_fisher_matrix(x: List[float], d: int = 2) -> np.ndarray:
    """
    時系列データから Fisher 情報行列 G を経験的に推定する。

    正規分布 N(μ, σ²) を仮定した場合の Fisher 行列（閉形式）：
      G = [[1/σ², 0     ],
           [0,     2/σ⁴]]

    パラメータ θ = (μ, σ²) に対する Fisher 行列。

    証拠分類：Derived（正規分布仮定のもとで）

    Parameters
    ----------
    x : 時系列データ
    d : パラメータ次元（正規分布の場合は 2）

    Returns
    -------
    G : np.ndarray  [d×d]
    """
    arr = np.asarray(x, dtype=float)
    if len(arr) < 2:
        return np.eye(d) * EPS
    mu_hat  = float(np.mean(arr))
    var_hat = float(np.var(arr, ddof=1))
    if var_hat < EPS:
        var_hat = EPS
    # 正規分布の Fisher 情報行列（2×2）
    G = np.array([
        [1.0 / var_hat,         0.0          ],
        [0.0,           2.0 / (var_hat ** 2) ]
    ])
    if d != 2:
        # d 次元の単位行列でスケーリング（近似）
        G = np.eye(d) * (np.trace(G) / 2.0)
    return G


def calc_T_temperature(G: np.ndarray) -> float:
    """
    Fisher 行列から情報的温度 T を計算する。

    T = tr(G⁻¹) / (d · k_B)  [K]

    理論的背景：
      Jaynes の最大エントロピー原理により、
      ボルツマン分布と指数型分布族は数学的に同型。
      k_B·T ↔ tr(G⁻¹)/d（系の揺らぎのエネルギースケール）
      の対応が厳密に導かれる。

    証拠分類：Derived

    Returns
    -------
    T : float  [K]
    """
    try:
        G_inv = np.linalg.inv(G + np.eye(G.shape[0]) * EPS)
        d = G.shape[0]
        return float(np.trace(G_inv) / (d * K_B))
    except np.linalg.LinAlgError:
        return 0.0


def calc_U_internal_energy(G: np.ndarray) -> float:
    """
    Fisher 行列から内部エネルギー U を計算する。

    U = tr(G) / d  [J相当]

    理論的背景：
      等分配則により各自由度に k_B·T/2 のエネルギーが分配される。
      U ~ k_B·T·tr(G) と Fisher 行列のトレースが対応する。
      G が大きい（感度が高い）= エネルギーランドスケープが急峻 = U が大きい。

    証拠分類：Derived

    Returns
    -------
    U : float  [J相当]
    """
    d = G.shape[0]
    return float(np.trace(G) / d)


def calc_F_free_energy(G: np.ndarray, H_bit: float) -> float:
    """
    Fisher 行列とシャノンエントロピーから自由エネルギー F を計算する。

    F = U - T·S
      = tr(G)/d - ln2·H_bit·tr(G⁻¹)/d
      [J]

    証拠分類：Derived

    Returns
    -------
    F : float  [J]
    """
    d = G.shape[0]
    U = calc_U_internal_energy(G)
    try:
        G_inv = np.linalg.inv(G + np.eye(d) * EPS)
        T_S = LN2 * H_bit * float(np.trace(G_inv)) / d
    except np.linalg.LinAlgError:
        T_S = 0.0
    return float(U - T_S)


def calc_V_volume(G: np.ndarray, theta_range: float = 1.0) -> float:
    """
    統計多様体の体積 V を計算する。

    V = ∫_Ω √det(G(θ)) dθ  [m³相当]

    実用的な近似：
      V ≈ √det(G) · θ_range^d
      （θ の積分範囲を θ_range で近似）

    理論的背景：
      Jeffreys 事前分布密度 √det(G) をパラメータ空間で積分した量。
      Fisher 行列が大きい（感度が高い）= 多様体の体積が大きい。

    証拠分類：Derived（θ_range は Assumed）

    Returns
    -------
    V : float  [m³相当]
    """
    try:
        det_G = max(float(np.linalg.det(G)), EPS)
        d = G.shape[0]
        return float(math.sqrt(det_G) * (theta_range ** d))
    except np.linalg.LinAlgError:
        return EPS


# ============================================================
# Layer 2: 𝓣 ベクトルの各成分
# ============================================================

def _shared_bins_distribution(
    data_A: List[float],
    data_B: List[float],
    bins: int = 50
) -> Tuple[np.ndarray, np.ndarray]:
    """共通ビンで両レジームの経験分布を推定する（内部関数）"""
    a = np.asarray(data_A, dtype=float)
    b = np.asarray(data_B, dtype=float)
    all_data = np.concatenate([a, b])
    bin_edges = np.linspace(all_data.min() - EPS, all_data.max() + EPS, bins + 1)
    p_A, _ = np.histogram(a, bins=bin_edges)
    p_B, _ = np.histogram(b, bins=bin_edges)
    p_A = p_A.astype(float) / (p_A.sum() + EPS)
    p_B = p_B.astype(float) / (p_B.sum() + EPS)
    return p_A, p_B


def calc_T1(data_A: List[float], data_B: List[float], bins: int = 50) -> Tuple[float, str]:
    """
    T₁：対数分配関数差（スケール差）

    T₁ = |ψ_A - ψ_B|

    ψ(θ) = log Z(θ) は対数分配関数（ポテンシャル関数）。
    ψ ≈ -E[log p(x)] で近似（負の対数尤度の経験平均）。

    意味：2つのレジームの状態空間のスケール・自由度の乖離。
    原典対応：|log(L_A/L_B)|（ゲージ理論の特性長比）
    証拠分類：Derived
    """
    a = np.asarray(data_A, dtype=float)
    b = np.asarray(data_B, dtype=float)

    def empirical_psi(arr: np.ndarray, bins: int = bins) -> float:
        p, _ = np.histogram(arr, bins=bins, density=False)
        p = p.astype(float) / (p.sum() + EPS)
        p = np.clip(p, EPS, None)
        log_p = np.log(p)
        indices = np.clip(
            ((arr - arr.min()) / ((arr.max() - arr.min() + EPS) / bins)).astype(int),
            0, bins - 1
        )
        return float(-np.mean(log_p[indices]))

    psi_A = empirical_psi(a)
    psi_B = empirical_psi(b)
    T1 = abs(psi_A - psi_B)
    ev = Evidence.DERIVED if (len(data_A) >= 10 and len(data_B) >= 10) else Evidence.ASSUMED
    return float(T1), ev


def calc_T2(data_A: List[float], data_B: List[float], bins: int = 50) -> Tuple[float, str]:
    """
    T₂：Bhattacharyya 距離（状態点の幾何学的ずれ）

    T₂ = d_B(p_A, p_B) = -log ∫√(p_A(x)·p_B(x)) dx = -log BC(p_A, p_B)

    Fisher-Rao 測地距離との関係：
      d_FR(p_A, p_B) ≈ 2√T₂  （単調近似）

    T₂ と T₃ の独立性：
      T₂ は対称（BC ベース）、T₃ は非対称（α≠0）
      α≠0 を守る限り T₂ と T₃ は独立した情報を持つ

    意味：統計多様体上での2点間の幾何学的距離。
          化学ポテンシャル μ = 2√T₂ に直結。
    原典対応：κ·ΔTopo（トポロジー差）
    証拠分類：Derived
    """
    p_A, p_B = _shared_bins_distribution(data_A, data_B, bins)
    BC = float(np.sum(np.sqrt(p_A * p_B)))
    BC = max(BC, EPS)
    T2 = float(-math.log(BC))
    ev = Evidence.DERIVED if (len(data_A) >= 10 and len(data_B) >= 10) else Evidence.ASSUMED
    return T2, ev


def calc_T3(
    data_A: List[float],
    data_B: List[float],
    alpha: float = 1.0,
    bins: int = 50
) -> Tuple[float, str]:
    """
    T₃：α-ダイバージェンス（非対称な意味的差）

    T₃ = D^(α)(p_A ‖ p_B) = 4/(1-α²) · (1 - ∫p_A^((1-α)/2)·p_B^((1+α)/2) dx)

    ⚠️ α = 0 は禁止（T₂ の Bhattacharyya と重複するため）

    α の意味：
      α = +1：KL ダイバージェンス D_KL(p_A‖p_B)  （A の視点から B の意外さ）
      α = -1：逆 KL ダイバージェンス               （B の視点から A の意外さ）

    意味：意味変換の非対称な差。「どちらの視点から見た差か」を α で制御。
    原典対応：λ·ΔOp（演算子差）
    証拠分類：Derived
    """
    assert abs(alpha) > 1e-6, (
        "α=0 は禁止（T₂ との重複が発生）。α=1 または α=-1 を使用してください。"
    )
    p_A, p_B = _shared_bins_distribution(data_A, data_B, bins)
    p_A = np.clip(p_A, EPS, None)
    p_B = np.clip(p_B, EPS, None)

    if abs(alpha - 1.0) < 1e-6:
        T3 = float(np.sum(p_A * np.log(p_A / p_B)))
    elif abs(alpha + 1.0) < 1e-6:
        T3 = float(np.sum(p_B * np.log(p_B / p_A)))
    else:
        integral = float(np.sum(p_A ** ((1 - alpha) / 2) * p_B ** ((1 + alpha) / 2)))
        T3 = float(4 / (1 - alpha ** 2) * (1 - integral))

    T3 = max(0.0, T3)
    ev = Evidence.DERIVED if (len(data_A) >= 10 and len(data_B) >= 10) else Evidence.ASSUMED
    return T3, ev


def calc_T4(
    data_A: List[float],
    data_B: List[float],
    alpha: float = 1.0
) -> Tuple[float, str]:
    """
    T₄：統計的ホロノミー（多様体の曲率差）

    T₄ = ‖Hol^(α)(γ) - I‖

    指数型分布族での曲率テンソル近似：
      R_{ijkl} = (1-α²)/4 · (g_il·g_jk - g_ik·g_jl)

    正規分布近似（1次元）での近似推定：
      G_A = 1/var(A)（スカラー Fisher 情報）
      R_A = (1-α²)/4 · G_A²
      T₄ ≈ |R_A - R_B| / max(R_A, R_B)

    意味：複数のレジームを巡る閉ループで整合性が保たれるか。
    原典対応：η·‖hol(W_cycle)-I‖（四元数ホロノミー）
    証拠分類：Assumed（近似推定）
    """
    a = np.asarray(data_A, dtype=float)
    b = np.asarray(data_B, dtype=float)
    var_A = max(float(np.var(a)), EPS)
    var_B = max(float(np.var(b)), EPS)
    G_A = 1.0 / var_A
    G_B = 1.0 / var_B
    R_A = (1 - alpha ** 2) / 4 * G_A ** 2
    R_B = (1 - alpha ** 2) / 4 * G_B ** 2
    T4 = float(abs(R_A - R_B) / (max(R_A, R_B) + EPS))
    T4 = max(0.0, min(1.0, T4))
    return T4, Evidence.ASSUMED


def calc_T5(data_A: List[float], data_B: List[float], bins: int = 50) -> Tuple[float, str]:
    """
    T₅：通信路損失率（旧社会層の統合）

    T₅ = L(e) = 1 - I(X_A; X_B) / H(X_A)
    I(X_A; X_B) = H(X_A) + H(X_B) - H(X_A, X_B)

    L = 0：完全に伝わる（損失ゼロ）
    L = 1：全く伝わらない

    旧社会層との対応：
      b_norm（規範的障壁）→ チャネルが制度的に閉じている
      b_cog（認知的障壁） → ノイズによる意味の劣化
      α_align（方向不一致）→ 符号化・復号化のミスマッチ

    意味：エッジ通過時の情報損失率。
    原典対応：b_norm + b_cog + max(0, 1-α_align)
    証拠分類：Derived（同時観測データがある場合）/ Assumed（推定の場合）
    """
    p_A, p_B = _shared_bins_distribution(data_A, data_B, bins)
    p_A = np.clip(p_A, EPS, None)
    p_B = np.clip(p_B, EPS, None)
    H_A = float(-np.sum(p_A * np.log(p_A)))
    H_B = float(-np.sum(p_B * np.log(p_B)))

    if len(data_A) == len(data_B):
        a = np.asarray(data_A, dtype=float)
        b = np.asarray(data_B, dtype=float)
        corr = float(np.corrcoef(a, b)[0, 1]) if len(a) > 1 else 0.0
        I_AB = max(0.0, -0.5 * math.log(max(1 - corr ** 2, EPS)))
        ev = Evidence.DERIVED
    else:
        I_AB = 0.0
        ev = Evidence.ASSUMED

    T5 = float(max(0.0, min(1.0, 1.0 - I_AB / (H_A + EPS))))
    return T5, ev


def calc_regime_tension_vector(
    data_A: List[float],
    data_B: List[float],
    alpha_T3: float = 1.0,
    bins: int = 50
) -> RegimeTensionVector:
    """
    𝓣 = (T₁, T₂, T₃, T₄, T₅) ∈ ℝ⁵ を計算して返す。

    設計原則：
      ・ハイパーパラメータ（重み係数）を持たない
      ・スカラーに集約しない
      ・各成分を独立した次元として保持する
      ・加算・重み付け合成を行わない

    Parameters
    ----------
    data_A   : レジームA の時系列データ
    data_B   : レジームB の時系列データ
    alpha_T3 : T₃ の α 値（α≠0 必須）
    bins     : 経験分布推定のビン数
    """
    assert abs(alpha_T3) > 1e-6, "alpha_T3 ≠ 0 が必須です"
    T1, ev1 = calc_T1(data_A, data_B, bins)
    T2, ev2 = calc_T2(data_A, data_B, bins)
    T3, ev3 = calc_T3(data_A, data_B, alpha_T3, bins)
    T4, ev4 = calc_T4(data_A, data_B, alpha_T3)
    T5, ev5 = calc_T5(data_A, data_B, bins)
    return RegimeTensionVector(
        T1=T1, T2=T2, T3=T3, T4=T4, T5=T5,
        T1_evidence=ev1, T2_evidence=ev2, T3_evidence=ev3,
        T4_evidence=ev4, T5_evidence=ev5,
        T3_alpha=alpha_T3
    )


# ============================================================
# Layer 2: グランドポテンシャル Φ
# ============================================================

def calc_grand_potential(F: float, mu: float, N: float) -> float:
    """
    グランドポテンシャル Φ を計算する。

    Φ = F - μ·N  [J]
      = tr(G)/d - ln2·H_bit·tr(G⁻¹)/d - 2√T₂·N

    理論的背景：
      開放系（粒子数 N が変化できる系）のポテンシャル。
      化学ポテンシャル μ = ∂F/∂N が粒子（人）の移動コストを表す。
      μ_A > μ_B のとき、人は A から B へ自発的に移動する。

    放置すれば減少する理由：
      第1項 -S·dT < 0：散逸により T が上昇
      第2項 -P·dV < 0：P < 0（散逸時）かつ dV < 0（散逸時）
      第3項 -N·dμ < 0：放置により T₂ が増大し μ が上昇

    証拠分類：Derived（F が Derived かつ N が Observed の場合）

    Parameters
    ----------
    F  : 自由エネルギー [J]
    mu : 化学ポテンシャル = 2√T₂ [J/人]
    N  : 移動する人数 [人]

    Returns
    -------
    Phi : float  [J]
    """
    return float(F - mu * N)


def calc_mu(T2: float) -> float:
    """
    化学ポテンシャル μ を計算する。

    μ = 2√T₂  [J/人]

    理論的背景：
      Fisher-Rao 測地距離 d_FR ≈ 2√T₂（Bhattacharyya 近似）
      μ は「1人がレジーム間を移動するコスト」に対応する。
      T₂ が小さい（2レジームが近い）ほど移動コストが低く、
      人が移動しやすい状態を表す。

    証拠分類：Derived（T₂ が Derived の場合）
    """
    return float(2.0 * math.sqrt(max(T2, 0.0)))


def calc_norm_dPhi_G(
    grad_Phi: np.ndarray,
    G_inv: np.ndarray
) -> float:
    """
    Φ の Fisher 計量でのノルム ‖dΦ‖_G を計算する。

    ‖dΦ‖_G = √(∇Φᵀ · G⁻¹ · ∇Φ)

    これは P の大きさ（絶対値）に対応する：
      P = sign(dΦ) · ‖dΦ‖_G

    証拠分類：Derived（G と ∇Φ が Derived の場合）

    Parameters
    ----------
    grad_Phi : Φ の θ に対する勾配 ∇_θΦ
    G_inv    : Fisher 行列の逆行列 G⁻¹

    Returns
    -------
    norm_dPhi_G : float  [J]
    """
    val = float(grad_Phi @ G_inv @ grad_Phi)
    return float(math.sqrt(max(val, 0.0)))


# ============================================================
# Layer 3: 位相判定
# ============================================================

def classify_phase(S_val: float, D_val: float, H_val: float) -> Tuple[str, str, str]:
    """
    S/D/H から現在の位相を判定する。

    5相の定義：
      Living System  ：中S（≤0.6）・中D（1.2〜1.7）・中〜高H（≥0.45）
                       動的平衡・健全な自己修復能力を保持
      Frozen Order   ：低S（≤0.4）・低D（≤1.4）・高H（≥0.55）
                       安定しているが変化が起きにくい硬直状態
      Runaway Growth ：中S・高D（≥1.6）・高H（≥0.6）
                       過成長・バブルの兆候
      Noise Dominant ：高S（≥0.65）・高D（≥1.5）・低H（≤0.45）
                       構造崩壊の前兆
      Collapse       ：高S（≥0.8）・不安定D・低H（≤0.3）
                       持続不能な状態

    Returns
    -------
    (phase_label, phase_css, explanation) : Tuple[str, str, str]
    """
    if S_val <= 0.6 and 1.2 <= D_val <= 1.7 and H_val >= 0.45:
        return (
            "Living System", "phase-living",
            "動的平衡状態。健全な自己修復能力を保持しています。"
        )
    if S_val <= 0.4 and D_val <= 1.4 and H_val >= 0.55:
        return (
            "Frozen Order", "phase-frozen",
            "安定しているが変化が起きにくい硬直状態です。"
        )
    if D_val >= 1.6 and H_val >= 0.6:
        return (
            "Runaway Growth", "phase-runaway",
            "過成長・バブルの兆候があります。急激な変化に注意が必要です。"
        )
    if S_val >= 0.8 and H_val <= 0.3:
        return (
            "Collapse", "phase-collapse",
            "持続不能な状態です。緊急の介入が必要です。"
        )
    if S_val >= 0.65 and H_val <= 0.45:
        return (
            "Noise Dominant", "phase-noise",
            "構造崩壊の前兆が見られます。早急な対応が必要です。"
        )
    return (
        "Living System", "phase-living",
        "現在のデータでは明確な位相判定が困難です。追加データを収集してください。"
    )


def detect_phase_transition(
    S_series: List[float],
    D_series: List[float],
    H_series: List[float],
    h_drop: float = -0.03,
    d_var_thresh: float = 0.01,
    s_rise: float = 0.03,
    lookback: int = 3
) -> Tuple[bool, str]:
    """
    相転移予兆を検出する。

    3条件が全て揃った場合のみアラートを出す（単独条件での断定禁止）：
      条件1：dH/dt < 0  （H の低下 = 持続性が失われ始めた）
      条件2：Var(D) 増大 （D の揺らぎ = 構造が不安定化）
      条件3：dS/dt > 0  （S の上昇 = 散逸が進行中）

    理論的背景：
      H の低下が最も早い予兆（先行指標）。
      S の上昇は崩壊が可視化された後（遅行指標）。
      3条件の同時成立で「強いアラート」を発する。

    Returns
    -------
    (alert, alert_text) : Tuple[bool, str]
    """
    if len(H_series) < lookback + 1:
        return False, ""
    dH = H_series[-1] - H_series[-lookback - 1]
    D_recent = D_series[-lookback:]
    D_var = statistics.variance(D_recent) if len(D_recent) >= 2 else 0.0
    dS = S_series[-1] - S_series[-lookback - 1]
    if dH < h_drop and D_var > d_var_thresh and dS > s_rise:
        return True, (
            f"⚠️ 相転移予兆：H低下（Δ{dH:.3f}）・"
            f"D不安定（Var={D_var:.4f}）・"
            f"S上昇（Δ{dS:.3f}）の3条件が同時成立。"
        )
    return False, ""


# ============================================================
# Layer 3: KPI ステータス判定
# ============================================================

def get_kpi_status(variable: str, value: float) -> Tuple[str, str]:
    """
    KPI 変数の値からステータス CSS クラスと説明テキストを返す。

    Returns
    -------
    (css_class, status_text) : Tuple[str, str]
    css_class ∈ {good, warn, alert, mid}
    """
    rules = {
        'S': [
            (0.8, 'alert', '高 ⚠'),
            (0.6, 'warn',  '中高 △'),
            (0.4, 'mid',   '中'),
            (0.0, 'good',  '低 ✓'),
        ],
        'D': [
            (1.8, 'alert', '高複雑 ⚠'),
            (1.6, 'warn',  '中高 △'),
            (1.4, 'mid',   '中'),
            (0.0, 'good',  '安定 ✓'),
        ],
        'H': [
            (0.6,  'good',  '高持続 ✓'),
            (0.45, 'mid',   '中'),
            (0.3,  'warn',  '低下中 △'),
            (0.0,  'alert', '低 ⚠'),
        ],
        'T1': [(0.5, 'alert', '高 ⚠'), (0.2, 'warn', '中 △'), (0.0, 'good', '低 ✓')],
        'T2': [(1.0, 'alert', '遠 ⚠'), (0.3, 'warn', '中 △'), (0.0, 'good', '近 ✓')],
        'T3': [(0.5, 'alert', '非対称大 ⚠'), (0.1, 'warn', '中 △'), (0.0, 'good', '低 ✓')],
        'T4': [(0.5, 'alert', '曲率大 ⚠'), (0.2, 'warn', '中 △'), (0.0, 'good', '整合 ✓')],
        'T5': [(0.7, 'alert', '損失大 ⚠'), (0.3, 'warn', '中 △'), (0.0, 'good', '低損失 ✓')],
    }
    # H は高い方が good（逆順）
    if variable == 'H':
        for threshold, css, text in rules['H']:
            if value >= threshold:
                return css, text
        return 'alert', '低 ⚠'
    for threshold, css, text in rules.get(variable, [(0.0, 'mid', '—')]):
        if value >= threshold:
            return css, text
    return 'mid', '—'


# ============================================================
# Layer 3: 介入最適化
# ============================================================

def calc_intervention_priority(
    Phi_current: float,
    Phi_prev: float,
    dt: float,
    tau: float,
    Phi_threshold: float,
    T2: float,
    N: float,
    cost_N: float = 1.0,
    cost_T2: float = 1.0,
    cost_theta: float = 10.0,
    G: Optional[np.ndarray] = None,
    H_bit: float = 0.0
) -> InterventionResult:
    """
    グランドポテンシャル Φ に基づく介入優先順位を計算する。

    アルゴリズム（設計方針）：
      Step 1: dΦ/dt を推定
      Step 2: Φ(t+τ) を予測
      Step 3: 介入必要性の判定
      Step 4: 各操作変数の感度計算
              ∂Φ/∂N  = -μ = -2√T₂
              ∂Φ/∂T₂ = -N/√T₂
              ∂Φ/∂θ  ≈ -S·∂T/∂θ（近似）
      Step 5: priority_i = |∂Φ/∂u_i| / c_i
      Step 6: 最小インパクト Δu* = ΔΦ_need / (∂Φ/∂u*)

    操作変数：
      u₁ = N  （移動する人数）         直接操作可能
      u₂ = T₂ （Bhattacharyya 距離）   直接操作可能
      u₃ = θ  （状態パラメータ）        間接操作

    ⚠️ v0.0.4 では設計方針の定式化のみ。
       Python による実際の最適化ループは v0.1 以降で実装予定。

    Parameters
    ----------
    Phi_current   : 現在の Φ [J]
    Phi_prev      : 前時点の Φ [J]
    dt            : 時間ステップ [任意単位]
    tau           : 予測ホライズン [任意単位]
    Phi_threshold : 介入閾値 Φ₀ [J]
    T2            : 現在の T₂ [J]
    N             : 現在の人数 [人]
    cost_N        : u₁=N の介入コスト [Assumed]
    cost_T2       : u₂=T₂ の介入コスト [Assumed]
    cost_theta    : u₃=θ の介入コスト [Assumed]
    G             : Fisher 情報行列（θ 方向の感度計算に使用）
    H_bit         : シャノンエントロピー（θ 方向の感度計算に使用）

    Returns
    -------
    InterventionResult
    """
    # Step 1：dΦ/dt を推定
    dPhi_dt = (Phi_current - Phi_prev) / (dt + EPS)

    # Step 2：Φ(t+τ) を予測
    Phi_predicted = Phi_current + dPhi_dt * tau

    # Step 3：介入必要性の判定
    needs_intervention = Phi_predicted < Phi_threshold
    delta_phi_need = max(0.0, Phi_threshold - Phi_predicted)

    # Step 4：各操作変数の感度計算
    mu = calc_mu(T2)
    # ∂Φ/∂N = -μ = -2√T₂
    sens_N  = -mu
    # ∂Φ/∂T₂ = -N/√T₂
    sens_T2 = -N / (math.sqrt(max(T2, EPS)))
    # ∂Φ/∂θ ≈ -S·∂T/∂θ（近似：Fisher 行列の変化に比例）
    if G is not None and H_bit > 0:
        d = G.shape[0]
        # 近似：∂T/∂θ ≈ tr(G⁻²·∂G)/d·k_B（ここでは単純スカラー近似）
        try:
            G_inv = np.linalg.inv(G + np.eye(d) * EPS)
            approx_dT_dtheta = float(np.trace(G_inv @ G_inv)) / (d * K_B)
        except np.linalg.LinAlgError:
            approx_dT_dtheta = 0.0
        S_val = K_B * LN2 * H_bit
        sens_theta = -S_val * approx_dT_dtheta
    else:
        sens_theta = 0.0

    # Step 5：priority_i = |∂Φ/∂u_i| / c_i
    eff_N     = abs(sens_N)  / (cost_N  + EPS)
    eff_T2    = abs(sens_T2) / (cost_T2 + EPS)
    eff_theta = abs(sens_theta) / (cost_theta + EPS)

    # Step 6：最高効率の操作変数と必要操作量
    efficiencies = {
        'N':     (eff_N,     sens_N),
        'T2':    (eff_T2,    sens_T2),
        'theta': (eff_theta, sens_theta),
    }
    best_var = max(efficiencies, key=lambda k: efficiencies[k][0])
    best_sens = efficiencies[best_var][1]
    delta_u = delta_phi_need / (best_sens + EPS) if abs(best_sens) > EPS else 0.0

    return InterventionResult(
        phi_current=Phi_current,
        phi_predicted=Phi_predicted,
        phi_threshold=Phi_threshold,
        needs_intervention=needs_intervention,
        delta_phi_need=delta_phi_need,
        sensitivity_N=sens_N,
        sensitivity_T2=sens_T2,
        sensitivity_theta=sens_theta,
        efficiency_N=eff_N,
        efficiency_T2=eff_T2,
        efficiency_theta=eff_theta,
        best_variable=best_var,
        delta_u_optimal=delta_u,
        priority_N=eff_N,
        priority_T2=eff_T2,
        priority_theta=eff_theta,
    )


# ============================================================
# SVG 座標計算（レポート生成用）
# ============================================================

def calc_sd_position(S_val: float, D_val: float) -> Tuple[int, int, int]:
    """S-D 位相空間図の現在位置座標を返す。"""
    SD_X = int(40 + max(0.0, min(1.0, S_val)) * 240)
    SD_Y = int(240 - (max(1.0, min(2.0, D_val)) - 1.0) * 220)
    return SD_X, SD_Y, SD_Y - 14


def calc_h_polyline(H_series: List[float]) -> Tuple[str, int, int, int]:
    """H 時系列ポリラインの座標文字列と現在値マーカー座標を返す。"""
    n = len(H_series)
    if n == 0:
        return "40,130", 40, 130, 116
    points = []
    for i, h in enumerate(H_series):
        x = int(40 + i * 240 / max(n - 1, 1))
        y = int(240 - max(0.0, min(1.0, h)) * 220)
        points.append(f"{x},{y}")
    H_CURRENT_X = 280
    H_CURRENT_Y = int(240 - max(0.0, min(1.0, H_series[-1])) * 220)
    return " ".join(points), H_CURRENT_X, H_CURRENT_Y, H_CURRENT_Y - 10


# ============================================================
# 統合解析関数
# ============================================================

def analyze_single_regime(
    x: List[float],
    N: float = 0.0,
    T2: float = 0.0,
    bins_entropy: int = 32,
    bins_tension: int = 50,
) -> SDFTState:
    """
    単一レジームの時系列データから全 SDFT 変数を計算する。

    Parameters
    ----------
    x              : 時系列データ
    N              : 移動する人数（Observed）
    T2             : Bhattacharyya 距離（別途計算して渡す）
    bins_entropy   : エントロピー計算のビン数
    bins_tension   : 𝓣 計算のビン数

    Returns
    -------
    SDFTState
    """
    state = SDFTState()

    # S/D/H
    state.H_bit = calc_H_bit(x, bins_entropy)
    state.S     = calc_S(x, bins_entropy)
    state.D     = calc_D(x)
    state.H     = calc_H(x)
    state.evidence_S = Evidence.DERIVED
    state.evidence_D = Evidence.DERIVED
    state.evidence_H = Evidence.DERIVED

    # Fisher 行列と導出変数
    G = calc_fisher_matrix(x)
    state.G     = G
    state.G_inv = np.linalg.inv(G + np.eye(G.shape[0]) * EPS)
    state.T     = calc_T_temperature(G)
    state.U     = calc_U_internal_energy(G)
    state.F     = calc_F_free_energy(G, state.H_bit)
    state.V     = calc_V_volume(G)
    state.evidence_T = Evidence.DERIVED
    state.evidence_U = Evidence.DERIVED
    state.evidence_F = Evidence.DERIVED
    state.evidence_V = Evidence.DERIVED

    # グランドポテンシャル
    state.N   = N
    state.mu  = calc_mu(T2)
    state.Phi = calc_grand_potential(state.F, state.mu, state.N)
    state.evidence_N   = Evidence.OBSERVED
    state.evidence_mu  = Evidence.DERIVED
    state.evidence_Phi = Evidence.DERIVED

    return state


# ============================================================
# メイン実行例
# ============================================================

if __name__ == "__main__":
    import random
    random.seed(42)

    print("=" * 60)
    print("SDFT revised v0.0.4 参照実装 — 動作確認")
    print("=" * 60)

    # サンプル時系列データ
    n = 300
    regime_A = [500 + random.gauss(0, 50) + i * 0.3 for i in range(n)]
    regime_B = [350 + random.gauss(0, 80) + i * 0.2 for i in range(n)]

    print("\n--- Layer 1: 構造場変数 ---")
    H_bit = calc_H_bit(regime_A)
    S     = calc_S(regime_A)
    D     = calc_D(regime_A)
    H     = calc_H(regime_A)
    print(f"H_bit = {H_bit:.4f}  [bit]   Derived")
    print(f"S     = {S:.4e} [J/K]  Derived")
    print(f"D     = {D:.4f}  [1-2]  Derived")
    print(f"H     = {H:.4f}  [0-1]  Derived")

    print("\n--- Layer 2: Fisher 行列由来 ---")
    G = calc_fisher_matrix(regime_A)
    T = calc_T_temperature(G)
    U = calc_U_internal_energy(G)
    F = calc_F_free_energy(G, H_bit)
    V = calc_V_volume(G)
    print(f"G     =\n{G}")
    print(f"T     = {T:.4e} [K]     Derived")
    print(f"U     = {U:.4f}  [J相当] Derived")
    print(f"F     = {F:.4f}  [J]     Derived")
    print(f"V     = {V:.4e} [m³相当] Derived")

    print("\n--- Layer 2: 𝓣 ベクトル ---")
    tension = calc_regime_tension_vector(regime_A, regime_B, alpha_T3=1.0)
    print(tension.summary())

    print("\n--- Layer 2: グランドポテンシャル Φ ---")
    N_people = 100.0
    mu_val   = calc_mu(tension.T2)
    Phi_val  = calc_grand_potential(F, mu_val, N_people)
    print(f"N   = {N_people:.1f}   [人]   Observed")
    print(f"μ   = {mu_val:.4f}  [J/人] Derived  (= 2√T₂)")
    print(f"Φ   = {Phi_val:.4f}  [J]    Derived")

    print("\n--- Layer 3: 位相判定 ---")
    label, css, explanation = classify_phase(S, D, H)
    print(f"位相: {label} ({css})")
    print(f"説明: {explanation}")

    print("\n--- Layer 3: 介入最適化（設計方針） ---")
    Phi_prev  = Phi_val * 0.98  # 仮の前時点（2%減少を仮定）
    result = calc_intervention_priority(
        Phi_current=Phi_val,
        Phi_prev=Phi_prev,
        dt=1.0,
        tau=5.0,
        Phi_threshold=Phi_val * 0.95,
        T2=tension.T2,
        N=N_people,
        G=G,
        H_bit=H_bit,
    )
    print(f"介入必要: {result.needs_intervention}")
    print(f"Φ 現在値: {result.phi_current:.4f}")
    print(f"Φ 予測値: {result.phi_predicted:.4f}")
    print(f"ΔΦ 必要量: {result.delta_phi_need:.4f}")
    print(f"最優先操作変数: {result.best_variable}")
    print(f"必要操作量: {result.delta_u_optimal:.4f}")
    print(f"\npriority_score:")
    print(f"  N  = {result.priority_N:.4f}  (感度 {result.sensitivity_N:.4f})")
    print(f"  T₂ = {result.priority_T2:.4f}  (感度 {result.sensitivity_T2:.4f})")
    print(f"  θ  = {result.priority_theta:.4f}  (感度 {result.sensitivity_theta:.4e})")

    print("\n" + "=" * 60)
    print("動作確認完了")
    print("=" * 60)
