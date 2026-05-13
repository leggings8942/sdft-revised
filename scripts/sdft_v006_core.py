"""
sdft_v006_core.py
SDFT revised v0.0.6 — コア参照実装

このモジュールは SDFT（Structural Dynamics Field Theory）の Layer 1 と
Layer 2 の基本実装を提供する。

【Layer 1：構造場変数】
  H_bit = -Σ p_i log2(p_i)        [bit]   Derived
  S     = k_B · ln2 · H_bit       [J/K]   Derived
  D     = higuchi_fd(x)            [1-2]   Derived
  H     = hurst(x)                 [0-1]   Derived

【Layer 2：Fisher 行列とその導出量】
  G  = Fisher 情報行列（標準化データから計算）  Derived
  T  = tr(G⁻¹)/(d·k_B)            [K]      Derived
  U  = tr(G)/d                    [J相当]  Derived
  F  = U - T·S                    [J]      Derived
  V  = √det(G) · θ_range^d        [m³相当] Derived

【Layer 2：レジーム間テンション 𝓣】
  𝓣 = (T₁, T₂, T₃, T₄, T₅) ∈ ℝ⁵
    T₁ = |ψ_A - ψ_B|                  対数分配関数差
    T₂ = -log BC(p_A, p_B)            Bhattacharyya 距離
    T₃ = D^(α≠0)(p_A‖p_B)            α-ダイバージェンス（α=0禁止）
    T₄ = |log det G_A - log det G_B|  統計的ホロノミー（近似）
    T₅ = 1 - I(X_A;X_B)/H(X_A)        通信路損失率

【Layer 3：グランドポテンシャル】
  μ  = 2√T₂                       [J/人]   Derived
  Φ  = F - μN                     [J]      Derived

設計原則：
  1. データを標準化してから G を計算する（スケール依存性の排除）
  2. T₃ の α は 0 以外を使う（T₂ との重複を避ける）
  3. T₄ は非標準化分散の対数行列式差で近似する
  4. 𝓣 をスカラーに集約しない
"""

from __future__ import annotations

import math
import statistics
from dataclasses import dataclass
from typing import List, Tuple, Optional, NamedTuple

import numpy as np

# ============================================================
# 物理定数・設定
# ============================================================

K_B: float = 1.380649e-23   # ボルツマン定数 [J/K]
LN2: float = math.log(2)
EPS: float = 1e-12


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
    T3_evidence: str = Evidence.DERIVED
    T4_evidence: str = Evidence.ASSUMED
    T5_evidence: str = Evidence.DERIVED
    T3_alpha: float = 1.0

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


# ============================================================
# Layer 1: S（シャノンエントロピー → 熱力学的エントロピー）
# ============================================================

def calc_H_bit(x: List[float], bins: int = 32) -> float:
    """
    シャノンエントロピー H_bit を計算する。

    H_bit = -Σ p_i log₂(p_i)  [bit]

    証拠分類：Derived
    """
    arr = np.asarray(x, dtype=float)
    if len(arr) < 2:
        return 0.0
    counts, _ = np.histogram(arr, bins=bins)
    total = counts.sum()
    if total == 0:
        return 0.0
    probs = counts[counts > 0] / total
    return float(-np.sum(probs * np.log2(probs + EPS)))


def calc_S(x: List[float], bins: int = 32) -> float:
    """
    熱力学的エントロピー S を計算する。

    S = k_B · ln(2) · H_bit  [J/K]

    証拠分類：Derived
    """
    return float(K_B * LN2 * calc_H_bit(x, bins))


# ============================================================
# Layer 1: D（フラクタル次元・ヒグチ法）
# ============================================================

def calc_D(x: List[float], k_max: int = 10) -> float:
    """
    ヒグチ法によるフラクタル次元 D を計算する。

    D ∈ [1.0, 2.0]
    証拠分類：Derived
    """
    arr = np.asarray(x, dtype=float)
    n = len(arr)
    if n < k_max + 2:
        return 1.0
    lk = []
    for k in range(1, k_max + 1):
        lengths = []
        for m in range(1, k + 1):
            idxs = np.arange(m - 1, n, k)
            if len(idxs) < 2:
                continue
            length = np.sum(np.abs(np.diff(arr[idxs])))
            norm   = (n - 1) / ((len(idxs) - 1) * k)
            lengths.append(float(length * norm))
        if lengths:
            lk.append((math.log(k + EPS), math.log(np.mean(lengths) + EPS)))
    if len(lk) < 2:
        return 1.0
    xs = [p[0] for p in lk]
    ys = [p[1] for p in lk]
    coeffs = np.polyfit(xs, ys, 1)
    fd = abs(float(coeffs[0]))
    return float(max(1.0, min(2.0, fd)))


# ============================================================
# Layer 1: H（Hurst 指数）
# ============================================================

def calc_H(x: List[float], max_lag: int = 20) -> float:
    """
    ラグ別標準偏差法による Hurst 指数 H を計算する。

    std(x_{t+lag} - x_t) ~ lag^H
    → log(std) = H * log(lag) + const
    → H = slope of log(std) vs log(lag)  [/2 は不要・std を使うため]

    H ∈ [0.0, 1.0]
      H > 0.5：持続性（トレンドが続く傾向）
      H = 0.5：ランダムウォーク
      H < 0.5：反持続性（反転しやすい）

    ⚠️ 短い系列（n<200）では推定精度が低い。ペンディング項目あり。

    証拠分類：Derived
    """
    arr = np.asarray(x, dtype=float)
    n   = len(arr)
    if n < max_lag + 2:
        return 0.5
    log_lags, log_stds = [], []
    for lag in range(2, min(max_lag, n // 2) + 1):
        diffs = arr[lag:] - arr[:-lag]
        std   = float(np.std(diffs))
        if std > EPS:
            log_lags.append(math.log(lag))
            log_stds.append(math.log(std))
    if len(log_lags) < 2:
        return 0.5
    coeffs = np.polyfit(log_lags, log_stds, 1)
    # std ~ lag^H なので slope = H（/2 は不要）
    return float(max(0.0, min(1.0, coeffs[0])))


# ============================================================
# Layer 2: Fisher 情報行列
# ============================================================

def calc_fisher_matrix(x: List[float], d: int = 2) -> np.ndarray:
    """
    時系列データから Fisher 情報行列 G を推定する。

    正規分布 N(μ, σ²) の場合（閉形式）：
      G = [[1/σ², 0    ],
           [0,    2/σ⁴ ]]

    ⚠️ 重要：データを標準化してから G を計算する。
       標準化しないと σ² がデータのスケールに依存し、
       T = tr(G⁻¹)/(d·k_B) が天文学的な値になる。
       標準化後は σ² ≈ 1 となり G ≈ diag(1, 2)、
       T·S の積が合理的なスケールになる。

    証拠分類：Derived（正規分布仮定のもとで）
    """
    arr = np.asarray(x, dtype=float)
    if len(arr) < 2:
        return np.eye(d) * EPS

    # データを標準化（平均0・標準偏差1）
    mean = float(arr.mean())
    std  = float(arr.std(ddof=1))
    if std < EPS:
        std = EPS
    x_std = (arr - mean) / std

    var_hat = float(np.var(x_std, ddof=1))
    if var_hat < EPS:
        var_hat = EPS

    # 正規分布の Fisher 情報行列（2×2）
    G = np.array([
        [1.0 / var_hat,           0.0           ],
        [0.0,           2.0 / (var_hat ** 2)    ]
    ])
    if d != 2:
        # d≠2 の場合は対角スケーリングで近似
        G = np.eye(d) * (np.trace(G) / 2.0)
    return G


def calc_T_temperature(G: np.ndarray) -> float:
    """
    情報的温度 T を計算する。

    T = tr(G⁻¹) / (d · k_B)  [K]

    Jaynes の最大エントロピー原理による対応：
      k_B · T ↔ tr(G⁻¹)/d（系の揺らぎのエネルギースケール）

    証拠分類：Derived
    """
    d = G.shape[0]
    try:
        G_inv = np.linalg.inv(G + np.eye(d) * EPS)
        return float(np.trace(G_inv) / (d * K_B))
    except np.linalg.LinAlgError:
        return 0.0


def calc_U_internal_energy(G: np.ndarray) -> float:
    """
    内部エネルギー U を計算する。

    U = tr(G) / d  [J相当]

    等分配則との対応：各自由度に k_BT/2 のエネルギーが分配される。

    証拠分類：Derived
    """
    return float(np.trace(G) / G.shape[0])


def calc_F_free_energy(G: np.ndarray, H_bit: float) -> float:
    """
    自由エネルギー F を計算する。

    F = U - T · S
      = tr(G)/d - ln2 · H_bit · tr(G⁻¹)/d
      [J]

    k_B がキャンセルされる：
      T · S = [tr(G⁻¹)/(d·k_B)] · [k_B·ln2·H_bit]
            = ln2 · H_bit · tr(G⁻¹)/d

    証拠分類：Derived
    """
    d = G.shape[0]
    U = calc_U_internal_energy(G)
    try:
        G_inv = np.linalg.inv(G + np.eye(d) * EPS)
        T_S   = LN2 * H_bit * float(np.trace(G_inv)) / d
    except np.linalg.LinAlgError:
        T_S = 0.0
    return float(U - T_S)


def calc_V_volume(G: np.ndarray, theta_range: float = 1.0) -> float:
    """
    統計多様体の体積 V を計算する（近似）。

    V ≈ √det(G) · theta_range^d  [m³相当・近似]

    theta_range：積分範囲の近似（Assumed）。

    証拠分類：Derived（theta_range は Assumed）
    """
    d = G.shape[0]
    try:
        det_G = max(float(np.linalg.det(G)), EPS)
        return float(math.sqrt(det_G) * (theta_range ** d))
    except Exception:
        return EPS


# ============================================================
# Layer 2: 𝓣 ベクトルの各成分
# ============================================================

def _empirical_distribution(
    data_A: List[float],
    data_B: List[float],
    bins: int = 50,
) -> Tuple[np.ndarray, np.ndarray]:
    """共通ビンで2つの経験分布を推定する（内部関数）。"""
    a = np.asarray(data_A, dtype=float)
    b = np.asarray(data_B, dtype=float)
    all_data  = np.concatenate([a, b])
    bin_edges = np.linspace(all_data.min() - EPS, all_data.max() + EPS, bins + 1)
    p_A, _ = np.histogram(a, bins=bin_edges)
    p_B, _ = np.histogram(b, bins=bin_edges)
    p_A = p_A.astype(float) / (p_A.sum() + EPS)
    p_B = p_B.astype(float) / (p_B.sum() + EPS)
    return p_A, p_B


def calc_T1(
    data_A: List[float],
    data_B: List[float],
    bins: int = 50,
) -> Tuple[float, str]:
    """
    T₁：対数分配関数差（スケール差）

    T₁ = |ψ_A - ψ_B|,  ψ ≈ -E[log p(x)]

    意味：2レジームのスケール・自由度の乖離。
    証拠分類：Derived
    """
    a = np.asarray(data_A, dtype=float)
    b = np.asarray(data_B, dtype=float)

    def _psi(arr: np.ndarray) -> float:
        p, _ = np.histogram(arr, bins=bins)
        p    = p.astype(float) / (p.sum() + EPS)
        p    = np.clip(p, EPS, None)
        # 各観測値に対応する確率を割り当て
        bin_edges = np.linspace(arr.min() - EPS, arr.max() + EPS, bins + 1)
        idx = np.clip(
            ((arr - arr.min()) / ((arr.max() - arr.min() + EPS) / bins)).astype(int),
            0, bins - 1
        )
        return float(-np.mean(np.log(p[idx])))

    T1 = abs(_psi(a) - _psi(b))
    ev = Evidence.DERIVED if (len(data_A) >= 10 and len(data_B) >= 10) else Evidence.ASSUMED
    return float(T1), ev


def calc_T2(
    data_A: List[float],
    data_B: List[float],
    bins: int = 50,
) -> Tuple[float, str]:
    """
    T₂：Bhattacharyya 距離（状態点の幾何学的ずれ）

    T₂ = d_B(p_A, p_B) = -log BC(p_A, p_B)
    BC = ∫√(p_A · p_B) dx

    化学ポテンシャルとの関係：μ = 2√T₂
    証拠分類：Derived
    """
    p_A, p_B = _empirical_distribution(data_A, data_B, bins)
    BC = float(np.sum(np.sqrt(p_A * p_B)))
    BC = max(BC, EPS)
    T2 = float(-math.log(BC))
    ev = Evidence.DERIVED if (len(data_A) >= 10 and len(data_B) >= 10) else Evidence.ASSUMED
    return T2, ev


def calc_T3(
    data_A: List[float],
    data_B: List[float],
    alpha: float = 1.0,
    bins: int = 50,
) -> Tuple[float, str]:
    """
    T₃：α-ダイバージェンス（非対称な意味的差）

    T₃ = D^(α)(p_A ‖ p_B),  α ≠ 0

    α=+1：KL(p_A‖p_B)（Aの視点からBの意外さ）
    α=-1：KL(p_B‖p_A)（Bの視点からAの意外さ）

    ⚠️ α=0 は禁止（T₂ との重複）
    証拠分類：Derived
    """
    if abs(alpha) < EPS:
        raise ValueError("T₃ の α=0 は禁止です。α=1 または α=-1 を使用してください。")

    p_A, p_B = _empirical_distribution(data_A, data_B, bins)
    p_A = np.clip(p_A, EPS, None)
    p_B = np.clip(p_B, EPS, None)

    if abs(alpha - 1.0) < 1e-6:
        T3 = float(np.sum(p_A * np.log(p_A / p_B)))
    elif abs(alpha + 1.0) < 1e-6:
        T3 = float(np.sum(p_B * np.log(p_B / p_A)))
    else:
        integral = float(np.sum(p_A ** ((1 - alpha) / 2) * p_B ** ((1 + alpha) / 2)))
        T3 = float(4.0 / (1 - alpha ** 2) * (1 - integral))

    T3 = max(0.0, T3)
    ev = Evidence.DERIVED if (len(data_A) >= 10 and len(data_B) >= 10) else Evidence.ASSUMED
    return T3, ev


def calc_T4(
    data_A: List[float],
    data_B: List[float],
    alpha: float = 1.0,
) -> Tuple[float, str]:
    """
    T₄：統計的ホロノミー（多様体の曲率差）

    近似実装：非標準化分散の対数行列式差
    T₄ = |log det G_A - log det G_B| / max(|log det G_A|, |log det G_B|)

    det G（非標準化）= 2/σ⁶ → log det G = log(2) - 3·log(σ²)

    ⚠️ 標準化データで計算すると G_A = G_B となり T₄ が常に 0 になる。
       ここでは元の分散を使って曲率差を推定する。

    証拠分類：Assumed（近似推定）
    """
    a = np.asarray(data_A, dtype=float)
    b = np.asarray(data_B, dtype=float)

    var_A = float(np.var(a, ddof=1))
    var_B = float(np.var(b, ddof=1))

    if var_A < EPS or var_B < EPS:
        return 0.0, Evidence.ASSUMED

    # det G = 2/σ⁶ → log det G = log(2) - 3·log(σ²)
    log_det_A = math.log(2.0) - 3.0 * math.log(var_A + EPS)
    log_det_B = math.log(2.0) - 3.0 * math.log(var_B + EPS)
    denom     = max(abs(log_det_A), abs(log_det_B), EPS)
    T4        = float(abs(log_det_A - log_det_B) / denom)
    return float(max(0.0, min(1.0, T4))), Evidence.ASSUMED


def calc_T5(
    data_A: List[float],
    data_B: List[float],
    bins: int = 50,
) -> Tuple[float, str]:
    """
    T₅：通信路損失率

    T₅ = L(e) = 1 - I(X_A; X_B) / H(X_A)
    I(X_A; X_B) = H(X_A) + H(X_B) - H(X_A, X_B)

    T₅ = 0：完全に伝わる（損失ゼロ）
    T₅ = 1：全く伝わらない

    証拠分類：Derived（同時観測時）/ Assumed（推定時）
    """
    p_A, p_B = _empirical_distribution(data_A, data_B, bins)
    p_A = np.clip(p_A, EPS, None)
    p_B = np.clip(p_B, EPS, None)
    H_A = float(-np.sum(p_A * np.log(p_A)))
    H_B = float(-np.sum(p_B * np.log(p_B)))

    if len(data_A) == len(data_B):
        a = np.asarray(data_A, dtype=float)
        b = np.asarray(data_B, dtype=float)
        if len(a) > 1:
            corr = float(np.corrcoef(a, b)[0, 1])
        else:
            corr = 0.0
        I_AB = max(0.0, -0.5 * math.log(max(1 - corr ** 2, EPS)))
        ev   = Evidence.DERIVED
    else:
        I_AB = 0.0
        ev   = Evidence.ASSUMED

    T5 = float(max(0.0, min(1.0, 1.0 - I_AB / (H_A + EPS))))
    return T5, ev


def calc_regime_tension_vector(
    data_A:  List[float],
    data_B:  List[float],
    alpha_T3: float = 1.0,
    bins:    int = 50,
) -> RegimeTensionVector:
    """
    𝓣 = (T₁, T₂, T₃, T₄, T₅) を計算する。

    設計原則：
      ・ハイパーパラメータなし
      ・スカラー集約しない
      ・各成分は独立した次元
    """
    if abs(alpha_T3) < EPS:
        raise ValueError("alpha_T3 ≠ 0 が必須です")
    T1, ev1 = calc_T1(data_A, data_B, bins)
    T2, ev2 = calc_T2(data_A, data_B, bins)
    T3, ev3 = calc_T3(data_A, data_B, alpha_T3, bins)
    T4, ev4 = calc_T4(data_A, data_B, alpha_T3)
    T5, ev5 = calc_T5(data_A, data_B, bins)
    return RegimeTensionVector(
        T1=T1, T2=T2, T3=T3, T4=T4, T5=T5,
        T1_evidence=ev1, T2_evidence=ev2,
        T3_evidence=ev3, T4_evidence=ev4,
        T5_evidence=ev5, T3_alpha=alpha_T3,
    )


# ============================================================
# Layer 2: グランドポテンシャル Φ
# ============================================================

def calc_mu(T2: float) -> float:
    """化学ポテンシャル μ = 2√T₂  [J/人]  Derived"""
    return float(2.0 * math.sqrt(max(T2, 0.0)))


def calc_grand_potential(F: float, mu: float, N: float) -> float:
    """グランドポテンシャル Φ = F - μN  [J]  Derived"""
    return float(F - mu * N)


def calc_norm_dPhi_G(
    grad_Phi: np.ndarray,
    G_inv:    np.ndarray,
) -> float:
    """
    ‖dΦ‖_G = √(∇Φᵀ · G⁻¹ · ∇Φ)  [J]  Derived
    P = sign(dΦ) · ‖dΦ‖_G
    """
    val = float(grad_Phi @ G_inv @ grad_Phi)
    return float(math.sqrt(max(val, 0.0)))


# ============================================================
# Layer 3: 位相判定
# ============================================================

def classify_phase(
    S_norm: float,
    D:      float,
    H:      float,
) -> Tuple[str, str, str]:
    """
    S_norm（0-1 正規化済み）・D・H から現在の位相を判定する。

    5相：Living System / Frozen Order / Runaway Growth /
         Noise Dominant / Collapse

    Returns: (phase_label, phase_css, explanation)
    """
    if S_norm <= 0.6 and 1.2 <= D <= 1.7 and H >= 0.45:
        return ("Living System",  "phase-living",
                "動的平衡状態。健全な自己修復能力を保持しています。")
    if S_norm <= 0.4 and D <= 1.4 and H >= 0.55:
        return ("Frozen Order",   "phase-frozen",
                "安定しているが変化が起きにくい硬直状態です。")
    if D >= 1.6 and H >= 0.6:
        return ("Runaway Growth", "phase-runaway",
                "過成長・バブルの兆候があります。")
    if S_norm >= 0.8 and H <= 0.3:
        return ("Collapse",       "phase-collapse",
                "持続不能な状態です。緊急の介入が必要です。")
    if S_norm >= 0.65 and H <= 0.45:
        return ("Noise Dominant", "phase-noise",
                "構造崩壊の前兆が見られます。")
    return ("Living System",  "phase-living",
            "現在のデータでは明確な位相判定が困難です。追加データを収集してください。")


def detect_phase_transition(
    S_series: List[float],
    D_series: List[float],
    H_series: List[float],
    h_drop:         float = -0.03,
    d_var_thresh:   float = 0.01,
    s_rise:         float = 0.03,
    lookback:       int   = 3,
) -> Tuple[bool, str]:
    """
    相転移予兆を検出する。

    3条件が全て揃った場合のみアラート：
      dH/dt < 0（H の低下）
      Var(D) 増大（D の不安定化）
      dS/dt > 0（S の上昇）
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
# KPI ステータス判定
# ============================================================

def get_kpi_status(variable: str, value: float) -> Tuple[str, str]:
    """KPI 変数の値からステータス CSS クラスと説明テキストを返す。"""
    if variable == 'H':
        if value >= 0.6:   return 'good', '高持続 ✓'
        if value >= 0.45:  return 'mid',  '中'
        if value >= 0.3:   return 'warn', '低下中 △'
        return 'alert', '低 ⚠'

    rules = {
        'S':  [(0.8,'alert','高 ⚠'), (0.6,'warn','中高 △'),
               (0.4,'mid','中'), (0.0,'good','低 ✓')],
        'D':  [(1.8,'alert','高複雑 ⚠'), (1.6,'warn','中高 △'),
               (1.4,'mid','中'), (0.0,'good','安定 ✓')],
        'T1': [(0.5,'alert','高 ⚠'), (0.2,'warn','中 △'), (0.0,'good','低 ✓')],
        'T2': [(1.0,'alert','遠 ⚠'), (0.3,'warn','中 △'), (0.0,'good','近 ✓')],
        'T3': [(0.5,'alert','非対称大 ⚠'), (0.1,'warn','中 △'), (0.0,'good','低 ✓')],
        'T4': [(0.5,'alert','曲率大 ⚠'), (0.2,'warn','中 △'), (0.0,'good','整合 ✓')],
        'T5': [(0.7,'alert','損失大 ⚠'), (0.3,'warn','中 △'), (0.0,'good','低損失 ✓')],
    }
    for threshold, css, text in rules.get(variable, []):
        if value >= threshold:
            return css, text
    return 'mid', '—'


# ============================================================
# SVG 座標計算
# ============================================================

def calc_sd_position(
    S_norm: float,
    D:      float,
) -> Tuple[int, int, int]:
    """S-D 位相空間図の現在位置座標を返す。"""
    SD_X = int(40 + max(0.0, min(1.0, S_norm)) * 240)
    SD_Y = int(240 - (max(1.0, min(2.0, D)) - 1.0) * 220)
    return SD_X, SD_Y, SD_Y - 14


def calc_h_polyline(
    H_series: List[float],
) -> Tuple[str, int, int, int]:
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
    x:    List[float],
    x_B:  List[float],
    N:    float = 0.0,
    alpha_T3: float = 1.0,
    bins: int = 32,
) -> dict:
    """
    1つのレジームに対して全 SDFT 変数を計算する。

    Returns
    -------
    dict：全変数の値と証拠分類ラベル
    """
    H_bit = calc_H_bit(x, bins)
    S     = calc_S(x, bins)
    D     = calc_D(x)
    H     = calc_H(x)
    G     = calc_fisher_matrix(x)
    T     = calc_T_temperature(G)
    U     = calc_U_internal_energy(G)
    F     = calc_F_free_energy(G, H_bit)
    V     = calc_V_volume(G)
    ten   = calc_regime_tension_vector(x, x_B, alpha_T3=alpha_T3)
    mu    = calc_mu(ten.T2)
    Phi   = calc_grand_potential(F, mu, N)

    try:
        G_inv = np.linalg.inv(G + np.eye(G.shape[0]) * EPS)
    except Exception:
        G_inv = np.eye(G.shape[0])

    S_norm = H_bit / math.log2(bins + EPS)
    phase_label, phase_css, phase_expl = classify_phase(S_norm, D, H)

    return {
        'H_bit': H_bit, 'S': S, 'D': D, 'H': H,
        'G': G, 'G_inv': G_inv, 'T': T, 'U': U, 'F': F, 'V': V,
        'tension': ten,
        'N': N, 'mu': mu, 'Phi': Phi,
        'S_norm': S_norm,
        'phase_label': phase_label,
        'phase_css': phase_css,
        'phase_explanation': phase_expl,
    }


# ============================================================
# 動作確認
# ============================================================

if __name__ == "__main__":
    import random
    random.seed(42)
    n = 300

    # AR(0.7) — 持続傾向（H > 0.5 に近いはず）
    x_A = [0.0]
    for _ in range(n - 1):
        x_A.append(0.7 * x_A[-1] + random.gauss(0, 1))

    # AR(-0.6) — 反転傾向（H < 0.5）
    x_B = [0.0]
    for _ in range(n - 1):
        x_B.append(-0.6 * x_B[-1] + random.gauss(0, 1))

    print("=" * 60)
    print("SDFT v0.0.6 コア実装 — 動作確認")
    print("=" * 60)

    result = analyze_single_regime(x_A, x_B, N=100.0, alpha_T3=1.0)

    print(f"\n構造場変数:")
    print(f"  H_bit = {result['H_bit']:.4f} bit   [Derived]")
    print(f"  S     = {result['S']:.4e} J/K  [Derived]")
    print(f"  D     = {result['D']:.4f}      [Derived]")
    print(f"  H     = {result['H']:.4f}      [Derived]  (>0.5 で持続傾向)")

    print(f"\nFisher 行列由来:")
    print(f"  G     =\n{result['G']}")
    print(f"  T     = {result['T']:.4e} K   [Derived]")
    print(f"  U     = {result['U']:.4f}      [Derived]")
    print(f"  F     = {result['F']:.4f}      [Derived]")
    print(f"  S*T   = {result['S']*result['T']:.4f} J  （単位確認）")

    ten = result['tension']
    print(f"\n{ten.summary()}")

    print(f"\nグランドポテンシャル:")
    print(f"  μ = 2√T₂ = {result['mu']:.4f} J/人  [Derived]")
    print(f"  Φ = F-μN = {result['Phi']:.4f} J     [Derived]")

    print(f"\n位相判定: {result['phase_label']}")
    print(f"  {result['phase_explanation']}")

    # H の検証：AR+ > AR- になるか
    H_anti = calc_H(x_B)
    print(f"\nH 検証:")
    print(f"  AR(+0.7) H = {result['H']:.4f}  (期待 > {H_anti:.4f})")
    print(f"  AR(-0.6) H = {H_anti:.4f}")
    print(f"  AR+ > AR-: {'✅ OK' if result['H'] > H_anti else '❌ NG'}")
