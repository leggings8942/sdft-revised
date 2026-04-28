"""
sdft_revised_v003_reference.py
SDFT revised v0.0.3 参照実装

主要変更（v0.0.2 から）:
  - 𝓣 をスカラーから 5次元ベクトル RegimeTensionVector に変更
  - T₂：Fisher行列スペクトル差 → Bhattacharyya距離（Fisher-Rao近似）
  - T₃：α=0（JSD）禁止 → α≠0 に固定（KL または逆KL）
  - 重みパラメータ（κ,λ,η,μ）廃止

対応変数: S, D, H, F, 𝓣（5次元ベクトル）
削除済み: q, W, ε, μ, σ, v, Z, α_maxwell, τ, E, B, φ, H_eff, F_lorentz,
          A, P_behavior, 6成分, 障壁類
"""

import math
import statistics
from dataclasses import dataclass, field
from typing import List, Optional, Tuple, NamedTuple
import numpy as np


EPS = 1e-9


# ============================================================
# 𝓣 のベクトル型定義
# ============================================================

class RegimeTensionVector(NamedTuple):
    """
    レジーム間テンション 𝓣 の5次元ベクトル表現。
    集約・スカラー化は行わない。各成分を個別に保持する。
    """
    T1: float  # 対数分配関数差（スケール差）
    T2: float  # Bhattacharyya距離（状態点の幾何学的ずれ・Fisher-Rao近似）
    T3: float  # α-ダイバージェンス（非対称な情報的差・α≠0）
    T4: float  # 統計的ホロノミー（曲率差）
    T5: float  # 通信路損失率（旧社会層）

    # 証拠分類ラベル
    T1_evidence: str = "Assumed"
    T2_evidence: str = "Assumed"
    T3_evidence: str = "Assumed"
    T4_evidence: str = "Assumed"
    T5_evidence: str = "Assumed"

    # T₃ の α 値（α≠0 を強制）
    T3_alpha: float = 1.0

    def as_array(self) -> np.ndarray:
        """数値5成分をndarrayとして返す"""
        return np.array([self.T1, self.T2, self.T3, self.T4, self.T5])

    def summary(self) -> str:
        """レポート表示用サマリー"""
        return (
            f"𝓣 = (T₁={self.T1:.4f}[{self.T1_evidence}], "
            f"T₂={self.T2:.4f}[{self.T2_evidence}], "
            f"T₃={self.T3:.4f}[{self.T3_evidence},α={self.T3_alpha}], "
            f"T₄={self.T4:.4f}[{self.T4_evidence}], "
            f"T₅={self.T5:.4f}[{self.T5_evidence}])"
        )


# ============================================================
# Layer 0: 入力データ構造（v0.0.1 から継承・変更なし）
# ============================================================

@dataclass
class NodeState:
    node_id: str
    regime: str = "default"
    S: float = 0.0
    D: float = 1.0
    H: float = 0.5
    U: float = 1.0
    Theta: float = 0.5
    evidence: str = "Assumed"


@dataclass
class EdgeState:
    src: str
    dst: str
    weight: float = 1.0
    evidence: str = "Assumed"


# ============================================================
# Layer 1: S/D/H の計算（v0.0.1 から変更なし）
# ============================================================

def calc_entropy(x: List[float], bins: int = 32) -> float:
    """Shannon エントロピー（正規化済み）。証拠分類: Derived"""
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
    entropy = 0.0
    for c in counts:
        if c > 0:
            p = c / n
            entropy -= p * math.log(p + EPS)
    max_entropy = math.log(bins)
    return float(entropy / max_entropy) if max_entropy > 0 else 0.0


def higuchi_fd(x: List[float], k_max: int = 10) -> float:
    """ヒグチフラクタル次元。証拠分類: Derived"""
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
            length = sum(abs(x[idxs[i]] - x[idxs[i-1]]) for i in range(1, len(idxs)))
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


def hurst_exponent(x: List[float], max_lag: int = 20) -> float:
    """Hurst指数。証拠分類: Derived"""
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
# Layer 2: 𝓣 ベクトルの各成分計算（v0.0.3 新規）
# ============================================================

def _empirical_distribution(data: List[float], bins: int = 50) -> np.ndarray:
    """経験分布（正規化済みヒストグラム）を返す"""
    arr = np.asarray(data, dtype=float)
    all_data = arr
    hist, _ = np.histogram(arr, bins=bins)
    p = hist.astype(float)
    p = p / (p.sum() + EPS)
    return p


def _empirical_distribution_shared_bins(
    data_A: List[float],
    data_B: List[float],
    bins: int = 50
) -> Tuple[np.ndarray, np.ndarray]:
    """共通ビン幅で両レジームの経験分布を推定する"""
    a = np.asarray(data_A, dtype=float)
    b = np.asarray(data_B, dtype=float)
    all_data = np.concatenate([a, b])
    bin_edges = np.linspace(all_data.min() - EPS, all_data.max() + EPS, bins + 1)
    p_A, _ = np.histogram(a, bins=bin_edges)
    p_B, _ = np.histogram(b, bins=bin_edges)
    p_A = p_A.astype(float) / (p_A.sum() + EPS)
    p_B = p_B.astype(float) / (p_B.sum() + EPS)
    return p_A, p_B


def calc_T1(data_A: List[float], data_B: List[float]) -> Tuple[float, str]:
    """
    T₁：対数分配関数差（スケール差）
    ψ = -E[log p(x)] ≈ 負の対数尤度の経験平均
    証拠分類: Derived（両レジームの時系列が存在する場合）
    """
    psi_A = -np.mean(np.log(
        np.clip(_empirical_distribution(data_A, bins=50), EPS, None)[[
            min(int((v - min(data_A)) / ((max(data_A) - min(data_A) + EPS) / 50)), 49)
            for v in data_A
        ]]
    ))
    psi_B = -np.mean(np.log(
        np.clip(_empirical_distribution(data_B, bins=50), EPS, None)[[
            min(int((v - min(data_B)) / ((max(data_B) - min(data_B) + EPS) / 50)), 49)
            for v in data_B
        ]]
    ))
    T1 = float(abs(psi_A - psi_B))
    evidence = "Derived" if (len(data_A) >= 10 and len(data_B) >= 10) else "Assumed"
    return T1, evidence


def calc_T2(data_A: List[float], data_B: List[float], bins: int = 50) -> Tuple[float, str]:
    """
    T₂：Bhattacharyya距離（状態点の幾何学的ずれ・Fisher-Rao測地距離の近似）

    d_B(p_A, p_B) = -log ∫√(p_A(x)·p_B(x))dx = -log BC(p_A, p_B)

    Fisher-Rao との関係：d_FR ≈ 2√(d_B)
    分布仮定不要・双峰データにも対応。

    ⚠️ T₃（α-ダイバージェンス）との独立性：
    T₂ は対称（BC ベース）、T₃ は非対称（α≠0）。α≠0 を守る限り独立。

    証拠分類: Derived（経験分布から直接計算可能）
    """
    p_A, p_B = _empirical_distribution_shared_bins(data_A, data_B, bins)
    BC = float(np.sum(np.sqrt(p_A * p_B)))
    BC = max(BC, EPS)
    T2 = float(-np.log(BC))
    evidence = "Derived" if (len(data_A) >= 10 and len(data_B) >= 10) else "Assumed"
    return T2, evidence


def calc_T3(
    data_A: List[float],
    data_B: List[float],
    alpha: float = 1.0,
    bins: int = 50
) -> Tuple[float, str]:
    """
    T₃：α-ダイバージェンス（非対称な情報的差）

    D^(α)(p_A ‖ p_B) = 4/(1-α²) · (1 - ∫p_A^((1-α)/2) · p_B^((1+α)/2) dx)

    ⚠️ α = 0 は禁止（T₂ の Bhattacharyya と重複するため）
    推奨：α = 1（KL: A→B）または α = -1（逆KL: B→A）

    α の意味：
      α = +1：A の視点から B の「意外さ」（KL ダイバージェンス）
      α = -1：B の視点から A の「意外さ」（逆KL ダイバージェンス）

    証拠分類: Derived（経験分布から計算可能）
    """
    assert abs(alpha) > 1e-6, (
        "α = 0 は禁止です。T₂（Bhattacharyya）と重複します。"
        "α = 1（KL）または α = -1（逆KL）を使用してください。"
    )
    p_A, p_B = _empirical_distribution_shared_bins(data_A, data_B, bins)
    p_A = np.clip(p_A, EPS, None)
    p_B = np.clip(p_B, EPS, None)

    if abs(alpha - 1.0) < 1e-6:
        # KL ダイバージェンス（α=1 の極限）
        T3 = float(np.sum(p_A * np.log(p_A / p_B)))
    elif abs(alpha + 1.0) < 1e-6:
        # 逆 KL（α=-1 の極限）
        T3 = float(np.sum(p_B * np.log(p_B / p_A)))
    else:
        integral = float(np.sum(p_A**((1-alpha)/2) * p_B**((1+alpha)/2)))
        T3 = float(4 / (1 - alpha**2) * (1 - integral))

    T3 = max(0.0, T3)
    evidence = "Derived" if (len(data_A) >= 10 and len(data_B) >= 10) else "Assumed"
    return T3, evidence


def calc_T4_approx(
    data_A: List[float],
    data_B: List[float],
    alpha: float = 1.0
) -> Tuple[float, str]:
    """
    T₄：統計的ホロノミー（曲率差）の近似推定

    Hol^(α)(γ) = P·exp(∮ Γ^(α)_{ijk} dθ^k)
    指数型分布族での曲率テンソル：R_{ijkl} = (1-α²)/4 · (g_il·g_jk - g_ik·g_jl)

    v0.0.3 では経験的 Fisher 行列の固有値から近似推定する。

    証拠分類: Assumed（近似推定のため）
    """
    # 経験的 Fisher 行列のスカラー近似（1次元の場合）
    a = np.asarray(data_A, dtype=float)
    b = np.asarray(data_B, dtype=float)
    # 分散の逆数を Fisher 情報の近似とする（正規分布近似）
    var_A = max(float(np.var(a)), EPS)
    var_B = max(float(np.var(b)), EPS)
    G_A = 1.0 / var_A  # スカラー Fisher 情報（正規分布）
    G_B = 1.0 / var_B
    # 曲率テンソルの近似：R ≈ (1-α²)/4 · G²
    R_A = (1 - alpha**2) / 4 * G_A**2
    R_B = (1 - alpha**2) / 4 * G_B**2
    # ホロノミーの近似：‖Hol - I‖ ≈ |R_A - R_B| / max(R_A, R_B)
    T4 = float(abs(R_A - R_B) / (max(R_A, R_B) + EPS))
    T4 = max(0.0, min(1.0, T4))  # [0,1] に正規化
    return T4, "Assumed"


def calc_T5(
    data_A: List[float],
    data_B: List[float],
    bins: int = 50
) -> Tuple[float, str]:
    """
    T₅：通信路損失率（旧社会層の統合）

    L(e) = 1 - I(X_A; X_B) / H(X_A)
    I(X_A; X_B) = H(X_A) + H(X_B) - H(X_A, X_B)

    L = 0：完全に伝わる
    L = 1：全く伝わらない

    旧社会層との対応：
      b_norm（規範的障壁）→ チャネルが制度的に閉じている
      b_cog（認知的障壁）→ ノイズによる意味の劣化
      α_align（方向不一致）→ 符号化・復号化のミスマッチ

    証拠分類: Derived（同時観測データがある場合）/ Assumed（推定の場合）
    """
    p_A, p_B = _empirical_distribution_shared_bins(data_A, data_B, bins)
    p_A = np.clip(p_A, EPS, None)
    p_B = np.clip(p_B, EPS, None)

    H_A = float(-np.sum(p_A * np.log(p_A)))
    H_B = float(-np.sum(p_B * np.log(p_B)))

    # 同時分布の近似（独立仮定の場合は I=0、実際はデータから推定が必要）
    # 簡易実装：A と B が同じ時点のペアデータとして扱えない場合は近似
    if len(data_A) == len(data_B):
        # ペアデータとして同時分布を推定
        a = np.asarray(data_A, dtype=float)
        b = np.asarray(data_B, dtype=float)
        H_AB = float(-np.sum(p_A * np.log(p_A)) + (-np.sum(p_B * np.log(p_B))))
        # 相互情報量の下限（独立時は0）
        corr = float(np.corrcoef(a, b)[0, 1]) if len(a) > 1 else 0.0
        # 相関から相互情報量を近似（ガウス近似）
        I_AB = max(0.0, -0.5 * math.log(max(1 - corr**2, EPS)))
        evidence = "Derived"
    else:
        # データ長が異なる場合は独立仮定
        I_AB = 0.0
        evidence = "Assumed"

    T5 = float(max(0.0, min(1.0, 1.0 - I_AB / (H_A + EPS))))
    return T5, evidence


# ============================================================
# Layer 2: 𝓣 ベクトルの統合計算
# ============================================================

def calc_regime_tension_vector(
    data_A: List[float],
    data_B: List[float],
    alpha_T3: float = 1.0,
    bins: int = 50
) -> RegimeTensionVector:
    """
    𝓣 = (T₁, T₂, T₃, T₄, T₅) ∈ ℝ⁵ を計算して返す。

    パラメータ:
      data_A: レジームA の時系列データ
      data_B: レジームB の時系列データ
      alpha_T3: T₃ の α 値（α≠0 必須。推奨: 1.0 または -1.0）
      bins: 経験分布推定のビン数

    注意:
      - 集約（スカラー化）は行わない
      - 各成分を個別に保持して返す
      - 重みパラメータ（κ,λ,η,μ）は使用しない
    """
    assert abs(alpha_T3) > 1e-6, "alpha_T3 ≠ 0 が必須です"

    T1, ev1 = calc_T1(data_A, data_B)
    T2, ev2 = calc_T2(data_A, data_B, bins=bins)
    T3, ev3 = calc_T3(data_A, data_B, alpha=alpha_T3, bins=bins)
    T4, ev4 = calc_T4_approx(data_A, data_B, alpha=alpha_T3)
    T5, ev5 = calc_T5(data_A, data_B, bins=bins)

    return RegimeTensionVector(
        T1=T1, T2=T2, T3=T3, T4=T4, T5=T5,
        T1_evidence=ev1,
        T2_evidence=ev2,
        T3_evidence=ev3,
        T4_evidence=ev4,
        T5_evidence=ev5,
        T3_alpha=alpha_T3
    )


# ============================================================
# Layer 3: 位相判定・相転移アラート（v0.0.1 から変更なし）
# ============================================================

def classify_phase(S: float, D: float, H: float) -> Tuple[str, str, str]:
    if S <= 0.6 and 1.2 <= D <= 1.7 and H >= 0.45:
        return ("Living System", "phase-living", "動的平衡状態。健全な自己修復能力を保持しています。")
    if S <= 0.4 and D <= 1.4 and H >= 0.55:
        return ("Frozen Order", "phase-frozen", "安定しているが変化が起きにくい硬直状態です。")
    if D >= 1.6 and H >= 0.6:
        return ("Runaway Growth", "phase-runaway", "過成長・バブルの兆候があります。")
    if S >= 0.8 and H <= 0.3:
        return ("Collapse", "phase-collapse", "持続不能な状態です。緊急の介入が必要です。")
    if S >= 0.65 and H <= 0.45:
        return ("Noise Dominant", "phase-noise", "構造崩壊の前兆が見られます。")
    return ("Living System", "phase-living", "現在のデータでは明確な位相判定が困難です。")


def detect_phase_transition(
    S_series: List[float], D_series: List[float], H_series: List[float],
    h_drop: float = -0.03, d_var_thresh: float = 0.01,
    s_rise: float = 0.03, lookback: int = 3
) -> Tuple[bool, str]:
    if len(H_series) < lookback + 1:
        return False, ""
    dH = H_series[-1] - H_series[-lookback - 1]
    D_var = statistics.variance(D_series[-lookback:]) if len(D_series[-lookback:]) >= 2 else 0.0
    dS = S_series[-1] - S_series[-lookback - 1]
    if dH < h_drop and D_var > d_var_thresh and dS > s_rise:
        return True, (
            f"⚠️ 相転移予兆：H低下（Δ{dH:.3f}）・D不安定（Var={D_var:.4f}）・"
            f"S上昇（Δ{dS:.3f}）の3条件が同時成立。"
        )
    return False, ""


def get_kpi_status(variable: str, value: float) -> Tuple[str, str]:
    rules = {
        'S': [(0.8,'alert','高 ⚠'),(0.6,'warn','中高 △'),(0.4,'mid','中'),(0.0,'good','低 ✓')],
        'D': [(1.8,'alert','高複雑 ⚠'),(1.6,'warn','中高 △'),(1.4,'mid','中'),(0.0,'good','安定 ✓')],
        'H': [(0.6,'good','高持続 ✓'),(0.45,'mid','中'),(0.3,'warn','低下中 △'),(0.0,'alert','低 ⚠')],
    }
    if variable == 'H':
        for threshold, css, text in rules.get('H', []):
            if value >= threshold:
                return css, text
        return 'alert', '低 ⚠'
    for threshold, css, text in rules.get(variable, []):
        if value >= threshold:
            return css, text
    return 'mid', '—'


def get_tension_component_status(value: float, component: str) -> Tuple[str, str]:
    """
    𝓣 の各成分のステータスを返す。
    スカラー集約は行わず、成分ごとに個別判定する。
    """
    thresholds = {
        'T1': [(0.5,'alert','高 ⚠'),(0.2,'warn','中 △'),(0.0,'good','低 ✓')],
        'T2': [(1.0,'alert','遠 ⚠'),(0.3,'warn','中 △'),(0.0,'good','近 ✓')],
        'T3': [(0.5,'alert','非対称大 ⚠'),(0.1,'warn','中 △'),(0.0,'good','低 ✓')],
        'T4': [(0.5,'alert','曲率大 ⚠'),(0.2,'warn','中 △'),(0.0,'good','整合 ✓')],
        'T5': [(0.7,'alert','損失大 ⚠'),(0.3,'warn','中 △'),(0.0,'good','低損失 ✓')],
    }
    for threshold, css, text in thresholds.get(component, [(0.0,'mid','—')]):
        if value >= threshold:
            return css, text
    return 'mid', '—'


# ============================================================
# SVG 座標計算（v0.0.1 から変更なし）
# ============================================================

def calc_sd_position(S: float, D: float) -> Tuple[int, int, int]:
    SD_X = int(40 + max(0.0, min(1.0, S)) * 240)
    SD_Y = int(240 - (max(1.0, min(2.0, D)) - 1.0) * 220)
    return SD_X, SD_Y, SD_Y - 14


def calc_h_polyline(H_series: List[float]) -> Tuple[str, int, int, int]:
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
# メイン実行例
# ============================================================

if __name__ == "__main__":
    import random
    random.seed(42)

    # サンプルデータ：2つのレジーム
    regime_A = [500 + random.gauss(0, 50) + i * 0.5 for i in range(200)]
    regime_B = [350 + random.gauss(0, 80) + i * 0.3 for i in range(200)]

    print("=== SDFT revised v0.0.3 参照実装 ===\n")

    # S/D/H 計算
    S = calc_entropy(regime_A)
    D = higuchi_fd(regime_A)
    H = hurst_exponent(regime_A)
    print(f"S = {S:.4f} [Derived]")
    print(f"D = {D:.4f} [Derived]")
    print(f"H = {H:.4f} [Derived]")
    print()

    # 𝓣 ベクトル計算（α=1: KL ダイバージェンス）
    tension = calc_regime_tension_vector(regime_A, regime_B, alpha_T3=1.0)
    print("=== 𝓣（レジーム間テンション・5次元ベクトル）===")
    print(tension.summary())
    print()
    print("各成分のステータス：")
    for comp, val in [('T1', tension.T1), ('T2', tension.T2),
                      ('T3', tension.T3), ('T4', tension.T4), ('T5', tension.T5)]:
        css, text = get_tension_component_status(val, comp)
        print(f"  {comp} = {val:.4f} → {text} ({css})")
    print()

    # ⚠️ 集約しない設計の確認
    print("※ 𝓣 はスカラーに集約しません。5成分を個別に表示します。")
    print()

    # 位相判定
    label, css, explanation = classify_phase(S, D, H)
    print(f"位相判定: {label} ({css})")
    print(f"説明: {explanation}")
