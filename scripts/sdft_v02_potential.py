"""
sdft_v02_potential.py
SDFT revised v0.2 — グランドポテンシャル Φ と位相判定

【方向 A：分解構造の保存】

  Φ = F_geom - μ · N

    F_geom = λ₂ （グラフ Laplacian のスペクトルギャップ）  Derived
    μ      = T₂ （測地距離 = 化学ポテンシャル）           Derived
    N      = 移動する人数                                Observed

【放置時に Φ が減少することの確認】
  放置 → 構造が緩む → λ₂ 低下 → F_geom 減少
  放置 → レジーム乖離 → μ 増大 → μN 増大
  両項とも Φ を押し下げる → Φ は単調減少

【位相分類（5 相）】
  S_g（軌道体積エントロピー）、D（相関次元）、λ₁（リアプノフ指数）
  の 3 変数で判定する。

【相転移予兆検出】
  dλ₁/dt の急上昇 ∧ Var(D) の増大 ∧ dS_g/dt の上昇
  の 3 条件 AND で検出する（v0.0.6 と同じ構造）。
"""

from __future__ import annotations

import math
import statistics
from typing import List, Tuple

import numpy as np

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
# 化学ポテンシャル
# ============================================================

def calc_mu(T2: float) -> float:
    """
    化学ポテンシャル μ = T₂（測地距離）

    v0.0.6 では μ = 2√T₂ だったが、v0.2 では T₂ が既に
    測地距離なのでそのまま使用する。

    証拠分類: Derived
    """
    return float(max(T2, 0.0))


# ============================================================
# グランドポテンシャル
# ============================================================

def calc_grand_potential(
    F_geom: float,
    mu: float,
    N: float,
) -> float:
    """
    Φ = F_geom - μ · N

    F_geom = λ₂ (spectral gap)
    μ      = T₂ (geodesic distance)
    N      = 人数 (Observed)

    証拠分類: Derived (N が Observed、他は Derived)
    """
    return float(F_geom - mu * N)


# ============================================================
# 位相分類（5相）
# ============================================================

def _normalize_S_g(
    S_g: float,
    S_g_ref_low: float = -5.0,
    S_g_ref_high: float = 5.0,
) -> float:
    """
    S_g を [0, 1] に正規化する。

    S_g_ref_low/high は Assumed パラメータ。
    データ固有の範囲を使うか、固定値を使うかはドメイン依存。
    """
    rng = S_g_ref_high - S_g_ref_low
    if rng < EPS:
        return 0.5
    return float(max(0.0, min(1.0, (S_g - S_g_ref_low) / rng)))


def _normalize_D(
    D: float,
    d_embed: int = 3,
) -> float:
    """
    相関次元 D を埋め込み次元で正規化する。

    D_norm = D / d_embed ∈ [0, 1]
    """
    if d_embed < 1:
        return 0.5
    return float(max(0.0, min(1.0, D / d_embed)))


def classify_phase(
    S_g: float,
    D: float,
    lambda_1: float,
    d_embed: int = 3,
) -> Tuple[str, str, str]:
    """
    S_g（正規化済み）・D（正規化済み）・λ₁ から現在の位相を判定する。

    5 相の定義：

    | 位相             | S_norm | D_norm   | λ₁      | 意味                     |
    |------------------|--------|----------|---------|--------------------------|
    | Living System    | ≤ 0.6  | 0.3-0.7  | < 0     | 動的平衡・健全な自己修復 |
    | Frozen Order     | ≤ 0.3  | ≤ 0.3    | < -0.1  | 安定だが変化なし         |
    | Runaway Growth   | 中     | ≥ 0.7    | < 0.05  | 過成長の兆候             |
    | Noise Dominant   | ≥ 0.7  | ≥ 0.6    | ≥ 0     | 構造崩壊の前兆           |
    | Collapse         | ≥ 0.85 | 不安定    | > 0.1   | 持続不能                 |

    注意：S_norm は _normalize_S_g() で、D_norm は _normalize_D() で
    事前に正規化すること。

    Returns: (phase_label, phase_css, explanation)
    """
    S_norm = _normalize_S_g(S_g)
    D_norm = _normalize_D(D, d_embed)

    # Collapse（最も深刻 → 先に判定）
    if S_norm >= 0.85 and lambda_1 > 0.1:
        return ("Collapse", "phase-collapse",
                "持続不能な状態です。軌道が急速に発散しており、緊急の介入が必要です。")

    # Noise Dominant
    if S_norm >= 0.7 and D_norm >= 0.6 and lambda_1 >= 0:
        return ("Noise Dominant", "phase-noise",
                "構造崩壊の前兆が見られます。軌道がカオス的で予測不能です。")

    # Frozen Order（Living System より前に判定：D_norm ≤ 0.3 で区別）
    if S_norm <= 0.3 and D_norm <= 0.3 and lambda_1 < -0.1:
        return ("Frozen Order", "phase-frozen",
                "安定しているが変化が起きにくい硬直状態です。")

    # Living System
    if S_norm <= 0.6 and 0.3 <= D_norm <= 0.7 and lambda_1 < 0:
        return ("Living System", "phase-living",
                "動的平衡状態。健全な自己修復能力を保持しています。")

    # Runaway Growth
    if D_norm >= 0.7 and lambda_1 < 0.05:
        return ("Runaway Growth", "phase-runaway",
                "過成長・バブルの兆候があります。")

    # デフォルト
    return ("Living System", "phase-living",
            "現在のデータでは明確な位相判定が困難です。追加データを収集してください。")


# ============================================================
# 相転移予兆検出
# ============================================================

def detect_phase_transition(
    S_g_series:     List[float],
    D_series:       List[float],
    lambda_1_series: List[float],
    lam_rise:        float = 0.02,
    d_var_thresh:    float = 0.01,
    s_rise:          float = 0.3,
    lookback:        int = 3,
) -> Tuple[bool, str]:
    """
    相転移予兆を検出する。

    3 条件が全て揃った場合のみアラート：
      dλ₁/dt > 0（リアプノフ指数の上昇 → カオスへの接近）
      Var(D) 増大（相関次元の不安定化）
      dS_g/dt > 0（軌道体積の膨張）
    """
    if len(lambda_1_series) < lookback + 1:
        return False, ""

    d_lam = lambda_1_series[-1] - lambda_1_series[-lookback - 1]
    D_recent = D_series[-lookback:]
    D_var = statistics.variance(D_recent) if len(D_recent) >= 2 else 0.0
    d_S = S_g_series[-1] - S_g_series[-lookback - 1]

    if d_lam > lam_rise and D_var > d_var_thresh and d_S > s_rise:
        return True, (
            f"⚠️ 相転移予兆：λ₁上昇（Δ{d_lam:.3f}）・"
            f"D不安定（Var={D_var:.4f}）・"
            f"S_g上昇（Δ{d_S:.3f}）の3条件が同時成立。"
        )
    return False, ""


# ============================================================
# KPI ステータス判定
# ============================================================

def get_kpi_status(variable: str, value: float) -> Tuple[str, str]:
    """KPI 変数の値からステータス CSS クラスと説明テキストを返す。"""

    if variable == 'lambda_1':
        if value < -0.1:   return 'good', '安定 ✓'
        if value < 0.0:    return 'mid',  '中立'
        if value < 0.05:   return 'warn', '不安定化 △'
        return 'alert', 'カオス ⚠'

    if variable == 'S_g':
        # S_g は正規化済み値で判定
        if value >= 0.85:  return 'alert', '高散逸 ⚠'
        if value >= 0.65:  return 'warn', '散逸進行 △'
        if value >= 0.35:  return 'mid', '中'
        return 'good', '秩序 ✓'

    if variable == 'D':
        # D_norm で判定
        if value >= 0.8:   return 'alert', '高複雑 ⚠'
        if value >= 0.6:   return 'warn', '中高 △'
        if value >= 0.35:  return 'mid', '中'
        return 'good', '安定 ✓'

    rules = {
        'T1': [(0.5, 'alert', '高 ⚠'), (0.2, 'warn', '中 △'), (0.0, 'good', '低 ✓')],
        'T2': [(2.0, 'alert', '遠 ⚠'), (0.5, 'warn', '中 △'), (0.0, 'good', '近 ✓')],
        'T3': [(0.8, 'alert', '捻り大 ⚠'), (0.3, 'warn', '中 △'), (0.0, 'good', '整合 ✓')],
        'T4': [(0.5, 'alert', '曲率差大 ⚠'), (0.2, 'warn', '中 △'), (0.0, 'good', '整合 ✓')],
        'T5': [(0.7, 'alert', '歪み大 ⚠'), (0.3, 'warn', '中 △'), (0.0, 'good', '低歪み ✓')],
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
    D_norm: float,
) -> Tuple[int, int, int]:
    """S-D 位相空間図の現在位置座標を返す。"""
    SD_X = int(40 + max(0.0, min(1.0, S_norm)) * 240)
    SD_Y = int(240 - max(0.0, min(1.0, D_norm)) * 220)
    return SD_X, SD_Y, SD_Y - 14


def calc_lambda1_polyline(
    lambda_1_series: List[float],
    y_min: float = -0.5,
    y_max: float = 0.5,
) -> Tuple[str, int, int, int]:
    """λ₁ 時系列ポリラインの座標文字列と現在値マーカー座標を返す。"""
    n = len(lambda_1_series)
    if n == 0:
        return "40,130", 40, 130, 116
    rng = y_max - y_min
    if rng < EPS:
        rng = 1.0
    points = []
    for i, v in enumerate(lambda_1_series):
        x = int(40 + i * 240 / max(n - 1, 1))
        frac = (v - y_min) / rng
        y = int(240 - max(0.0, min(1.0, frac)) * 220)
        points.append(f"{x},{y}")
    cur_x = 280
    cur_frac = (lambda_1_series[-1] - y_min) / rng
    cur_y = int(240 - max(0.0, min(1.0, cur_frac)) * 220)
    return " ".join(points), cur_x, cur_y, cur_y - 10


# ============================================================
# 動作確認
# ============================================================

if __name__ == "__main__":
    print("=" * 60)
    print("SDFT v0.2 Potential — 動作確認")
    print("=" * 60)

    # Φ の計算
    F_geom = 0.15    # spectral gap
    mu     = 0.8     # geodesic distance
    N      = 100.0
    Phi = calc_grand_potential(F_geom, mu, N)
    print(f"\nΦ = {F_geom} - {mu} × {N} = {Phi:.4f}")
    print(f"  F_geom (λ₂) = {F_geom}")
    print(f"  μ (T₂)      = {mu}")
    print(f"  N            = {N}")

    # 位相判定
    test_cases = [
        ("Living System 期待",  -1.0, 0.5, -0.05, 3),
        ("Frozen Order 期待",   -3.0, 0.3, -0.2,  3),
        ("Noise Dominant 期待",  2.0, 2.5,  0.1,  3),
        ("Collapse 期待",        3.0, 2.8,  0.15, 3),
    ]
    print("\n--- 位相判定 ---")
    for name, S_g, D, lam, d_e in test_cases:
        label, css, expl = classify_phase(S_g, D, lam, d_e)
        print(f"  {name}: → {label}")

    # 相転移検出
    S_series   = [0.0, 0.1, 0.2, 0.5, 1.0]
    D_series   = [1.0, 1.1, 1.0, 1.5, 2.0]
    lam_series = [-0.1, -0.08, -0.05, 0.0, 0.05]
    alert, msg = detect_phase_transition(S_series, D_series, lam_series)
    print(f"\n--- 相転移検出 ---")
    print(f"  alert = {alert}")
    if msg:
        print(f"  {msg}")

    print("\n" + "=" * 60)
    print("動作確認完了")
