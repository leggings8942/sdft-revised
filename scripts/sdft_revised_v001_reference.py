"""
sdft_revised_v001_reference.py
SDFT revised v0.0.1 参照実装

対応変数: S, D, H, F, 𝓣（再定義中）
削除済み: q, W, ε, μ, σ, v, Z, α, τ, E, B, φ, H_eff, F_lorentz,
          A, P_behavior, 6成分, 障壁類
"""

import math
import statistics
from dataclasses import dataclass, field
from typing import List, Optional, Tuple


EPS = 1e-9


# ============================================================
# Layer 0: 入力データ構造
# ============================================================

@dataclass
class NodeState:
    """ノードの状態（v0.0.1: 構造場変数のみ）"""
    node_id: str
    regime: str = "default"
    # 構造場変数（時系列から算出、またはDataFrame列として入力）
    S: float = 0.0   # エントロピー
    D: float = 1.0   # フラクタル次元
    H: float = 0.5   # Hurst指数
    U: float = 1.0   # 内部エネルギー（Assumed）
    Theta: float = 0.5  # 温度様パラメータ（Assumed）
    # 𝓣は再定義中のためノードには持たせず、エッジに定義
    evidence: str = "Assumed"  # 証拠分類ラベル


@dataclass
class EdgeState:
    """エッジの状態（v0.0.1: 接続情報のみ）"""
    src: str
    dst: str
    weight: float = 1.0   # 接続強度
    regime_tension: float = 0.0  # 𝓣（Assumed・再定義中）
    evidence: str = "Assumed"


# ============================================================
# Layer 1: 直接計算される変数（時系列 → S/D/H）
# ============================================================

def calc_entropy(x: List[float], bins: int = 32) -> float:
    """
    Shannon エントロピーを計算する。
    証拠分類: Derived（時系列データが存在する場合）
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
    entropy = 0.0
    for c in counts:
        if c > 0:
            p = c / n
            entropy -= p * math.log(p + EPS)
    # 正規化（最大エントロピー log(bins) で割る）
    max_entropy = math.log(bins)
    return float(entropy / max_entropy) if max_entropy > 0 else 0.0


def higuchi_fd(x: List[float], k_max: int = 10) -> float:
    """
    ヒグチ法によるフラクタル次元を計算する。
    証拠分類: Derived（時系列データが存在する場合）
    返値範囲: 1.0 〜 2.0
    """
    n = len(x)
    if n < k_max + 2:
        return 1.0
    lk = []
    x_arr = x
    for k in range(1, k_max + 1):
        lmk = []
        for m in range(1, k + 1):
            idxs = list(range(m - 1, n, k))
            if len(idxs) < 2:
                continue
            length = sum(
                abs(x_arr[idxs[i]] - x_arr[idxs[i - 1]])
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
    slope = num / (den + EPS)
    # フラクタル次元は傾きの絶対値
    fd = abs(slope)
    return float(max(1.0, min(2.0, fd)))


def hurst_exponent(x: List[float], max_lag: int = 20) -> float:
    """
    Hurst指数を計算する。
    証拠分類: Derived（時系列データが存在する場合）
    返値範囲: 0.0 〜 1.0
    H > 0.5: 持続性（トレンドが続く）
    H < 0.5: 反持続性（反転しやすい）
    H ≈ 0.5: ランダムウォーク
    """
    n = len(x)
    if n < max_lag + 2:
        return 0.5
    lags = range(2, min(max_lag, n // 2) + 1)
    log_lags = []
    log_stds = []
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
    H = float(max(0.0, min(1.0, slope / 2.0)))
    return H


def calc_free_energy(U: float, Theta: float, S: float) -> Tuple[float, str]:
    """
    自由エネルギー F = U - Θ·S を計算する。
    証拠分類: U・Θ が Assumed なら Assumed、Hypothesis なら Hypothesis。
    """
    F = U - Theta * S
    return float(F), "Assumed"


# ============================================================
# Layer 1: 時系列解析（窓スライド）
# ============================================================

def analyze_sdh_series(
    x: List[float],
    window: int = 128,
    step: int = 8
) -> List[dict]:
    """
    時系列を窓でスライドして S/D/H を時系列として計算する。
    返値: [{'t': int, 'S': float, 'D': float, 'H': float}, ...]
    """
    results = []
    n = len(x)
    for start in range(0, n - window + 1, step):
        seg = x[start:start + window]
        results.append({
            't': start + window,
            'S': calc_entropy(seg),
            'D': higuchi_fd(seg),
            'H': hurst_exponent(seg),
            'evidence': 'Derived'
        })
    return results


# ============================================================
# Layer 3: 位相判定
# ============================================================

PHASE_RULES = [
    # (label, css_class, S_max, D_min, D_max, H_min)
    ("Living System",  "phase-living",  0.6, 1.2, 1.7, 0.45),
    ("Frozen Order",   "phase-frozen",  0.4, 1.0, 1.4, 0.55),
    ("Runaway Growth", "phase-runaway", 0.6, 1.5, 2.0, 0.6),
    ("Noise Dominant", "phase-noise",   1.0, 1.5, 2.0, 0.0),
    ("Collapse",       "phase-collapse",1.0, 1.0, 2.0, 0.0),
]


def classify_phase(S: float, D: float, H: float) -> Tuple[str, str, str]:
    """
    S/D/H から位相ラベル・CSSクラス・説明を返す。
    返値: (label, css_class, explanation)
    """
    # Living System
    if S <= 0.6 and 1.2 <= D <= 1.7 and H >= 0.45:
        return ("Living System", "phase-living",
                "動的平衡状態。健全な自己修復能力を保持しています。")
    # Frozen Order
    if S <= 0.4 and D <= 1.4 and H >= 0.55:
        return ("Frozen Order", "phase-frozen",
                "安定しているが変化が起きにくい硬直状態です。")
    # Runaway Growth
    if D >= 1.6 and H >= 0.6:
        return ("Runaway Growth", "phase-runaway",
                "過成長・バブルの兆候があります。急激な変化に注意が必要です。")
    # Collapse（最悪ケース）
    if S >= 0.8 and H <= 0.3:
        return ("Collapse", "phase-collapse",
                "持続不能な状態です。緊急の介入が必要です。")
    # Noise Dominant
    if S >= 0.65 and H <= 0.45:
        return ("Noise Dominant", "phase-noise",
                "構造崩壊の前兆が見られます。早急な対応が必要です。")
    # デフォルト
    return ("Living System", "phase-living",
            "現在のデータでは明確な位相判定が困難です。追加データを収集してください。")


def detect_phase_transition(
    S_series: List[float],
    D_series: List[float],
    H_series: List[float],
    h_drop_threshold: float = -0.03,
    d_var_threshold: float = 0.01,
    s_rise_threshold: float = 0.03,
    lookback: int = 3
) -> Tuple[bool, str]:
    """
    相転移予兆を検出する。
    3条件が全て揃った時のみ True を返す（1条件だけでの断定は禁止）。
    条件: (dH/dt < 0) AND (Var(D)増大) AND (dS/dt > 0)
    """
    if len(H_series) < lookback + 1:
        return False, ""

    # dH/dt: 直近 lookback 期間の H の変化
    dH = H_series[-1] - H_series[-lookback - 1]
    # Var(D): 直近 lookback 期間の D の分散
    D_recent = D_series[-lookback:]
    D_var = statistics.variance(D_recent) if len(D_recent) >= 2 else 0.0
    # dS/dt: 直近 lookback 期間の S の変化
    dS = S_series[-1] - S_series[-lookback - 1]

    cond_H = dH < h_drop_threshold
    cond_D = D_var > d_var_threshold
    cond_S = dS > s_rise_threshold

    if cond_H and cond_D and cond_S:
        alert_text = (
            f"⚠️ 相転移予兆を検出：H低下（Δ{dH:.3f}）・D不安定（Var={D_var:.4f}）・"
            f"S上昇（Δ{dS:.3f}）の3条件が同時に成立しています。"
        )
        return True, alert_text
    return False, ""


# ============================================================
# Layer 3: KPI ステータス判定
# ============================================================

def get_kpi_status(variable: str, value: float) -> Tuple[str, str]:
    """
    KPI の値から STATUS_CSS クラスと状態テキストを返す。
    返値: (css_class, status_text)
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
            # H は高い方が良い（持続性）
            (0.6, 'good',  '高持続 ✓'),
            (0.45,'mid',   '中'),
            (0.3, 'warn',  '低下中 △'),
            (0.0, 'alert', '低 ⚠'),
        ],
        'T': [
            (0.7, 'alert', '高テンション ⚠'),
            (0.4, 'warn',  '中 △'),
            (0.2, 'mid',   '低中'),
            (0.0, 'good',  '低 ✓'),
        ],
    }
    if variable == 'H':
        # H は逆順（高い方が good）
        for threshold, css, text in rules.get('H', []):
            if value >= threshold:
                return css, text
        return 'alert', '低 ⚠'
    for threshold, css, text in rules.get(variable, []):
        if value >= threshold:
            return css, text
    return 'mid', '—'


# ============================================================
# SVG 座標計算
# ============================================================

def calc_sd_position(S: float, D: float) -> Tuple[int, int, int]:
    """S-D位相空間図の現在位置座標を返す。"""
    SD_X = int(40 + max(0.0, min(1.0, S)) * 240)
    SD_Y = int(240 - (max(1.0, min(2.0, D)) - 1.0) * 220)
    SD_Y_LABEL = SD_Y - 14
    return SD_X, SD_Y, SD_Y_LABEL


def calc_h_polyline(H_series: List[float]) -> Tuple[str, int, int, int]:
    """H時系列ポリラインの座標文字列と現在値マーカー座標を返す。"""
    n = len(H_series)
    if n == 0:
        return "40,130", 40, 130, 116
    points = []
    for i, h in enumerate(H_series):
        x = int(40 + i * 240 / max(n - 1, 1))
        y = int(240 - max(0.0, min(1.0, h)) * 220)
        points.append(f"{x},{y}")
    polyline = " ".join(points)
    H_CURRENT_X = 280
    H_CURRENT_Y = int(240 - max(0.0, min(1.0, H_series[-1])) * 220)
    H_CURRENT_Y_LABEL = H_CURRENT_Y - 10
    return polyline, H_CURRENT_X, H_CURRENT_Y, H_CURRENT_Y_LABEL


# ============================================================
# メイン実行例
# ============================================================

if __name__ == "__main__":
    import random
    random.seed(42)

    # サンプル時系列データ（来店数の模擬データ）
    sample_data = [500 + random.gauss(0, 50) + i * 0.5 for i in range(200)]

    print("=== SDFT revised v0.0.1 参照実装 ===\n")

    # S/D/H 計算
    S = calc_entropy(sample_data)
    D = higuchi_fd(sample_data)
    H = hurst_exponent(sample_data)
    F, F_evidence = calc_free_energy(U=1.0, Theta=0.5, S=S)

    print(f"S（エントロピー）  = {S:.4f}  [Derived]")
    print(f"D（フラクタル次元）= {D:.4f}  [Derived]")
    print(f"H（Hurst指数）    = {H:.4f}  [Derived]")
    print(f"F（自由エネルギー）= {F:.4f}  [{F_evidence}]")
    print(f"𝓣（レジーム間テンション）= 再定義中  [Assumed]")
    print()

    # 位相判定
    label, css, explanation = classify_phase(S, D, H)
    print(f"位相判定: {label} ({css})")
    print(f"説明: {explanation}")
    print()

    # KPI ステータス
    for var, val in [('S', S), ('D', D), ('H', H)]:
        css_cls, text = get_kpi_status(var, val)
        print(f"KPI {var}: {val:.4f} → {text} ({css_cls})")
    print()

    # 時系列解析
    sdh_series = analyze_sdh_series(sample_data, window=64, step=16)
    print(f"時系列解析: {len(sdh_series)} 点")

    if len(sdh_series) >= 4:
        S_ser = [r['S'] for r in sdh_series]
        D_ser = [r['D'] for r in sdh_series]
        H_ser = [r['H'] for r in sdh_series]
        alert, alert_text = detect_phase_transition(S_ser, D_ser, H_ser)
        if alert:
            print(f"アラート: {alert_text}")
        else:
            print("相転移予兆: なし")
