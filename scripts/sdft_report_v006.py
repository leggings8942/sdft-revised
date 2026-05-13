"""
sdft_report_v006.py
SDFT revised v0.0.6 — レポート生成パイプライン

このモジュールは SDFT の最終出力層を提供する：
  ① §4-3 エッジ障壁ヒートマップ
     flow_ease = (1 - μ/μ_max)·(1 - T₅) を計算しテーブル + SVG で表示
  ② Φ₀（介入閾値）の自動設定
     Living System 相の境界値から導出
  ③ τ（予測ホライズン）の自動設定
     H（Hurst 指数）から導出
  ④ レポート生成パイプライン
     ReportContext → プレースホルダー置換 → HTML 保存

使い方の例：
    from sdft_v006_core      import analyze_single_regime
    from sdft_intervention_v006 import (
        StateSnapshot, generate_state_sequence,
        scan_lambda_tradeoff, find_optimal_lambda,
        run_intervention_loop,
    )
    from sdft_report_v006 import (
        ReportContext, generate_report,
        auto_phi_threshold, auto_tau,
    )

    # 状態列の生成
    states     = generate_state_sequence(x_series, x_B_series, N_series)

    # 自動設定
    phi_0, _   = auto_phi_threshold(states[-1])
    tau, _     = auto_tau(states[-1].H)

    # 介入最適化
    tradeoff   = scan_lambda_tradeoff(states, phi_0)
    opt_lam, _ = find_optimal_lambda(tradeoff)
    result     = run_intervention_loop(states, phi_0, tau=tau, lambda_count=opt_lam)

    # ReportContext の生成
    ctx = ReportContext.from_state(
        state=states[-1],
        regime_data={"A": data_A, "B": data_B},
        edges=[("A", "B"), ("B", "A")],
        intervention_result=result,
        opt_lambda=opt_lam,
        H_series=[s.H for s in states],
        title="SDFT 分析レポート",
        subject="分析対象",
    )

    # HTML レポートの生成
    html = generate_report(
        ctx=ctx,
        template_path="assets/report_template_v006.html",
        output_path="/mnt/user-data/outputs/report.html",
    )
"""

from __future__ import annotations

import math
import datetime
from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Optional
import numpy as np

from sdft_v006_core import (
    K_B, LN2, EPS,
    calc_H_bit, calc_S, calc_D, calc_H,
    calc_fisher_matrix, calc_T_temperature,
    calc_U_internal_energy, calc_F_free_energy, calc_V_volume,
    calc_T2, calc_T5, calc_regime_tension_vector,
    calc_mu, calc_grand_potential,
    classify_phase, detect_phase_transition,
    get_kpi_status,
    calc_sd_position, calc_h_polyline,
)
from sdft_intervention_v006 import (
    ActionConditionDiagnosis,
    calc_action_condition,
    StateSnapshot,
    generate_state_sequence,
    run_intervention_loop,
    find_optimal_lambda,
    scan_lambda_tradeoff,
)
from sdft_phi_graph_v006 import (
    PhiGraphData,
    build_phi_graph_data,
    render_phi_energy_graph_svg,
    render_phi_intervention_table_html,
)


# ============================================================
# ① エッジ障壁ヒートマップ
# ============================================================

@dataclass
class EdgeFlowData:
    """
    1本のエッジの流れやすさデータ。
    新しい変数は使用せず T₂・T₅ の読み替えのみ。
    """
    src:        str    # 送り元レジーム
    dst:        str    # 送り先レジーム
    T2:         float  # Bhattacharyya距離
    T5:         float  # 通信路損失率
    mu:         float  # = 2√T₂（移動コスト）
    flow_ease:  float  # = (1 - μ/μ_max) × (1 - T₅)  ∈ [0,1]
    status:     str    # "easy" / "moderate" / "hard"
    label:      str    # 表示ラベル

    @classmethod
    def compute(
        cls,
        src: str,
        dst: str,
        data_src: List[float],
        data_dst: List[float],
        mu_max: float = 3.0,
        bins: int = 50,
    ) -> "EdgeFlowData":
        """
        2つのレジームのデータからエッジの flow_ease を計算する。

        flow_ease = (1 - μ/μ_max) × (1 - T₅)

        理論的根拠：
          (1 - μ/μ_max)：移動コストが低いほど流れやすい
                         μ = 2√T₂（Fisher-Rao 近似・化学ポテンシャル）
          (1 - T₅)     ：通信路損失が低いほど流れやすい（AND 条件）
          積            ：両条件が揃って初めて流れる（乗算の正当化）

        全て Derived。外部パラメータは mu_max のみ（正規化用・Assumed）。
        """
        T2_val, _ = calc_T2(data_src, data_dst, bins=bins)
        T5_val, _ = calc_T5(data_src, data_dst, bins=bins)
        mu_val    = calc_mu(T2_val)

        a1 = max(0.0, 1.0 - mu_val / (mu_max + EPS))
        a2 = max(0.0, 1.0 - T5_val)
        ease = float(a1 * a2)

        if ease >= 0.6:
            status = "easy"
            label  = f"流れやすい ({ease:.2f})"
        elif ease >= 0.3:
            status = "moderate"
            label  = f"中程度 ({ease:.2f})"
        else:
            status = "hard"
            label  = f"流れにくい ({ease:.2f})"

        return cls(
            src=src, dst=dst,
            T2=T2_val, T5=T5_val,
            mu=mu_val, flow_ease=ease,
            status=status, label=label,
        )


def calc_edge_flows(
    regime_data: Dict[str, List[float]],
    edges: List[Tuple[str, str]],
    mu_max: float = 3.0,
    bins: int = 50,
) -> List[EdgeFlowData]:
    """
    複数のエッジについて flow_ease を一括計算する。

    Parameters
    ----------
    regime_data : {'レジーム名': 時系列データ}
    edges       : [(src, dst), ...] エッジのリスト
    mu_max      : μ の最大値（正規化用・Assumed）

    Returns
    -------
    List[EdgeFlowData]
    """
    results = []
    for src, dst in edges:
        if src in regime_data and dst in regime_data:
            ef = EdgeFlowData.compute(
                src=src, dst=dst,
                data_src=regime_data[src],
                data_dst=regime_data[dst],
                mu_max=mu_max, bins=bins,
            )
            results.append(ef)
    return results


def render_edge_heatmap_svg(
    edge_flows: List[EdgeFlowData],
    width: int = 600,
    height: int = 300,
) -> str:
    """
    エッジ障壁ヒートマップを SVG 文字列として生成する。

    各エッジを横棒グラフで表示：
      緑（flow_ease 高）→ 黄 → 赤（flow_ease 低）

    Returns
    -------
    SVG 文字列（<svg>...</svg>）
    """
    if not edge_flows:
        return (
            f'<svg viewBox="0 0 {width} {height}" '
            f'xmlns="http://www.w3.org/2000/svg">'
            f'<text x="{width//2}" y="{height//2}" text-anchor="middle" '
            f'font-size="14" fill="#aaa">エッジデータなし</text></svg>'
        )

    n = len(edge_flows)
    row_h   = min(40, (height - 60) // max(n, 1))
    bar_max = width - 200   # バーの最大幅

    lines = [
        f'<svg viewBox="0 0 {width} {height + 20}" '
        f'xmlns="http://www.w3.org/2000/svg" '
        f'style="background:#f8f9fb;border-radius:8px;'
        f'border:1px solid #e0e4ed;">',
        # タイトル
        f'<text x="10" y="22" font-size="12" font-weight="700" '
        f'fill="#1a3a5c">§4-3 エッジ障壁ヒートマップ</text>',
        f'<text x="10" y="36" font-size="9" fill="#888">'
        f'flow_ease = (1 - μ/μ_max) × (1 - T₅)　　'
        f'μ=2√T₂ [Derived]、T₅ [Derived]、μ_max [Assumed]</text>',
    ]

    # カラースケールの凡例
    legend_y = height + 5
    lines += [
        f'<text x="10" y="{legend_y}" font-size="9" fill="#888">低</text>',
        f'<rect x="30" y="{legend_y-10}" width="120" height="10" '
        f'fill="url(#grad_legend)"/>',
        f'<text x="155" y="{legend_y}" font-size="9" fill="#888">高</text>',
        '<defs>',
        '<linearGradient id="grad_legend" x1="0" y1="0" x2="1" y2="0">',
        '<stop offset="0%" stop-color="#e74c3c"/>',
        '<stop offset="50%" stop-color="#f8c471"/>',
        '<stop offset="100%" stop-color="#2ecc71"/>',
        '</linearGradient>',
        '</defs>',
    ]

    for i, ef in enumerate(edge_flows):
        y  = 50 + i * (row_h + 8)
        bw = int(ef.flow_ease * bar_max)

        # バーの色（flow_ease で補間）
        r, g = _ease_to_rgb(ef.flow_ease)
        color = f"rgb({r},{g},80)"

        # エッジ名
        lines.append(
            f'<text x="10" y="{y + row_h - 6}" '
            f'font-size="10" fill="#1a3a5c" font-weight="600">'
            f'{ef.src} → {ef.dst}</text>'
        )

        # バー背景
        lines.append(
            f'<rect x="140" y="{y}" width="{bar_max}" height="{row_h}" '
            f'rx="4" fill="#e0e4ed" opacity="0.5"/>'
        )
        # バー本体
        if bw > 0:
            lines.append(
                f'<rect x="140" y="{y}" width="{bw}" height="{row_h}" '
                f'rx="4" fill="{color}"/>'
            )
        # 数値ラベル
        lines.append(
            f'<text x="{140 + bw + 5}" y="{y + row_h - 6}" '
            f'font-size="10" fill="#555">'
            f'{ef.label}　T₂={ef.T2:.3f} T₅={ef.T5:.3f}</text>'
        )

    lines.append('</svg>')
    return "\n".join(lines)


def _ease_to_rgb(ease: float) -> Tuple[int, int]:
    """flow_ease → (R, G) の補間（赤→黄→緑）"""
    e = max(0.0, min(1.0, ease))
    if e < 0.5:
        r = 231
        g = int(196 * e * 2)
    else:
        r = int(231 * (1 - (e - 0.5) * 2))
        g = 204
    return r, g


def render_edge_heatmap_table_html(edge_flows: List[EdgeFlowData]) -> str:
    """
    エッジ障壁ヒートマップのテーブル HTML を生成する（§4-3 用）。
    """
    if not edge_flows:
        return '<p style="color:#aaa;font-size:.88rem;">エッジデータがありません。</p>'

    status_css = {"easy": "good", "moderate": "warn", "hard": "alert"}
    status_jp  = {"easy": "○ 流れやすい", "moderate": "△ 中程度", "hard": "× 流れにくい"}

    rows = []
    for ef in edge_flows:
        css = status_css.get(ef.status, "mid")
        txt = status_jp.get(ef.status, "—")
        bar_w = int(ef.flow_ease * 100)
        rows.append(
            f'<tr>'
            f'<td><strong>{ef.src} → {ef.dst}</strong></td>'
            f'<td>{ef.mu:.4f}</td>'
            f'<td>{ef.T5:.4f}</td>'
            f'<td>'
            f'<div style="display:flex;align-items:center;gap:8px;">'
            f'<div style="width:100px;background:#e0e4ed;border-radius:4px;height:8px;">'
            f'<div style="width:{bar_w}px;background:#1a3a5c;border-radius:4px;height:8px;"></div>'
            f'</div>'
            f'<span>{ef.flow_ease:.3f}</span>'
            f'</div>'
            f'</td>'
            f'<td><span class="kpi-status {css}">{txt}</span></td>'
            f'</tr>'
        )

    header = (
        '<table class="evidence-table" style="margin-top:10px;">'
        '<thead><tr>'
        '<th>エッジ</th>'
        '<th>μ=2√T₂（移動コスト）</th>'
        '<th>T₅（通信路損失）</th>'
        '<th>flow_ease</th>'
        '<th>状態</th>'
        '</tr></thead>'
        '<tbody>'
    )
    footer = (
        '</tbody></table>'
        '<p style="font-size:.8rem;color:#888;margin-top:8px;">'
        'flow_ease = (1 - μ/μ_max) × (1 - T₅)　'
        'μ=2√T₂ [Derived]　T₅ [Derived]　μ_max [Assumed]'
        '</p>'
    )
    return header + "\n".join(rows) + footer


# ============================================================
# ② Φ₀（介入閾値）の自動設定
# ============================================================

def auto_phi_threshold(
    state: StateSnapshot,
    method: str = "living_system_boundary",
) -> Tuple[float, str]:
    """
    グランドポテンシャルの介入閾値 Φ₀ を自動設定する。

    method = "living_system_boundary"（推奨）：
      Living System 相の境界条件（S≤0.6, D∈[1.2,1.7], H≥0.45）に
      対応する Φ の推定値を Φ₀ として設定する。

      現在の Φ = F - μN において、
      Living System 境界での S・D・H を使って F_boundary を推定し：
        Φ₀ = F_boundary - μ·N

    method = "ratio"：
      Φ₀ = Φ_current × ratio（デフォルト 0.95）
      簡易設定。理論的根拠は薄いが実用的。

    Returns
    -------
    (Phi_0, description) : Tuple[float, str]
    """
    if method == "living_system_boundary":
        # Living System 境界での目標エントロピー
        S_target   = 0.55 * K_B  # S ≤ 0.6 の中間点
        H_bit_tgt  = S_target / (K_B * LN2)

        # 現在の G から温度 T と体積 V は変わらないと仮定
        # F_boundary = U - T·S_boundary
        if state.G is not None:
            G  = state.G
            d  = G.shape[0]
            U  = float(np.trace(G) / d)
            try:
                G_inv = np.linalg.inv(G + np.eye(d) * EPS)
                T_val = float(np.trace(G_inv) / (d * K_B))
            except Exception:
                T_val = calc_T_temperature(G)
            F_boundary = U - T_val * S_target
        else:
            # G がない場合は現在の F をそのまま使用
            F_boundary = state.F

        Phi_0 = calc_grand_potential(F_boundary, state.mu, state.N)
        desc  = (
            f"Living System 境界からの自動設定 "
            f"（S_target={S_target:.2e} J/K、"
            f"F_boundary={F_boundary:.4f} J）"
        )

    elif method == "ratio":
        ratio = 0.95
        Phi_0 = state.Phi * ratio
        desc  = f"現在 Φ の {ratio*100:.0f}% を閾値として設定"

    elif method == "phase_transition_boundary":
        # 相転移予兆の条件（H < 0.45）に対応する Φ を推定
        # H が 0.45 を下回ったときの F を推定
        S_warn    = 0.65 * K_B   # Noise Dominant の境界
        if state.G is not None:
            G  = state.G
            d  = G.shape[0]
            U  = float(np.trace(G) / d)
            try:
                G_inv = np.linalg.inv(G + np.eye(d) * EPS)
                T_val = float(np.trace(G_inv) / (d * K_B))
            except Exception:
                T_val = calc_T_temperature(G)
            F_warn = U - T_val * S_warn
        else:
            F_warn = state.F
        Phi_0 = calc_grand_potential(F_warn, state.mu, state.N)
        desc  = "相転移予兆境界（Noise Dominant 手前）からの自動設定"

    else:
        Phi_0 = state.Phi * 0.95
        desc  = "デフォルト設定（現在 Φ の 95%）"

    return float(Phi_0), desc


def auto_tau(H: float) -> Tuple[float, str]:
    """
    予測ホライズン τ を Hurst 指数から自動設定する。

    H > 0.5（持続傾向）：トレンドが続くため長いホライズンが使える
    H ≈ 0.5（ランダム）：予測精度が低いため短いホライズンが適切
    H < 0.5（反持続）：反転しやすいため非常に短いホライズンが必要

    Returns
    -------
    (tau, description) : Tuple[float, str]
    """
    if H >= 0.6:
        tau  = 5.0
        desc = f"H={H:.3f}（高持続）→ τ=5 ステップ"
    elif H >= 0.45:
        tau  = 3.0
        desc = f"H={H:.3f}（中程度）→ τ=3 ステップ"
    elif H >= 0.3:
        tau  = 2.0
        desc = f"H={H:.3f}（低持続）→ τ=2 ステップ"
    else:
        tau  = 1.0
        desc = f"H={H:.3f}（反持続傾向）→ τ=1 ステップ（最短）"
    return tau, desc


# ============================================================
# ③ レポート生成パイプライン
# ============================================================

@dataclass
class ReportContext:
    """
    レポート生成に必要な全データをまとめたコンテキスト。
    Python からプレースホルダーを埋めて HTML を生成する。
    """
    # 基本情報
    title:            str = "SDFT 構造場分析レポート"
    subject:          str = "分析対象"
    purpose:          str = "構造場診断"

    # 構造場変数
    H_bit: float = 0.0
    S:     float = 0.0
    D:     float = 0.0
    H:     float = 0.5
    F:     float = 0.0
    T:     float = 0.0
    U:     float = 0.0
    V:     float = 0.0
    Phi:   float = 0.0
    N:     float = 0.0
    mu:    float = 0.0
    dPhi_dt: float = 0.0
    Phi_threshold: float = 0.0
    Phi_predicted: float = 0.0

    # 𝓣 ベクトル
    T1: float = 0.0; T1_evidence: str = "Derived"
    T2: float = 0.0; T2_evidence: str = "Derived"
    T3: float = 0.0; T3_evidence: str = "Derived"; T3_alpha: float = 1.0
    T4: float = 0.0; T4_evidence: str = "Assumed"
    T5: float = 0.0; T5_evidence: str = "Derived"

    # 位相判定
    phase_label:       str = "Living System"
    phase_css:         str = "phase-living"
    phase_explanation: str = ""
    phase_transition_alert: str = ""

    # SVG 座標
    sd_x:             int = 160
    sd_y:             int = 130
    sd_y_label:       int = 116
    h_polyline:       str = "40,130 280,130"
    h_current_x:      int = 280
    h_current_y:      int = 130
    h_current_y_label: int = 116

    # §1 エグゼクティブサマリー
    exec_situation: str = ""
    exec_problem:   str = ""
    exec_action_1:  str = ""
    exec_action_2:  str = ""
    exec_action_3:  str = ""
    exec_impact:    str = ""
    exec_top_risk:  str = ""

    # §3 介入設計
    intervention_needed: str = "不要"
    delta_phi_need:      float = 0.0
    sensitivity_N:       float = 0.0
    sensitivity_T2:      float = 0.0
    sensitivity_theta:   float = 0.0
    priority_N:          float = 0.0
    priority_T2:         float = 0.0
    priority_theta:      float = 0.0
    best_variable:       str = "N"
    delta_u_optimal:     float = 0.0
    cost_N:    float = 1.0
    cost_T2:   float = 1.0
    cost_theta: float = 10.0
    lambda_val: float = 1.0
    total_interventions: int = 0

    # §4-2 行動条件診断
    action_condition: Optional[ActionConditionDiagnosis] = None

    # §4-3 エッジ障壁ヒートマップ
    edge_flows: List[EdgeFlowData] = field(default_factory=list)

    # §3 Φ エネルギーグラフ
    phi_graph_data: Optional[PhiGraphData] = None

    # §4-4 証拠分類
    evidence_rows: str = ""

    # §5 シナリオ
    scenario_rows: str = ""

    # §6 総論
    uncertainty_desc:  str = ""
    additional_obs_rows: str = ""

    # §4-1 SDFT 写像（外部から設定。省略時はデフォルト値を使用）
    mapping_nodes:    str = ""
    mapping_links:    str = ""
    mapping_field:    str = ""
    mapping_boundary: str = ""
    mapping_regime:   str = ""

    @classmethod
    def from_state(
        cls,
        state: StateSnapshot,
        phi_prev: float = 0.0,
        dt: float = 1.0,
        regime_data: Optional[Dict[str, List[float]]] = None,
        edges: Optional[List[Tuple[str, str]]] = None,
        intervention_result=None,
        opt_lambda: float = 1.0,
        H_series: Optional[List[float]] = None,
        title: str = "SDFT 構造場分析レポート",
        subject: str = "分析対象",
        purpose: str = "構造場診断",
    ) -> "ReportContext":
        """StateSnapshot から ReportContext を生成する。"""
        ctx = cls()
        ctx.title   = title
        ctx.subject = subject
        ctx.purpose = purpose

        # 構造場変数
        ctx.H_bit = state.H_bit
        ctx.S     = state.S
        ctx.D     = state.D
        ctx.H     = state.H
        ctx.F     = state.F
        ctx.T     = state.T
        ctx.U     = state.U
        ctx.V     = state.V
        ctx.Phi   = state.Phi
        ctx.N     = state.N
        ctx.mu    = state.mu
        ctx.dPhi_dt = (state.Phi - phi_prev) / (dt + EPS)

        # Φ₀ の自動設定
        ctx.Phi_threshold, _ = auto_phi_threshold(state)
        ctx.Phi_predicted     = state.Phi + ctx.dPhi_dt * auto_tau(state.H)[0]

        # 𝓣
        if state.tension:
            t = state.tension
            ctx.T1 = t.T1; ctx.T1_evidence = t.T1_evidence
            ctx.T2 = t.T2; ctx.T2_evidence = t.T2_evidence
            ctx.T3 = t.T3; ctx.T3_evidence = t.T3_evidence; ctx.T3_alpha = t.T3_alpha
            ctx.T4 = t.T4; ctx.T4_evidence = t.T4_evidence
            ctx.T5 = t.T5; ctx.T5_evidence = t.T5_evidence

        # 位相判定
        S_norm = state.H_bit / math.log2(32 + EPS)
        lbl, css, expl = classify_phase(S_norm, state.D, state.H)
        ctx.phase_label       = lbl
        ctx.phase_css         = css
        ctx.phase_explanation = expl

        # SVG 座標
        ctx.sd_x, ctx.sd_y, ctx.sd_y_label = calc_sd_position(S_norm, state.D)
        h_series = H_series or [state.H]
        (ctx.h_polyline, ctx.h_current_x,
         ctx.h_current_y, ctx.h_current_y_label) = calc_h_polyline(h_series)

        # §3 介入設計
        if intervention_result and hasattr(intervention_result, 'steps'):
            steps = intervention_result.steps
            if steps:
                last = steps[-1]
                ctx.intervention_needed = "必要" if last.needs_intervention else "不要"
                ctx.delta_phi_need      = last.delta_phi_need
                ctx.sensitivity_N       = last.gradients.get('N', 0.0)
                ctx.sensitivity_T2      = last.gradients.get('T2', 0.0)
                ctx.sensitivity_theta   = last.gradients.get('theta', 0.0)
                ctx.priority_N          = last.priorities.get('N', 0.0)
                ctx.priority_T2         = last.priorities.get('T2', 0.0)
                ctx.priority_theta      = last.priorities.get('theta', 0.0)
                ctx.best_variable       = last.best_variable
                ctx.delta_u_optimal     = last.delta_u_optimal
            ctx.total_interventions = intervention_result.total_interventions
        ctx.lambda_val = opt_lambda

        # §4-2 行動条件診断
        ctx.action_condition = state.action_condition

        # §4-3 エッジ障壁ヒートマップ
        if regime_data and edges:
            ctx.edge_flows = calc_edge_flows(regime_data, edges)

        # §3 Φ エネルギーグラフ
        if intervention_result is not None:
            ctx.phi_graph_data = build_phi_graph_data(
                intervention_result,
                phi_threshold=float(ctx.Phi_threshold),
            )

        # §4-4 証拠分類表（既定の全変数）
        ctx.evidence_rows = _build_evidence_rows(ctx)

        return ctx


def _build_evidence_rows(ctx: "ReportContext") -> str:
    """証拠分類テーブルの HTML 行を生成する。"""
    items = [
        ("H_bit（シャノンエントロピー）", f"{ctx.H_bit:.4f} bit",
         "Derived", "時系列ヒストグラムから計算"),
        ("S（熱力学的エントロピー）", f"{ctx.S:.4e} J/K",
         "Derived", "S = k_B·ln2·H_bit"),
        ("D（フラクタル次元）", f"{ctx.D:.4f}",
         "Derived", "ヒグチ法"),
        ("H（Hurst指数）", f"{ctx.H:.4f}",
         "Derived", "ラグ別分散法"),
        ("T（情報的温度）", f"{ctx.T:.4e} K",
         "Derived", "tr(G⁻¹)/(d·k_B)"),
        ("U（内部エネルギー）", f"{ctx.U:.4f} J相当",
         "Derived", "tr(G)/d"),
        ("F（自由エネルギー）", f"{ctx.F:.4f} J",
         "Derived", "U - T·S"),
        ("N（人数）", f"{ctx.N:.0f} 人",
         "Observed", "直接計測"),
        ("μ（化学ポテンシャル）", f"{ctx.mu:.4f} J/人",
         "Derived", "2√T₂"),
        ("Φ（グランドポテンシャル）", f"{ctx.Phi:.4f} J",
         "Derived", "F - μN"),
        ("T₁（対数分配関数差）", f"{ctx.T1:.4f}",
         ctx.T1_evidence, "|ψ_A - ψ_B|"),
        ("T₂（Bhattacharyya距離）", f"{ctx.T2:.4f}",
         ctx.T2_evidence, "-log BC(p_A, p_B)"),
        (f"T₃（α-ダイバージェンス α={ctx.T3_alpha}）", f"{ctx.T3:.4f}",
         ctx.T3_evidence, f"D^(α={ctx.T3_alpha})(p_A‖p_B)"),
        ("T₄（統計的ホロノミー）", f"{ctx.T4:.4f}",
         ctx.T4_evidence, "‖Hol^(α)(γ) - I‖（近似）"),
        ("T₅（通信路損失率）", f"{ctx.T5:.4f}",
         ctx.T5_evidence, "1 - I(X_A;X_B)/H(X_A)"),
    ]
    ev_css = {
        "Observed": "ev-observed",
        "Derived":  "ev-derived-t",
        "Assumed":  "ev-assumed-t",
        "Hypothesis": "ev-hypothesis",
    }
    rows = []
    for name, val, ev, note in items:
        css = ev_css.get(ev, "")
        rows.append(
            f'<tr><td>{name}</td><td>{val}</td>'
            f'<td class="{css}">{ev}</td><td>{note}</td></tr>'
        )
    return "\n".join(rows)


def _build_action_condition_rows(
    diag: Optional[ActionConditionDiagnosis],
) -> str:
    """行動条件診断テーブルの HTML 行を生成する。"""
    if diag is None:
        return '<tr><td colspan="4" style="color:#aaa;">診断データなし</td></tr>'
    status_mark = {"good": "○", "warn": "△", "alert": "×", "mid": "－"}
    status_css  = {"good": "good", "warn": "warn", "alert": "alert", "mid": "mid"}
    rows = []
    for row in diag.table_rows():
        mark = status_mark.get(row["status"], "－")
        css  = status_css.get(row["status"], "mid")
        rows.append(
            f'<tr>'
            f'<td><span class="kpi-status {css}">{mark}</span> {row["item"]}</td>'
            f'<td style="font-size:.85rem;color:#555;">{row["variable"]}</td>'
            f'<td class="{css}" style="font-size:.85rem;">{row["status"].upper()}</td>'
            f'<td style="font-size:.85rem;">{row["interpretation"]}</td>'
            f'</tr>'
        )
    return "\n".join(rows)


def _priority_bar_pct(priorities: Dict[str, float], key: str) -> int:
    max_p = max(priorities.values()) if priorities else 1.0
    return int(priorities.get(key, 0.0) / (max_p + EPS) * 100)


def generate_report(
    ctx:           "ReportContext",
    template_path: str,
    output_path:   str = "",
) -> str:
    """
    テンプレートにプレースホルダーを埋め込んで HTML レポートを生成する。

    Parameters
    ----------
    ctx           : ReportContext（全データ）
    template_path : report_template_v006.html のパス
    output_path   : 保存先パス（省略時は返値のみ）

    Returns
    -------
    HTML 文字列
    """
    with open(template_path, 'r', encoding='utf-8') as f:
        html = f.read()

    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")

    # KPI ステータス
    from sdft_v006_core import get_kpi_status
    S_norm = ctx.H_bit / math.log2(32 + EPS)
    s_css, s_txt = get_kpi_status('S', S_norm)
    d_css, d_txt = get_kpi_status('D', ctx.D)
    h_css, h_txt = get_kpi_status('H', ctx.H)

    # F の状態（正負で判定）
    f_css  = "good" if ctx.F > 0 else "alert"
    f_txt  = "正 ✓" if ctx.F > 0 else "負 ⚠"

    # priority パーセント
    prios = {
        'N':     ctx.priority_N,
        'T2':    ctx.priority_T2,
        'theta': ctx.priority_theta,
    }
    max_prio = max(prios.values()) if max(prios.values()) > 0 else 1.0
    p_n_pct  = int(ctx.priority_N    / (max_prio + EPS) * 100)
    p_t2_pct = int(ctx.priority_T2   / (max_prio + EPS) * 100)
    p_th_pct = int(ctx.priority_theta / (max_prio + EPS) * 100)

    # best variable badge
    best_n_badge  = '<span class="best-badge">★ 推奨</span>' if ctx.best_variable == 'N'     else ''
    best_t2_badge = '<span class="best-badge">★ 推奨</span>' if ctx.best_variable == 'T2'    else ''
    best_th_badge = '<span class="best-badge">★ 推奨</span>' if ctx.best_variable == 'theta' else ''

    # 相転移アラート
    alert_html = ""
    if ctx.phase_transition_alert:
        alert_html = (
            f'<div class="phase-alert-box">'
            f'{ctx.phase_transition_alert}'
            f'</div>'
        )

    # エッジヒートマップ HTML
    edge_hm_html = render_edge_heatmap_table_html(ctx.edge_flows)
    edge_svg      = render_edge_heatmap_svg(ctx.edge_flows)

    # 行動条件診断 HTML
    acd_rows = _build_action_condition_rows(ctx.action_condition)
    acd_overall_css = ""
    acd_overall_status = "—"
    acd_overall_msg    = ""
    if ctx.action_condition:
        acd_overall_css    = ctx.action_condition.overall_status
        acd_overall_status = {
            "good": "○ 行動が起きやすい",
            "warn": "△ 一部に障壁あり",
            "alert": "× 行動が起きにくい",
            "mid": "― 中程度",
        }.get(ctx.action_condition.overall_status, "—")
        acd_overall_msg = ctx.action_condition.overall_message

    # エグゼクティブサマリーのデフォルト
    exec_situation = ctx.exec_situation or f"現在のフェーズは {ctx.phase_label}。"
    exec_problem   = ctx.exec_problem   or "詳細な分析を要します。"
    exec_action_1  = ctx.exec_action_1  or "データの継続的な観測を実施する。"
    exec_action_2  = ctx.exec_action_2  or "H（Hurst指数）の推移を監視する。"
    exec_action_3  = ctx.exec_action_3  or f"最優先介入変数（{ctx.best_variable}）への施策を検討する。"
    exec_impact    = ctx.exec_impact    or "Φ の維持・向上が期待される。"
    exec_top_risk  = ctx.exec_top_risk  or "相転移への移行リスク。"

    replacements = {
        # ヘッダー
        "{{REPORT_TITLE}}":       ctx.title,
        "{{ANALYSIS_DATE}}":      now,
        "{{SUBJECT}}":            ctx.subject,
        "{{PHASE_LABEL}}":        ctx.phase_label,
        "{{ANALYSIS_PURPOSE}}":   ctx.purpose,
        # §1
        "{{EXEC_SITUATION}}":     exec_situation,
        "{{EXEC_PROBLEM}}":       exec_problem,
        "{{EXEC_ACTION_1}}":      exec_action_1,
        "{{EXEC_ACTION_2}}":      exec_action_2,
        "{{EXEC_ACTION_3}}":      exec_action_3,
        "{{EXEC_IMPACT}}":        exec_impact,
        "{{EXEC_TOP_RISK}}":      exec_top_risk,
        # §2 KPI
        "{{PHASE_CSS}}":          ctx.phase_css,
        "{{PHASE_EXPLANATION}}":  ctx.phase_explanation,
        "{{S_VALUE}}":            f"{S_norm:.4f}",
        "{{S_STATUS_CSS}}":       s_css,
        "{{S_STATUS}}":           s_txt,
        "{{D_VALUE}}":            f"{ctx.D:.4f}",
        "{{D_STATUS_CSS}}":       d_css,
        "{{D_STATUS}}":           d_txt,
        "{{H_VALUE}}":            f"{ctx.H:.4f}",
        "{{H_STATUS_CSS}}":       h_css,
        "{{H_STATUS}}":           h_txt,
        "{{F_VALUE}}":            f"{ctx.F:.4f}",
        "{{F_STATUS_CSS}}":       f_css,
        "{{F_STATUS}}":           f_txt,
        "{{MU_VALUE}}":           f"{ctx.mu:.4f}",
        "{{PHI_VALUE}}":          f"{ctx.Phi:.4f}",
        "{{PHI_DELTA}}":          f"{ctx.dPhi_dt:.4f}",
        "{{PHI_THRESHOLD}}":      f"{ctx.Phi_threshold:.4f}",
        "{{PHI_PREDICTED}}":      f"{ctx.Phi_predicted:.4f}",
        "{{N_VALUE}}":            f"{ctx.N:.0f}",
        "{{T1_VALUE}}":           f"{ctx.T1:.4f}",
        "{{T1_EVIDENCE}}":        ctx.T1_evidence,
        "{{T2_VALUE}}":           f"{ctx.T2:.4f}",
        "{{T2_EVIDENCE}}":        ctx.T2_evidence,
        "{{T3_VALUE}}":           f"{ctx.T3:.4f}",
        "{{T3_EVIDENCE}}":        ctx.T3_evidence,
        "{{T3_ALPHA}}":           str(ctx.T3_alpha),
        "{{T4_VALUE}}":           f"{ctx.T4:.4f}",
        "{{T4_EVIDENCE}}":        ctx.T4_evidence,
        "{{T5_VALUE}}":           f"{ctx.T5:.4f}",
        "{{T5_EVIDENCE}}":        ctx.T5_evidence,
        "{{PHASE_TRANSITION_ALERT}}": alert_html,
        # §3 介入
        "{{INTERVENTION_NEEDED}}":  ctx.intervention_needed,
        "{{DELTA_PHI_NEED}}":       f"{ctx.delta_phi_need:.4f}",
        "{{SENSITIVITY_N}}":        f"{ctx.sensitivity_N:.4f}",
        "{{SENSITIVITY_T2}}":       f"{ctx.sensitivity_T2:.4f}",
        "{{SENSITIVITY_THETA}}":    f"{ctx.sensitivity_theta:.4e}",
        "{{PRIORITY_N}}":           f"{ctx.priority_N:.4f}",
        "{{PRIORITY_T2}}":          f"{ctx.priority_T2:.4f}",
        "{{PRIORITY_THETA}}":       f"{ctx.priority_theta:.4f}",
        "{{PRIORITY_N_PCT}}":       str(p_n_pct),
        "{{PRIORITY_T2_PCT}}":      str(p_t2_pct),
        "{{PRIORITY_THETA_PCT}}":   str(p_th_pct),
        "{{BEST_VARIABLE}}":        ctx.best_variable,
        "{{DELTA_U_OPTIMAL}}":      f"{ctx.delta_u_optimal:.4f}",
        "{{BEST_N_BADGE}}":         best_n_badge,
        "{{BEST_T2_BADGE}}":        best_t2_badge,
        "{{BEST_THETA_BADGE}}":     best_th_badge,
        "{{COST_N}}":               f"{ctx.cost_N:.1f}",
        "{{COST_T2}}":              f"{ctx.cost_T2:.1f}",
        "{{COST_THETA}}":           f"{ctx.cost_theta:.1f}",
        # §4-2 行動条件診断
        "{{ACTION_CONDITION_ROWS}}": acd_rows,
        "{{ACTION_OVERALL_CSS}}":    acd_overall_css,
        "{{ACTION_OVERALL_STATUS}}": acd_overall_status,
        "{{ACTION_OVERALL_MESSAGE}}": acd_overall_msg,
        # §4-3 エッジ障壁ヒートマップ
        "{{EDGE_HEATMAP_TABLE}}":   edge_hm_html,
        "{{EDGE_HEATMAP_SVG}}":     edge_svg,
        # §3 Φ エネルギーグラフ
        "{{PHI_ENERGY_GRAPH_SVG}}": render_phi_energy_graph_svg(ctx.phi_graph_data)
            if ctx.phi_graph_data else '<p style="color:#aaa;">データなし</p>',
        "{{PHI_INTERVENTION_TABLE}}": render_phi_intervention_table_html(ctx.phi_graph_data)
            if ctx.phi_graph_data else '<p style="color:#aaa;">データなし</p>',
        # §4-4 証拠分類
        "{{EVIDENCE_ROWS}}":        ctx.evidence_rows,
        # SVG 位相図
        "{{SD_X}}":                 str(ctx.sd_x),
        "{{SD_Y}}":                 str(ctx.sd_y),
        "{{SD_Y_LABEL}}":           str(ctx.sd_y_label),
        "{{H_POLYLINE}}":           ctx.h_polyline,
        "{{H_CURRENT_X}}":          str(ctx.h_current_x),
        "{{H_CURRENT_Y}}":          str(ctx.h_current_y),
        "{{H_CURRENT_Y_LABEL}}":    str(ctx.h_current_y_label),
        # §5 シナリオ
        "{{SCENARIO_ROWS}}":        ctx.scenario_rows or _default_scenario_rows(ctx),
        # §6 総論
        "{{UNCERTAINTY_DESC}}":     ctx.uncertainty_desc or "追加データによる精度向上が望ましい。",
        "{{ADDITIONAL_OBS_ROWS}}":  ctx.additional_obs_rows or _default_obs_rows(),
        # §4-1 SDFT 写像
        "{{MAPPING_NODES}}":    ctx.mapping_nodes    or "分析対象の構成要素（部門・ゾーン・セグメント）",
        "{{MAPPING_LINKS}}":    ctx.mapping_links    or "構成要素間の接続（強度・方向のみ）",
        "{{MAPPING_FIELD}}":    ctx.mapping_field    or "外部環境（市場・競合・マクロ環境）",
        "{{MAPPING_BOUNDARY}}": ctx.mapping_boundary or "レジームの境界（分析対象の範囲）",
        "{{MAPPING_REGIME}}":   ctx.mapping_regime   or "同一オペモードを持つ領域（A / B）",
    }

    for placeholder, value in replacements.items():
        html = html.replace(placeholder, str(value))

    if output_path:
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(html)
        print(f"レポートを保存しました: {output_path}")

    return html


def _default_scenario_rows(ctx: "ReportContext") -> str:
    """デフォルトのシナリオ表 HTML 行を生成する。"""
    rows = [
        ("H（Hurst指数）",
         f"{min(ctx.H + 0.1, 1.0):.2f}（上昇）",
         f"{ctx.H:.2f}（維持）",
         f"{max(ctx.H - 0.1, 0.0):.2f}（低下）"),
        ("S（エントロピー）",
         "低下（秩序回復）", "現状維持", "上昇（散逸進行）"),
        ("Φ（グランドポテンシャル）",
         "維持 or 上昇", f"{ctx.Phi:.2f}（現在）", "減少継続"),
    ]
    html_rows = []
    for label, opt, base, pess in rows:
        html_rows.append(
            f'<tr><td>{label}</td>'
            f'<td style="color:#1e8449;">{opt}</td>'
            f'<td>{base}</td>'
            f'<td style="color:#c0392b;">{pess}</td></tr>'
        )
    return "\n".join(html_rows)


def _default_obs_rows() -> str:
    return (
        '<tr><td>𝓣 の時系列</td><td>レジーム間テンションの推移把握</td><td>高</td></tr>'
        '<tr><td>N（人数）の実測</td><td>化学ポテンシャル精度向上</td><td>高</td></tr>'
        '<tr><td>両レジームの同時観測</td><td>T₅ の精度向上（KSG推定量）</td><td>中</td></tr>'
    )


# ============================================================
# メイン実行例
# ============================================================

if __name__ == "__main__":
    import os, random
    random.seed(42)

    print("=" * 60)
    print("SDFT revised v0.0.6 — レポート生成パイプライン 動作確認")
    print("=" * 60)

    # サンプルデータ
    n = 200
    regime_A = [500 + random.gauss(0, 60) + i * 0.2 for i in range(n)]
    regime_B = [350 + random.gauss(0, 80) for _ in range(n)]
    regime_C = [420 + random.gauss(0, 70) for _ in range(n)]

    # StateSnapshot を計算
    snap = StateSnapshot(t=0.0, x=regime_A, x_B=regime_B, N=100.0)
    snap.compute(phi_prev=0.0, dt=1.0, alpha_T3=1.0)

    print(f"\nΦ = {snap.Phi:.4f}")
    print(f"行動条件: {snap.action_condition.overall_status.upper()} "
          f"{snap.action_condition.overall_message}")

    # Φ₀ の自動設定
    phi_0, phi_0_desc = auto_phi_threshold(snap)
    print(f"\nΦ₀（自動設定）= {phi_0:.4f}")
    print(f"設定根拠: {phi_0_desc}")

    # τ の自動設定
    tau, tau_desc = auto_tau(snap.H)
    print(f"τ（自動設定）= {tau}")
    print(f"設定根拠: {tau_desc}")

    # エッジ障壁ヒートマップ
    regime_data = {"A": regime_A, "B": regime_B, "C": regime_C}
    edges       = [("A", "B"), ("B", "A"), ("A", "C"), ("C", "A")]
    edge_flows  = calc_edge_flows(regime_data, edges)

    print("\n--- エッジ障壁ヒートマップ ---")
    for ef in edge_flows:
        print(f"  {ef.src} → {ef.dst}: flow_ease={ef.flow_ease:.3f}  {ef.label}")

    # 介入ループ
    x_series   = [regime_A] * 10
    x_B_series = [regime_B] * 10
    N_series   = [100.0] * 10
    states     = generate_state_sequence(x_series, x_B_series, N_series)
    tradeoff   = scan_lambda_tradeoff(states, phi_0, lambda_values=[0.1, 1.0, 5.0])
    opt_lam, _ = find_optimal_lambda(tradeoff)
    result     = run_intervention_loop(states, phi_0, tau=tau, lambda_count=opt_lam)

    # ReportContext の生成
    ctx = ReportContext.from_state(
        state=snap,
        phi_prev=0.0,
        regime_data=regime_data,
        edges=edges,
        intervention_result=result,
        opt_lambda=opt_lam,
        H_series=[s.H for s in states],
        title="SDFT v0.0.6 動作確認レポート",
        subject="サンプル施設",
        purpose="構造場診断・介入最適化",
    )

    # レポート生成（テンプレートがある場合）
    template_path = os.path.join(
        os.path.dirname(__file__), '..', 'assets', 'report_template_v006.html'
    )
    output_path = "/mnt/user-data/outputs/sdft-report-v006-test.html"

    if os.path.exists(template_path):
        generate_report(ctx, template_path, output_path)
        print(f"\nレポート生成成功: {output_path}")
    else:
        print(f"\nテンプレートが見つかりません: {template_path}")
        print("ReportContext の生成は成功しています。")

    print("\n" + "=" * 60)
    print("動作確認完了")
    print("=" * 60)

# ============================================================
# §3 Φ エネルギー-時刻グラフ（追記）
# ============================================================

