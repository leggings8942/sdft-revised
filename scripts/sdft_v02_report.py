"""
sdft_v02_report.py
SDFT revised v0.2 — レポート生成パイプライン

v0.0.6 からの変更点：
  ・H_bit/S/H/T/U/V/G → S_g/D/λ₁ (Layer 1)
  ・F → F_geom (λ₂, スペクトルギャップ)
  ・μ = 2√T₂ → μ = T₂（測地距離）
  ・T₃ α-ダイバージェンス → PCA フレーム回転角（簡易版）
  ・T₄ ホロノミー → 曲率差（簡易版）
  ・EdgeFlowData: flow_ease = (1 - μ/μ_max) × (1 - T₅) の μ 定義更新
  ・位相空間図: S-D 軸 → S_g-D 軸、H 時系列 → λ₁ 時系列
  ・証拠分類表を多様体変数に更新
"""

from __future__ import annotations

import math
import datetime
from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Optional
import numpy as np

from sdft_v02_embedding import Evidence, EPS
from sdft_v02_potential import (
    calc_mu, calc_grand_potential, classify_phase,
    detect_phase_transition, get_kpi_status,
    calc_sd_position, calc_lambda1_polyline,
    _normalize_S_g, _normalize_D,
)
from sdft_v02_intervention import (
    ActionConditionDiagnosis,
    calc_action_condition,
    StateSnapshot,
    generate_state_sequence,
    run_intervention_loop,
    find_optimal_lambda,
    scan_lambda_tradeoff,
    auto_phi_threshold,
    auto_tau,
)
from sdft_v02_phi_graph import (
    PhiGraphData,
    build_phi_graph_data,
    render_phi_energy_graph_svg,
    render_phi_intervention_table_html,
)


# ============================================================
# エッジ障壁ヒートマップ
# ============================================================

@dataclass
class EdgeFlowData:
    """エッジの流れやすさデータ（多様体版）。"""
    src:       str
    dst:       str
    T2:        float   # 測地距離
    T5:        float   # 埋め込み歪み
    mu:        float   # = T₂（化学ポテンシャル）
    flow_ease: float   # (1 - μ/μ_max) × (1 - T₅)
    status:    str
    label:     str

    @classmethod
    def compute(
        cls,
        src: str, dst: str,
        embedded_src: np.ndarray,
        embedded_dst: np.ndarray,
        mu_max: float = 5.0,
        k: int = 10,
    ) -> "EdgeFlowData":
        from sdft_v02_geometry import calc_geodesic_distance
        from sdft_v02_tension import calc_T5

        T2_val = calc_geodesic_distance(embedded_src, embedded_dst, k=k)
        T5_val, _ = calc_T5(embedded_src, embedded_dst, k=min(k, 5))
        mu_val = T2_val  # μ = T₂

        a1 = max(0.0, 1.0 - mu_val / (mu_max + EPS))
        a2 = max(0.0, 1.0 - T5_val)
        ease = float(a1 * a2)

        if ease >= 0.6:
            status, label = "easy", f"流れやすい ({ease:.2f})"
        elif ease >= 0.3:
            status, label = "moderate", f"中程度 ({ease:.2f})"
        else:
            status, label = "hard", f"流れにくい ({ease:.2f})"

        return cls(src=src, dst=dst, T2=T2_val, T5=T5_val,
                   mu=mu_val, flow_ease=ease, status=status, label=label)


def calc_edge_flows(
    regime_embeddings: Dict[str, np.ndarray],
    edges: List[Tuple[str, str]],
    mu_max: float = 5.0,
    k: int = 10,
) -> List[EdgeFlowData]:
    results = []
    for src, dst in edges:
        if src in regime_embeddings and dst in regime_embeddings:
            ef = EdgeFlowData.compute(
                src=src, dst=dst,
                embedded_src=regime_embeddings[src],
                embedded_dst=regime_embeddings[dst],
                mu_max=mu_max, k=k,
            )
            results.append(ef)
    return results


def _ease_to_rgb(ease: float) -> Tuple[int, int]:
    e = max(0.0, min(1.0, ease))
    if e < 0.5:
        r, g = 231, int(196 * e * 2)
    else:
        r, g = int(231 * (1 - (e - 0.5) * 2)), 204
    return r, g


def render_edge_heatmap_svg(edge_flows: List[EdgeFlowData], width=600, height=300) -> str:
    if not edge_flows:
        return (
            f'<svg viewBox="0 0 {width} {height}" xmlns="http://www.w3.org/2000/svg">'
            f'<text x="{width//2}" y="{height//2}" text-anchor="middle" '
            f'font-size="14" fill="#aaa">エッジデータなし</text></svg>'
        )

    n = len(edge_flows)
    row_h = min(40, (height - 60) // max(n, 1))
    bar_max = width - 200

    lines = [
        f'<svg viewBox="0 0 {width} {height + 20}" xmlns="http://www.w3.org/2000/svg" '
        f'style="background:#f8f9fb;border-radius:8px;border:1px solid #e0e4ed;">',
        f'<text x="10" y="22" font-size="12" font-weight="700" fill="#1a3a5c">'
        f'§4-3 エッジ障壁ヒートマップ</text>',
        f'<text x="10" y="36" font-size="9" fill="#888">'
        f'flow_ease = (1 - μ/μ_max) × (1 - T₅)　μ=T₂（測地距離）</text>',
    ]

    legend_y = height + 5
    lines += [
        f'<text x="10" y="{legend_y}" font-size="9" fill="#888">低</text>',
        f'<rect x="30" y="{legend_y-10}" width="120" height="10" fill="url(#grad_legend)"/>',
        f'<text x="155" y="{legend_y}" font-size="9" fill="#888">高</text>',
        '<defs><linearGradient id="grad_legend" x1="0" y1="0" x2="1" y2="0">',
        '<stop offset="0%" stop-color="#e74c3c"/>',
        '<stop offset="50%" stop-color="#f8c471"/>',
        '<stop offset="100%" stop-color="#2ecc71"/>',
        '</linearGradient></defs>',
    ]

    for i, ef in enumerate(edge_flows):
        y = 50 + i * (row_h + 8)
        bw = int(ef.flow_ease * bar_max)
        r, g = _ease_to_rgb(ef.flow_ease)
        color = f"rgb({r},{g},80)"
        lines += [
            f'<text x="10" y="{y + row_h - 6}" font-size="10" fill="#1a3a5c" '
            f'font-weight="600">{ef.src} → {ef.dst}</text>',
            f'<rect x="140" y="{y}" width="{bar_max}" height="{row_h}" '
            f'rx="4" fill="#e0e4ed" opacity="0.5"/>',
        ]
        if bw > 0:
            lines.append(
                f'<rect x="140" y="{y}" width="{bw}" height="{row_h}" rx="4" fill="{color}"/>'
            )
        lines.append(
            f'<text x="{140 + bw + 5}" y="{y + row_h - 6}" font-size="10" fill="#555">'
            f'{ef.label}　T₂={ef.T2:.3f} T₅={ef.T5:.3f}</text>'
        )

    lines.append('</svg>')
    return "\n".join(lines)


def render_edge_heatmap_table_html(edge_flows: List[EdgeFlowData]) -> str:
    if not edge_flows:
        return '<p style="color:#aaa;font-size:.88rem;">エッジデータがありません。</p>'

    status_css = {"easy": "good", "moderate": "warn", "hard": "alert"}
    status_jp = {"easy": "○ 流れやすい", "moderate": "△ 中程度", "hard": "× 流れにくい"}

    rows = []
    for ef in edge_flows:
        css = status_css.get(ef.status, "mid")
        txt = status_jp.get(ef.status, "—")
        bar_w = int(ef.flow_ease * 100)
        rows.append(
            f'<tr><td><strong>{ef.src} → {ef.dst}</strong></td>'
            f'<td>{ef.mu:.4f}</td><td>{ef.T5:.4f}</td>'
            f'<td><div style="display:flex;align-items:center;gap:8px;">'
            f'<div style="width:100px;background:#e0e4ed;border-radius:4px;height:8px;">'
            f'<div style="width:{bar_w}px;background:#1a3a5c;border-radius:4px;height:8px;"></div>'
            f'</div><span>{ef.flow_ease:.3f}</span></div></td>'
            f'<td><span class="kpi-status {css}">{txt}</span></td></tr>'
        )

    header = (
        '<table class="evidence-table" style="margin-top:10px;">'
        '<thead><tr><th>エッジ</th><th>μ=T₂（測地距離）</th>'
        '<th>T₅（埋め込み歪み）</th><th>flow_ease</th><th>状態</th>'
        '</tr></thead><tbody>'
    )
    footer = (
        '</tbody></table>'
        '<p style="font-size:.8rem;color:#888;margin-top:8px;">'
        'flow_ease = (1 - μ/μ_max) × (1 - T₅)　μ=T₂ [Derived]　T₅ [Derived]　μ_max [Assumed]</p>'
    )
    return header + "\n".join(rows) + footer


# ============================================================
# ReportContext
# ============================================================

@dataclass
class ReportContext:
    title:   str = "SDFT v0.2 構造場分析レポート"
    subject: str = "分析対象"
    purpose: str = "構造場診断"

    # Layer 1
    S_g:      float = 0.0
    D:        float = 0.0
    lambda_1: float = 0.0
    d_embed:  int   = 3

    # Layer 2
    F_geom:    float = 0.0
    curvature: float = 0.0
    Phi:       float = 0.0
    N:         float = 0.0
    mu:        float = 0.0
    dPhi_dt:   float = 0.0
    Phi_threshold: float = 0.0
    Phi_predicted: float = 0.0

    # 𝓣
    T1: float = 0.0; T1_evidence: str = "Derived"
    T2: float = 0.0; T2_evidence: str = "Derived"
    T3: float = 0.0; T3_evidence: str = "Assumed"
    T4: float = 0.0; T4_evidence: str = "Assumed"
    T5: float = 0.0; T5_evidence: str = "Derived"

    # Phase
    phase_label:       str = "Living System"
    phase_css:         str = "phase-living"
    phase_explanation: str = ""
    phase_transition_alert: str = ""

    # SVG coordinates
    sd_x: int = 160; sd_y: int = 130; sd_y_label: int = 116
    lam_polyline: str = "40,130 280,130"
    lam_current_x: int = 280; lam_current_y: int = 130; lam_current_y_label: int = 116

    # §1 Exec summary
    exec_situation: str = ""
    exec_problem:   str = ""
    exec_action_1:  str = ""
    exec_action_2:  str = ""
    exec_action_3:  str = ""
    exec_impact:    str = ""
    exec_top_risk:  str = ""

    # §3 Intervention
    intervention_needed: str = "不要"
    delta_phi_need:      float = 0.0
    sensitivity_N:       float = 0.0
    sensitivity_T2:      float = 0.0
    sensitivity_theta:   float = 0.0
    priority_N:    float = 0.0
    priority_T2:   float = 0.0
    priority_theta: float = 0.0
    best_variable:   str = "N"
    delta_u_optimal: float = 0.0
    cost_N:    float = 1.0
    cost_T2:   float = 1.0
    cost_theta: float = 10.0
    lambda_val: float = 1.0
    total_interventions: int = 0

    # §4-2
    action_condition: Optional[ActionConditionDiagnosis] = None
    # §4-3
    edge_flows: List[EdgeFlowData] = field(default_factory=list)
    # §3 Phi graph
    phi_graph_data: Optional[PhiGraphData] = None
    # §4-4
    evidence_rows: str = ""
    # §5
    scenario_rows: str = ""
    # §6
    uncertainty_desc: str = ""
    additional_obs_rows: str = ""
    # §4-1
    mapping_nodes: str = ""
    mapping_links: str = ""
    mapping_field: str = ""
    mapping_boundary: str = ""
    mapping_regime: str = ""

    @classmethod
    def from_state(
        cls,
        state: StateSnapshot,
        phi_prev: Optional[float] = None,
        dt: float = 1.0,
        regime_embeddings: Optional[Dict[str, np.ndarray]] = None,
        edges: Optional[List[Tuple[str, str]]] = None,
        intervention_result=None,
        opt_lambda: float = 1.0,
        lambda_1_series: Optional[List[float]] = None,
        title: str = "SDFT v0.2 構造場分析レポート",
        subject: str = "分析対象",
        purpose: str = "構造場診断",
    ) -> "ReportContext":
        ctx = cls()
        ctx.title = title
        ctx.subject = subject
        ctx.purpose = purpose

        # Layer 1
        ctx.S_g = state.S_g
        ctx.D = state.D
        ctx.lambda_1 = state.lambda_1
        ctx.d_embed = state.d_embed

        # Layer 2
        ctx.F_geom = state.F_geom
        ctx.curvature = state.curvature
        ctx.Phi = state.Phi
        ctx.N = state.N
        ctx.mu = state.mu

        if phi_prev is None:
            phi_prev = state.Phi
        ctx.dPhi_dt = (state.Phi - phi_prev) / (dt + EPS)

        ctx.Phi_threshold, _ = auto_phi_threshold(state)
        tau_val, _ = auto_tau(state.lambda_1)
        ctx.Phi_predicted = state.Phi + ctx.dPhi_dt * tau_val

        # 𝓣
        if state.tension:
            t = state.tension
            ctx.T1 = t.T1; ctx.T1_evidence = t.T1_evidence
            ctx.T2 = t.T2; ctx.T2_evidence = t.T2_evidence
            ctx.T3 = t.T3; ctx.T3_evidence = t.T3_evidence
            ctx.T4 = t.T4; ctx.T4_evidence = t.T4_evidence
            ctx.T5 = t.T5; ctx.T5_evidence = t.T5_evidence

        # Phase
        lbl, css, expl = classify_phase(ctx.S_g, ctx.D, ctx.lambda_1, ctx.d_embed)
        ctx.phase_label = lbl
        ctx.phase_css = css
        ctx.phase_explanation = expl

        # SVG coordinates
        S_norm = _normalize_S_g(ctx.S_g)
        D_norm = _normalize_D(ctx.D, ctx.d_embed)
        ctx.sd_x, ctx.sd_y, ctx.sd_y_label = calc_sd_position(S_norm, D_norm)

        lam_series = lambda_1_series or [ctx.lambda_1]
        (ctx.lam_polyline, ctx.lam_current_x,
         ctx.lam_current_y, ctx.lam_current_y_label) = calc_lambda1_polyline(lam_series)

        # §3 Intervention
        if intervention_result and hasattr(intervention_result, 'steps'):
            steps = intervention_result.steps
            if steps:
                last = steps[-1]
                ctx.intervention_needed = "必要" if last.needs_intervention else "不要"
                ctx.delta_phi_need = last.delta_phi_need
                ctx.sensitivity_N = last.gradients.get('N', 0.0)
                ctx.sensitivity_T2 = last.gradients.get('T2', 0.0)
                ctx.sensitivity_theta = last.gradients.get('theta', 0.0)
                ctx.priority_N = last.priorities.get('N', 0.0)
                ctx.priority_T2 = last.priorities.get('T2', 0.0)
                ctx.priority_theta = last.priorities.get('theta', 0.0)
                ctx.best_variable = last.best_variable
                ctx.delta_u_optimal = last.delta_u_optimal
            ctx.total_interventions = intervention_result.total_interventions
        ctx.lambda_val = opt_lambda

        # §4-2
        ctx.action_condition = state.action_condition

        # §4-3
        if regime_embeddings and edges:
            ctx.edge_flows = calc_edge_flows(regime_embeddings, edges)

        # §3 Phi graph
        if intervention_result is not None:
            ctx.phi_graph_data = build_phi_graph_data(
                intervention_result, phi_threshold=float(ctx.Phi_threshold))

        # §4-4
        ctx.evidence_rows = _build_evidence_rows(ctx)

        return ctx


def _build_evidence_rows(ctx: "ReportContext") -> str:
    items = [
        ("S_g（軌道体積エントロピー）", f"{ctx.S_g:.4f}",
         "Derived", "k-NN ball 体積の平均対数"),
        ("D（相関次元）", f"{ctx.D:.4f}",
         "Derived", "Grassberger-Procaccia 法"),
        ("λ₁（最大リアプノフ指数）", f"{ctx.lambda_1:.4f}",
         "Derived", "Rosenstein 法"),
        ("F_geom = λ₂（スペクトルギャップ）", f"{ctx.F_geom:.6f}",
         "Derived", "正規化グラフ Laplacian Fiedler 値"),
        ("κ（Ollivier-Ricci 曲率）", f"{ctx.curvature:.4f}",
         "Derived", "W₁最適割当"),
        ("N（人数）", f"{ctx.N:.0f} 人",
         "Observed", "直接計測"),
        ("μ（化学ポテンシャル）", f"{ctx.mu:.4f}",
         "Derived", "T₂（測地距離）"),
        ("Φ（グランドポテンシャル）", f"{ctx.Phi:.4f}",
         "Derived", "F_geom - μN"),
        ("T₁（体積スケール差）", f"{ctx.T1:.4f}",
         ctx.T1_evidence, "|log Vol(M_A) - log Vol(M_B)|"),
        ("T₂（測地距離）", f"{ctx.T2:.4f}",
         ctx.T2_evidence, "Dijkstra on k-NN graph"),
        ("T₃（平行移動不整合）", f"{ctx.T3:.4f}",
         ctx.T3_evidence, "PCA フレーム回転角（簡易版）"),
        ("T₄（ホロノミー・曲率差）", f"{ctx.T4:.4f}",
         ctx.T4_evidence, "|κ_A - κ_B|/max(|κ_A|,|κ_B|)（簡易版）"),
        ("T₅（埋め込み歪み）", f"{ctx.T5:.4f}",
         ctx.T5_evidence, "1 - σ_min(J)/σ_max(J)"),
    ]
    ev_css = {
        "Observed": "ev-observed",
        "Derived": "ev-derived-t",
        "Assumed": "ev-assumed-t",
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


def _build_action_condition_rows(diag: Optional[ActionConditionDiagnosis]) -> str:
    if diag is None:
        return '<tr><td colspan="4" style="color:#aaa;">診断データなし</td></tr>'
    status_mark = {"good": "○", "warn": "△", "alert": "×", "mid": "－"}
    status_css = {"good": "good", "warn": "warn", "alert": "alert", "mid": "mid"}
    rows = []
    for row in diag.table_rows():
        mark = status_mark.get(row["status"], "－")
        css = status_css.get(row["status"], "mid")
        rows.append(
            f'<tr><td><span class="kpi-status {css}">{mark}</span> {row["item"]}</td>'
            f'<td style="font-size:.85rem;color:#555;">{row["variable"]}</td>'
            f'<td class="{css}" style="font-size:.85rem;">{row["status"].upper()}</td>'
            f'<td style="font-size:.85rem;">{row["interpretation"]}</td></tr>'
        )
    return "\n".join(rows)


def _default_scenario_rows(ctx: "ReportContext") -> str:
    rows = [
        ("λ₁（リアプノフ指数）",
         f"{ctx.lambda_1 - 0.05:.3f}（低下→安定化）",
         f"{ctx.lambda_1:.3f}（維持）",
         f"{ctx.lambda_1 + 0.05:.3f}（上昇→カオス化）"),
        ("S_g（体積エントロピー）",
         "低下（軌道コンパクト化）", "現状維持", "上昇（軌道拡散）"),
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
        '<tr><td>両レジームの同時観測</td><td>T₅ の精度向上</td><td>中</td></tr>'
    )


# ============================================================
# レポート生成
# ============================================================

def generate_report(
    ctx: "ReportContext",
    template_path: str,
    output_path: str = "",
) -> str:
    with open(template_path, 'r', encoding='utf-8') as f:
        html = f.read()

    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")

    # KPI status
    S_norm = _normalize_S_g(ctx.S_g)
    D_norm = _normalize_D(ctx.D, ctx.d_embed)
    s_css, s_txt = get_kpi_status('S_g', S_norm)
    d_css, d_txt = get_kpi_status('D', D_norm)
    lam_css, lam_txt = get_kpi_status('lambda_1', ctx.lambda_1)
    fg_css = "good" if ctx.F_geom > 0.01 else ("warn" if ctx.F_geom > 0.001 else "alert")
    fg_txt = f"λ₂={ctx.F_geom:.4f}" if ctx.F_geom > 0.001 else "低結合 ⚠"

    # Priority percentages
    prios = {'N': ctx.priority_N, 'T2': ctx.priority_T2, 'theta': ctx.priority_theta}
    max_prio = max(prios.values()) if max(prios.values()) > 0 else 1.0
    p_n_pct = int(ctx.priority_N / (max_prio + EPS) * 100)
    p_t2_pct = int(ctx.priority_T2 / (max_prio + EPS) * 100)
    p_th_pct = int(ctx.priority_theta / (max_prio + EPS) * 100)

    best_n_badge = '<span class="best-badge">★ 推奨</span>' if ctx.best_variable == 'N' else ''
    best_t2_badge = '<span class="best-badge">★ 推奨</span>' if ctx.best_variable == 'T2' else ''
    best_th_badge = '<span class="best-badge">★ 推奨</span>' if ctx.best_variable == 'theta' else ''

    alert_html = ""
    if ctx.phase_transition_alert:
        alert_html = f'<div class="phase-alert-box">{ctx.phase_transition_alert}</div>'

    edge_hm_html = render_edge_heatmap_table_html(ctx.edge_flows)
    edge_svg = render_edge_heatmap_svg(ctx.edge_flows)

    acd_rows = _build_action_condition_rows(ctx.action_condition)
    acd_overall_css = ""
    acd_overall_status = "—"
    acd_overall_msg = ""
    if ctx.action_condition:
        acd_overall_css = ctx.action_condition.overall_status
        acd_overall_status = {
            "good": "○ 行動が起きやすい", "warn": "△ 一部に障壁あり",
            "alert": "× 行動が起きにくい", "mid": "― 中程度",
        }.get(ctx.action_condition.overall_status, "—")
        acd_overall_msg = ctx.action_condition.overall_message

    exec_situation = ctx.exec_situation or f"現在のフェーズは {ctx.phase_label}。"
    exec_problem = ctx.exec_problem or "詳細な分析を要します。"
    exec_action_1 = ctx.exec_action_1 or "データの継続的な観測を実施する。"
    exec_action_2 = ctx.exec_action_2 or "λ₁（リアプノフ指数）の推移を監視する。"
    exec_action_3 = ctx.exec_action_3 or f"最優先介入変数（{ctx.best_variable}）への施策を検討する。"
    exec_impact = ctx.exec_impact or "Φ の維持・向上が期待される。"
    exec_top_risk = ctx.exec_top_risk or "相転移への移行リスク。"

    replacements = {
        "{{REPORT_TITLE}}": ctx.title,
        "{{ANALYSIS_DATE}}": now,
        "{{SUBJECT}}": ctx.subject,
        "{{PHASE_LABEL}}": ctx.phase_label,
        "{{ANALYSIS_PURPOSE}}": ctx.purpose,
        # §1
        "{{EXEC_SITUATION}}": exec_situation,
        "{{EXEC_PROBLEM}}": exec_problem,
        "{{EXEC_ACTION_1}}": exec_action_1,
        "{{EXEC_ACTION_2}}": exec_action_2,
        "{{EXEC_ACTION_3}}": exec_action_3,
        "{{EXEC_IMPACT}}": exec_impact,
        "{{EXEC_TOP_RISK}}": exec_top_risk,
        # §2 KPI
        "{{PHASE_CSS}}": ctx.phase_css,
        "{{PHASE_EXPLANATION}}": ctx.phase_explanation,
        "{{SG_VALUE}}": f"{ctx.S_g:.4f}",
        "{{SG_STATUS_CSS}}": s_css,
        "{{SG_STATUS}}": s_txt,
        "{{D_VALUE}}": f"{ctx.D:.4f}",
        "{{D_STATUS_CSS}}": d_css,
        "{{D_STATUS}}": d_txt,
        "{{LAMBDA1_VALUE}}": f"{ctx.lambda_1:.4f}",
        "{{LAMBDA1_STATUS_CSS}}": lam_css,
        "{{LAMBDA1_STATUS}}": lam_txt,
        "{{FGEOM_VALUE}}": f"{ctx.F_geom:.6f}",
        "{{FGEOM_STATUS_CSS}}": fg_css,
        "{{FGEOM_STATUS}}": fg_txt,
        "{{CURVATURE_VALUE}}": f"{ctx.curvature:.4f}",
        "{{MU_VALUE}}": f"{ctx.mu:.4f}",
        "{{PHI_VALUE}}": f"{ctx.Phi:.4f}",
        "{{PHI_DELTA}}": f"{ctx.dPhi_dt:.4f}",
        "{{PHI_THRESHOLD}}": f"{ctx.Phi_threshold:.4f}",
        "{{PHI_PREDICTED}}": f"{ctx.Phi_predicted:.4f}",
        "{{N_VALUE}}": f"{ctx.N:.0f}",
        "{{T1_VALUE}}": f"{ctx.T1:.4f}",
        "{{T1_EVIDENCE}}": ctx.T1_evidence,
        "{{T2_VALUE}}": f"{ctx.T2:.4f}",
        "{{T2_EVIDENCE}}": ctx.T2_evidence,
        "{{T3_VALUE}}": f"{ctx.T3:.4f}",
        "{{T3_EVIDENCE}}": ctx.T3_evidence,
        "{{T4_VALUE}}": f"{ctx.T4:.4f}",
        "{{T4_EVIDENCE}}": ctx.T4_evidence,
        "{{T5_VALUE}}": f"{ctx.T5:.4f}",
        "{{T5_EVIDENCE}}": ctx.T5_evidence,
        "{{PHASE_TRANSITION_ALERT}}": alert_html,
        # §3
        "{{INTERVENTION_NEEDED}}": ctx.intervention_needed,
        "{{DELTA_PHI_NEED}}": f"{ctx.delta_phi_need:.4f}",
        "{{SENSITIVITY_N}}": f"{ctx.sensitivity_N:.4f}",
        "{{SENSITIVITY_T2}}": f"{ctx.sensitivity_T2:.4f}",
        "{{SENSITIVITY_THETA}}": f"{ctx.sensitivity_theta:.4e}",
        "{{PRIORITY_N}}": f"{ctx.priority_N:.4f}",
        "{{PRIORITY_T2}}": f"{ctx.priority_T2:.4f}",
        "{{PRIORITY_THETA}}": f"{ctx.priority_theta:.4f}",
        "{{PRIORITY_N_PCT}}": str(p_n_pct),
        "{{PRIORITY_T2_PCT}}": str(p_t2_pct),
        "{{PRIORITY_THETA_PCT}}": str(p_th_pct),
        "{{BEST_VARIABLE}}": ctx.best_variable,
        "{{DELTA_U_OPTIMAL}}": f"{ctx.delta_u_optimal:.4f}",
        "{{BEST_N_BADGE}}": best_n_badge,
        "{{BEST_T2_BADGE}}": best_t2_badge,
        "{{BEST_THETA_BADGE}}": best_th_badge,
        "{{COST_N}}": f"{ctx.cost_N:.1f}",
        "{{COST_T2}}": f"{ctx.cost_T2:.1f}",
        "{{COST_THETA}}": f"{ctx.cost_theta:.1f}",
        # §4-2
        "{{ACTION_CONDITION_ROWS}}": acd_rows,
        "{{ACTION_OVERALL_CSS}}": acd_overall_css,
        "{{ACTION_OVERALL_STATUS}}": acd_overall_status,
        "{{ACTION_OVERALL_MESSAGE}}": acd_overall_msg,
        # §4-3
        "{{EDGE_HEATMAP_TABLE}}": edge_hm_html,
        "{{EDGE_HEATMAP_SVG}}": edge_svg,
        # §3 Phi graph
        "{{PHI_ENERGY_GRAPH_SVG}}": render_phi_energy_graph_svg(ctx.phi_graph_data)
            if ctx.phi_graph_data else '<p style="color:#aaa;">データなし</p>',
        "{{PHI_INTERVENTION_TABLE}}": render_phi_intervention_table_html(ctx.phi_graph_data)
            if ctx.phi_graph_data else '<p style="color:#aaa;">データなし</p>',
        # §4-4
        "{{EVIDENCE_ROWS}}": ctx.evidence_rows,
        # SVG phase diagrams
        "{{SD_X}}": str(ctx.sd_x),
        "{{SD_Y}}": str(ctx.sd_y),
        "{{SD_Y_LABEL}}": str(ctx.sd_y_label),
        "{{LAM_POLYLINE}}": ctx.lam_polyline,
        "{{LAM_CURRENT_X}}": str(ctx.lam_current_x),
        "{{LAM_CURRENT_Y}}": str(ctx.lam_current_y),
        "{{LAM_CURRENT_Y_LABEL}}": str(ctx.lam_current_y_label),
        # §5
        "{{SCENARIO_ROWS}}": ctx.scenario_rows or _default_scenario_rows(ctx),
        # §6
        "{{UNCERTAINTY_DESC}}": ctx.uncertainty_desc or "追加データによる精度向上が望ましい。",
        "{{ADDITIONAL_OBS_ROWS}}": ctx.additional_obs_rows or _default_obs_rows(),
        # §4-1
        "{{MAPPING_NODES}}": ctx.mapping_nodes or "分析対象の構成要素",
        "{{MAPPING_LINKS}}": ctx.mapping_links or "構成要素間の接続",
        "{{MAPPING_FIELD}}": ctx.mapping_field or "外部環境",
        "{{MAPPING_BOUNDARY}}": ctx.mapping_boundary or "レジームの境界",
        "{{MAPPING_REGIME}}": ctx.mapping_regime or "同一オペモードを持つ領域（A / B）",
    }

    for placeholder, value in replacements.items():
        html = html.replace(placeholder, str(value))

    if output_path:
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(html)

    return html
