"""
sdft_v02_phi_graph.py
SDFT revised v0.2 — §3 Φ エネルギー-時刻グラフ

v0.0.6 からの変更点：
  ・Y軸ラベルを「Φ（グランドポテンシャル）」に変更（J 単位を削除）
  ・操作変数ラベルを多様体版に更新（T₂→測地距離）
  ・それ以外の SVG 描画ロジックはそのまま維持
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Dict, Tuple, Optional

EPS = 1e-12

_VAR_COLORS = {'N': '#27ae60', 'T2': '#e67e22', 'theta': '#8e44ad'}
_VAR_LABELS = {'N': 'N', 'T2': 'T₂', 'theta': 'θ'}
_VAR_LABELS_JP = {
    'N':     'u₁ = N（人数）',
    'T2':    'u₂ = T₂（測地距離）',
    'theta': 'u₃ = θ（構造）',
}


@dataclass
class PhiGraphData:
    phi_trajectory:      List[float]
    phi_threshold:       float
    total_interventions: int
    intervention_steps:  List[dict]


def build_phi_graph_data(opt_result, phi_threshold: float) -> PhiGraphData:
    steps_detail = []
    for step in opt_result.steps:
        if step.needs_intervention:
            steps_detail.append({
                't':              step.t,
                'phi_before':     step.phi_before,
                'phi_after':      step.phi_after,
                'phi_jump':       step.phi_after - step.phi_before,
                'phi_predicted':  step.phi_predicted,
                'delta_u':        step.delta_u_optimal,
                'best_var':       step.best_variable,
                'var_label':      _VAR_LABELS.get(step.best_variable, step.best_variable),
                'priority_N':     step.priorities.get('N', 0.0),
                'priority_T2':    step.priorities.get('T2', 0.0),
                'priority_theta': step.priorities.get('theta', 0.0),
                'delta_phi_need': step.delta_phi_need,
            })
    return PhiGraphData(
        phi_trajectory=opt_result.phi_trajectory,
        phi_threshold=phi_threshold,
        total_interventions=opt_result.total_interventions,
        intervention_steps=steps_detail,
    )


def render_phi_energy_graph_svg(
    graph_data: PhiGraphData,
    width: int = 700,
    height: int = 340,
) -> str:
    traj = graph_data.phi_trajectory
    n = len(traj)

    if n == 0:
        return (
            f'<svg viewBox="0 0 {width} {height}" xmlns="http://www.w3.org/2000/svg">'
            f'<text x="{width//2}" y="{height//2}" text-anchor="middle" '
            f'fill="#aaa" font-size="13">データなし</text></svg>'
        )

    PAD_L, PAD_R, PAD_T, PAD_B = 72, 24, 48, 52
    W = width - PAD_L - PAD_R
    H = height - PAD_T - PAD_B

    all_phi = list(traj) + [graph_data.phi_threshold]
    for s in graph_data.intervention_steps:
        all_phi += [s['phi_before'], s['phi_after']]

    phi_min_raw = min(all_phi)
    phi_max_raw = max(all_phi)
    phi_span = phi_max_raw - phi_min_raw
    if phi_span < EPS:
        phi_span = 1.0
    phi_min = phi_min_raw - phi_span * 0.12
    phi_max = phi_max_raw + phi_span * 0.18

    def tx(t):
        return PAD_L + int((t / max(n - 1, 1)) * W)

    def ty(phi):
        frac = (phi - phi_min) / (phi_max - phi_min + EPS)
        return PAD_T + H - int(frac * H)

    lines = [
        f'<svg viewBox="0 0 {width} {height}" xmlns="http://www.w3.org/2000/svg" '
        f'style="background:#fff;border-radius:10px;border:1px solid #e0e4ed;">',
        '<defs>',
        '<marker id="arr_up" markerWidth="8" markerHeight="8" refX="4" refY="4" orient="auto">',
        '<path d="M1,7 L4,1 L7,7" fill="none" stroke="#27ae60" stroke-width="1.5"/>',
        '</marker>',
        '</defs>',
    ]

    # Grid + Y axis
    n_yticks = 5
    for i in range(n_yticks + 1):
        phi_g = phi_min + (phi_max - phi_min) * i / n_yticks
        y_g = ty(phi_g)
        lines.append(
            f'<line x1="{PAD_L}" y1="{y_g}" x2="{PAD_L+W}" y2="{y_g}" '
            f'stroke="#edf0f5" stroke-width="1"/>'
        )
        lines.append(
            f'<text x="{PAD_L-6}" y="{y_g+4}" text-anchor="end" '
            f'font-size="9" fill="#888">{phi_g:.1f}</text>'
        )

    # Axes
    lines += [
        f'<line x1="{PAD_L}" y1="{PAD_T}" x2="{PAD_L}" y2="{PAD_T+H}" '
        f'stroke="#7f8c8d" stroke-width="1.5"/>',
        f'<line x1="{PAD_L}" y1="{PAD_T+H}" x2="{PAD_L+W}" y2="{PAD_T+H}" '
        f'stroke="#7f8c8d" stroke-width="1.5"/>',
        f'<text x="11" y="{PAD_T+H//2}" text-anchor="middle" font-size="10" '
        f'fill="#555" transform="rotate(-90,11,{PAD_T+H//2})">'
        f'Φ（グランドポテンシャル）</text>',
        f'<text x="{PAD_L+W//2}" y="{height-8}" text-anchor="middle" '
        f'font-size="10" fill="#555">時刻 t（ステップ）</text>',
    ]

    tick_step = max(1, n // 8)
    for i in range(0, n, tick_step):
        x_i = tx(i)
        lines.append(
            f'<text x="{x_i}" y="{PAD_T+H+14}" text-anchor="middle" '
            f'font-size="9" fill="#888">{i}</text>'
        )

    # No-intervention estimate (gray dashed)
    if n >= 2:
        dPhi_init = traj[1] - traj[0]
        no_intv = [traj[0] + dPhi_init * i for i in range(n)]
        pts_gray = " ".join(
            f"{tx(i)},{ty(max(phi_min, min(phi_max, v)))}"
            for i, v in enumerate(no_intv)
        )
        lines.append(
            f'<polyline points="{pts_gray}" fill="none" stroke="#bdc3c7" '
            f'stroke-width="1.5" stroke-dasharray="6,3" opacity="0.8"/>'
        )

    # Threshold (red dashed)
    y_thresh = ty(graph_data.phi_threshold)
    if PAD_T <= y_thresh <= PAD_T + H:
        lines += [
            f'<line x1="{PAD_L}" y1="{y_thresh}" x2="{PAD_L+W}" y2="{y_thresh}" '
            f'stroke="#e74c3c" stroke-width="1.8" stroke-dasharray="10,4"/>',
            f'<rect x="{PAD_L+5}" y="{y_thresh-16}" width="80" height="14" '
            f'rx="3" fill="#e74c3c" opacity="0.9"/>',
            f'<text x="{PAD_L+10}" y="{y_thresh-5}" font-size="9" '
            f'font-weight="700" fill="white">Φ₀ 介入閾値</text>',
        ]

    # Actual trajectory (blue)
    pts_blue = " ".join(
        f"{tx(i)},{ty(max(phi_min, min(phi_max, v)))}"
        for i, v in enumerate(traj)
    )
    lines.append(
        f'<polygon points="{tx(0)},{PAD_T+H} {pts_blue} {tx(n-1)},{PAD_T+H}" '
        f'fill="#2d5f8a" opacity="0.07"/>'
    )
    lines.append(
        f'<polyline points="{pts_blue}" fill="none" stroke="#2d5f8a" '
        f'stroke-width="2.5" stroke-linejoin="round"/>'
    )

    # Intervention points
    for idx, s in enumerate(graph_data.intervention_steps):
        t_i = int(s['t'])
        x_i = tx(t_i)
        y_before = ty(max(phi_min, min(phi_max, s['phi_before'])))
        y_after = ty(max(phi_min, min(phi_max, s['phi_after'])))
        jump = s['phi_jump']
        color = _VAR_COLORS.get(s['best_var'], '#27ae60')
        v_lbl = _VAR_LABELS.get(s['best_var'], s['best_var'])

        lines.append(
            f'<line x1="{x_i}" y1="{PAD_T}" x2="{x_i}" y2="{PAD_T+H}" '
            f'stroke="{color}" stroke-width="1.2" stroke-dasharray="3,3" opacity="0.5"/>'
        )
        lines.append(
            f'<circle cx="{x_i}" cy="{y_before}" r="5" fill="white" '
            f'stroke="{color}" stroke-width="2"/>'
        )
        if abs(y_after - y_before) > 4:
            arrow_x = x_i + 8
            y_top = min(y_before, y_after)
            y_bot = max(y_before, y_after)
            lines.append(
                f'<line x1="{arrow_x}" y1="{y_bot-2}" '
                f'x2="{arrow_x}" y2="{y_top+8}" '
                f'stroke="{color}" stroke-width="2.5" '
                f'marker-end="url(#arr_up)"/>'
            )
        lines.append(
            f'<circle cx="{x_i}" cy="{y_after}" r="6" fill="{color}" opacity="0.9"/>'
        )

        box_w = 120
        box_h = 36
        lx = min(x_i + 14, width - box_w - PAD_R)
        ly = max(PAD_T + 4, min(y_after - box_h - 4, PAD_T + H - box_h - 4))
        jump_str = f"ΔΦ={'+' if jump >= 0 else ''}{jump:.3f}"
        du_str = f"Δ{v_lbl}={abs(s['delta_u']):.3f}"

        lines += [
            f'<rect x="{lx}" y="{ly}" width="{box_w}" height="{box_h}" '
            f'rx="5" fill="{color}" opacity="0.92"/>',
            f'<text x="{lx+6}" y="{ly+13}" font-size="9" font-weight="700" '
            f'fill="white">介入{idx+1}：{du_str}</text>',
            f'<text x="{lx+6}" y="{ly+26}" font-size="9" fill="white">{jump_str}</text>',
        ]

    # Current value marker
    cx_now = tx(n - 1)
    cy_now = ty(max(phi_min, min(phi_max, traj[-1])))
    lines += [
        f'<circle cx="{cx_now}" cy="{cy_now}" r="6" fill="#1a3a5c"/>',
        f'<text x="{cx_now-30}" y="{cy_now-10}" font-size="9" '
        f'font-weight="700" fill="#1a3a5c">現在 {traj[-1]:.2f}</text>',
    ]

    # Legend
    ley = PAD_T + 4
    legend = [
        ('#2d5f8a', '実際の Φ', False),
        ('#bdc3c7', '介入なし推定', True),
        ('#e74c3c', '介入閾値 Φ₀', True),
    ]
    for li, (col, lbl, dashed) in enumerate(legend):
        lx2 = PAD_L + li * 140
        dash = 'stroke-dasharray="6,3"' if dashed else ''
        lines += [
            f'<line x1="{lx2}" y1="{ley+5}" x2="{lx2+18}" y2="{ley+5}" '
            f'stroke="{col}" stroke-width="2" {dash}/>',
            f'<text x="{lx2+22}" y="{ley+9}" font-size="9" fill="#555">{lbl}</text>',
        ]
    var_leg_y = ley + 18
    for vi, (var, col) in enumerate(_VAR_COLORS.items()):
        lx3 = PAD_L + vi * 130
        lbl3 = _VAR_LABELS_JP.get(var, var)
        lines += [
            f'<circle cx="{lx3+5}" cy="{var_leg_y+4}" r="5" fill="{col}" opacity="0.9"/>',
            f'<text x="{lx3+14}" y="{var_leg_y+8}" font-size="9" fill="#555">{lbl3}</text>',
        ]

    lines.append('</svg>')
    return "\n".join(lines)


def render_phi_intervention_table_html(graph_data: PhiGraphData) -> str:
    if not graph_data.intervention_steps:
        return (
            '<div style="background:#eafaf1;border:1px solid #a9dfbf;'
            'border-radius:8px;padding:14px 18px;font-size:.9rem;color:#1e8449;">'
            '✅ 介入不要：Φ が閾値 Φ₀ を上回っています。'
            '</div>'
        )

    rows = []
    for i, s in enumerate(graph_data.intervention_steps):
        col = _VAR_COLORS.get(s['best_var'], '#555')
        v_jp = _VAR_LABELS_JP.get(s['best_var'], s['best_var'])
        jump = s['phi_jump']
        jcol = '#27ae60' if jump >= 0 else '#e74c3c'
        jstr = f"+{jump:.4f}" if jump >= 0 else f"{jump:.4f}"

        prios = {
            'N': s['priority_N'], 'T2': s['priority_T2'],
            'theta': s['priority_theta'],
        }
        max_p = max(prios.values()) if max(prios.values()) > 0 else 1.0

        def pbar(k):
            pct = int(prios[k] / max_p * 72)
            c = _VAR_COLORS.get(k, '#555')
            lbl = _VAR_LABELS.get(k, k)
            return (
                f'<div style="display:flex;align-items:center;gap:3px;margin-bottom:2px;">'
                f'<span style="width:14px;font-size:.76rem;color:{c};">{lbl}</span>'
                f'<div style="width:72px;background:#e0e4ed;border-radius:3px;height:5px;">'
                f'<div style="width:{pct}px;background:{c};border-radius:3px;height:5px;"></div>'
                f'</div>'
                f'<span style="font-size:.74rem;color:#888;min-width:30px;">{prios[k]:.2f}</span>'
                f'</div>'
            )

        rows.append(
            f'<tr>'
            f'<td style="text-align:center;font-weight:700;color:#1a3a5c;">t = {int(s["t"])}</td>'
            f'<td><span style="background:{col};color:white;border-radius:4px;'
            f'padding:2px 8px;font-size:.82rem;font-weight:700;">{v_jp}</span></td>'
            f'<td style="font-size:.88rem;text-align:right;">{abs(s["delta_u"]):.4f}</td>'
            f'<td style="font-size:.88rem;text-align:right;color:{jcol};font-weight:700;">{jstr}</td>'
            f'<td style="font-size:.82rem;text-align:right;">{s["phi_before"]:.4f}</td>'
            f'<td style="font-size:.82rem;text-align:right;color:{col};font-weight:600;">{s["phi_after"]:.4f}</td>'
            f'<td>{pbar("N")}{pbar("T2")}{pbar("theta")}</td>'
            f'</tr>'
        )

    header = (
        '<table class="evidence-table" style="margin-top:12px;">'
        '<thead><tr>'
        '<th style="text-align:center;">介入時刻</th>'
        '<th>推奨操作変数（★最優先）</th>'
        '<th style="text-align:right;">最小 Δu*</th>'
        '<th style="text-align:right;">ΔΦ</th>'
        '<th style="text-align:right;">Φ（介入前）</th>'
        '<th style="text-align:right;">Φ（介入後）</th>'
        '<th>priority_score（相対）</th>'
        '</tr></thead><tbody>'
    )
    footer = (
        '</tbody></table>'
        '<p style="font-size:.8rem;color:#888;margin-top:6px;">'
        '最小 Δu*：Φ₀ を回復するために必要な最小操作量。'
        'ΔΦ：介入によるグランドポテンシャルの変化量。</p>'
    )
    return header + "\n".join(rows) + footer
