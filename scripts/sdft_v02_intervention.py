"""
sdft_v02_intervention.py
SDFT revised v0.2 — 介入最適化ループと行動条件診断

v0.0.6 からの変更点：
  ① λ₁（リアプノフ指数）が H（Hurst 指数）を置換
  ② μ = T₂（測地距離）が μ = 2√T₂ を置換
  ③ F_geom = λ₂（スペクトルギャップ）が F = U - TS を置換
  ④ λ スキャン no-op バグ修正：介入発火閾値に λ を反映
  ⑤ auto_phi_threshold ratio バグ修正：符号不変式に変更
  ⑥ dΦ/dt 初期値バグ修正：最初のステップは phi_prev = Phi で初期化

【行動条件診断の6項目】
  T₂       → 到達可能性（測地距離 = 移動コスト μ の低さ）
  T₅       → 伝達効率（埋め込み歪みの小ささ）
  T₃       → 方向整合性（平行移動不整合の小ささ）
  T₁       → スケール適合（体積スケール差の小ささ）
  dΦ/dt    → 変化の好機（Φ が減少中＝変化の機が熟している）
  λ₁       → 文脈の持続性（λ₁ < 0 で軌道が安定・予測可能）

【介入最適化の定式化】
  min  Σ‖u_k‖² + λK
  s.t. Φ(t) ≥ Φ₀  ∀t

  操作変数：
    u₁ = N    ∂Φ/∂N  = -μ = -T₂
    u₂ = T₂   ∂Φ/∂T₂ = -N
    u₃ = θ    ∂Φ/∂θ  = ∂F_geom/∂θ ≈ 数値微分
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Optional

import numpy as np

from sdft_v02_embedding import Evidence, EPS
from sdft_v02_tension import RegimeTensionVector

# ============================================================
# 行動条件診断
# ============================================================

@dataclass
class ActionConditionDiagnosis:
    """§4-2 行動条件診断。既存変数の読み替えのみで構成。"""
    T1:      float
    T2:      float
    T3:      float
    T5:      float
    mu:      float
    lambda_1: float
    dPhi_dt: float

    status_reachability:  str = ""
    status_transmission:  str = ""
    status_alignment:     str = ""
    status_scale:         str = ""
    status_timing:        str = ""
    status_persistence:   str = ""

    interp_reachability:  str = ""
    interp_transmission:  str = ""
    interp_alignment:     str = ""
    interp_scale:         str = ""
    interp_timing:        str = ""
    interp_persistence:   str = ""

    overall_status:  str = ""
    overall_message: str = ""

    def __post_init__(self):
        self._diagnose()

    def _diagnose(self):
        # ① 到達可能性（μ = T₂ 測地距離）
        if self.mu < 0.5:
            self.status_reachability = "good"
            self.interp_reachability = f"測地距離 μ={self.mu:.3f} — 近い。レジーム間の移動が容易。"
        elif self.mu < 1.5:
            self.status_reachability = "warn"
            self.interp_reachability = f"測地距離 μ={self.mu:.3f} — 中程度。一定の距離がある。"
        else:
            self.status_reachability = "alert"
            self.interp_reachability = f"測地距離 μ={self.mu:.3f} — 遠い。多様体上の距離が大きい。"

        # ② 伝達効率（T₅ 埋め込み歪み）
        if self.T5 < 0.3:
            self.status_transmission = "good"
            self.interp_transmission = f"埋め込み歪み T₅={self.T5:.3f} — 低い。写像が等長に近い。"
        elif self.T5 < 0.6:
            self.status_transmission = "warn"
            self.interp_transmission = f"埋め込み歪み T₅={self.T5:.3f} — 中程度。情報の一部が潰れている。"
        else:
            self.status_transmission = "alert"
            self.interp_transmission = f"埋め込み歪み T₅={self.T5:.3f} — 高い。写像のランク落ちが大きい。"

        # ③ 方向整合性（T₃ 平行移動不整合）
        if self.T3 < 0.3:
            self.status_alignment = "good"
            self.interp_alignment = f"フレーム回転 T₃={self.T3:.3f} — 小さい。双方向に自然な対応。"
        elif self.T3 < 0.8:
            self.status_alignment = "warn"
            self.interp_alignment = f"フレーム回転 T₃={self.T3:.3f} — 中程度。局所構造に捻れがある。"
        else:
            self.status_alignment = "alert"
            self.interp_alignment = f"フレーム回転 T₃={self.T3:.3f} — 大きい。局所フレームが大きく回転している。"

        # ④ スケール適合（T₁ 体積スケール差）
        if self.T1 < 0.3:
            self.status_scale = "good"
            self.interp_scale = f"体積差 T₁={self.T1:.3f} — 小さい。スケールが近い。"
        elif self.T1 < 1.0:
            self.status_scale = "warn"
            self.interp_scale = f"体積差 T₁={self.T1:.3f} — 中程度。スケールのズレあり。"
        else:
            self.status_scale = "alert"
            self.interp_scale = f"体積差 T₁={self.T1:.3f} — 大きい。空間の大きさが大幅に異なる。"

        # ⑤ 変化の好機（dΦ/dt）
        if self.dPhi_dt < -0.001:
            self.status_timing = "good"
            self.interp_timing = f"dΦ/dt={self.dPhi_dt:.4f} — Φ が減少中。変化の機が熟している。"
        elif self.dPhi_dt < 0.0:
            self.status_timing = "warn"
            self.interp_timing = f"dΦ/dt={self.dPhi_dt:.4f} — わずかに減少。好機の兆し。"
        else:
            self.status_timing = "mid"
            self.interp_timing = f"dΦ/dt={self.dPhi_dt:.4f} — Φ が安定または上昇中。"

        # ⑥ 文脈の持続性（λ₁ リアプノフ指数）
        if self.lambda_1 < -0.05:
            self.status_persistence = "good"
            self.interp_persistence = f"λ₁={self.lambda_1:.3f} — 安定（軌道が収束）。予測可能で持続的。"
        elif self.lambda_1 < 0.01:
            self.status_persistence = "warn"
            self.interp_persistence = f"λ₁={self.lambda_1:.3f} — 中立（限周期的）。やや不安定。"
        else:
            self.status_persistence = "alert"
            self.interp_persistence = f"λ₁={self.lambda_1:.3f} — カオス的（軌道が発散）。文脈が維持されない。"

        # 総合判定
        statuses = [
            self.status_reachability, self.status_transmission,
            self.status_alignment, self.status_scale,
            self.status_timing, self.status_persistence,
        ]
        alert_count = statuses.count("alert")
        good_count  = statuses.count("good")

        if alert_count == 0 and good_count >= 4:
            self.overall_status  = "good"
            self.overall_message = "行動が起きやすい状態。全指標が良好。"
        elif alert_count >= 3:
            self.overall_status  = "alert"
            self.overall_message = "行動が起きにくい状態。複数の障壁が同時に存在。"
        elif alert_count >= 1:
            self.overall_status  = "warn"
            alert_items = []
            if self.status_reachability == "alert":
                alert_items.append("測地距離が大きい（T₂）")
            if self.status_transmission == "alert":
                alert_items.append("埋め込み歪みが大きい（T₅）")
            if self.status_alignment == "alert":
                alert_items.append("フレーム回転が大きい（T₃）")
            if self.status_scale == "alert":
                alert_items.append("体積スケール差が大きい（T₁）")
            if self.status_persistence == "alert":
                alert_items.append("軌道がカオス的（λ₁）")
            self.overall_message = "一部に障壁あり。優先課題：" + "・".join(alert_items)
        else:
            self.overall_status  = "mid"
            self.overall_message = "行動条件は中程度。特定の介入で改善の余地がある。"

    def table_rows(self) -> List[Dict]:
        return [
            {"item": "到達可能性",   "variable": f"μ = T₂ = {self.mu:.4f}",
             "status": self.status_reachability, "interpretation": self.interp_reachability},
            {"item": "伝達効率",     "variable": f"T₅ = {self.T5:.4f}",
             "status": self.status_transmission, "interpretation": self.interp_transmission},
            {"item": "方向整合性",   "variable": f"T₃ = {self.T3:.4f}",
             "status": self.status_alignment, "interpretation": self.interp_alignment},
            {"item": "スケール適合", "variable": f"T₁ = {self.T1:.4f}",
             "status": self.status_scale, "interpretation": self.interp_scale},
            {"item": "変化の好機",   "variable": f"dΦ/dt = {self.dPhi_dt:.4f}",
             "status": self.status_timing, "interpretation": self.interp_timing},
            {"item": "文脈の持続性", "variable": f"λ₁ = {self.lambda_1:.4f}",
             "status": self.status_persistence, "interpretation": self.interp_persistence},
        ]

    def summary(self) -> str:
        lines = [f"=== 行動条件診断 ===",
                 f"総合：{self.overall_status.upper()}  {self.overall_message}", ""]
        for row in self.table_rows():
            mark = {"good": "○", "warn": "△", "alert": "×", "mid": "－"}.get(row["status"], "－")
            lines.append(f"  {mark} {row['item']:10s}  {row['variable']:20s}")
            lines.append(f"    {row['interpretation']}")
        return "\n".join(lines)


def calc_action_condition(
    T1: float, T2: float, T3: float, T5: float,
    lambda_1: float,
    phi_current: float, phi_prev: float, dt: float = 1.0,
) -> ActionConditionDiagnosis:
    from sdft_v02_potential import calc_mu
    mu = calc_mu(T2)
    dPhi_dt = (phi_current - phi_prev) / (dt + EPS)
    return ActionConditionDiagnosis(
        T1=T1, T2=T2, T3=T3, T5=T5, mu=mu,
        lambda_1=lambda_1, dPhi_dt=dPhi_dt,
    )


# ============================================================
# StateSnapshot
# ============================================================

@dataclass
class StateSnapshot:
    """時刻 t における SDFT v0.2 変数の完全スナップショット。"""
    t:            float
    embedded_A:   np.ndarray
    embedded_B:   np.ndarray
    N:            float = 0.0

    # Layer 1
    S_g:       float = 0.0
    D:         float = 0.0
    lambda_1:  float = 0.0
    d_embed:   int   = 3

    # Layer 2
    F_geom:    float = 0.0
    curvature: float = 0.0

    # 𝓣
    tension: Optional[RegimeTensionVector] = None

    # Φ
    mu:  float = 0.0
    Phi: float = 0.0

    # 行動条件診断
    action_condition: Optional[ActionConditionDiagnosis] = None

    def compute(self, phi_prev: Optional[float] = None, dt: float = 1.0) -> "StateSnapshot":
        """全変数を計算して自己を返す。"""
        from sdft_v02_embedding import calc_volume_entropy, calc_correlation_dim, calc_lyapunov
        from sdft_v02_geometry import calc_spectral_gap, calc_ollivier_ricci
        from sdft_v02_tension import calc_regime_tension_vector
        from sdft_v02_potential import calc_mu, calc_grand_potential

        # Layer 1（レジーム A の軌道構造）
        self.S_g      = calc_volume_entropy(self.embedded_A)
        self.D        = calc_correlation_dim(self.embedded_A)
        self.lambda_1 = calc_lyapunov(self.embedded_A)
        self.d_embed  = self.embedded_A.shape[1]

        # Layer 2（レジーム A の幾何構造）
        self.F_geom    = calc_spectral_gap(self.embedded_A)
        self.curvature = calc_ollivier_ricci(self.embedded_A, n_edges=100)

        # 𝓣
        self.tension = calc_regime_tension_vector(self.embedded_A, self.embedded_B)

        # Φ
        self.mu  = calc_mu(self.tension.T2)
        self.Phi = calc_grand_potential(self.F_geom, self.mu, self.N)

        # 行動条件診断
        # ⑥ dΦ/dt バグ修正：phi_prev が None なら Phi で初期化（dΦ/dt = 0）
        if phi_prev is None:
            phi_prev = self.Phi

        self.action_condition = calc_action_condition(
            T1=self.tension.T1, T2=self.tension.T2,
            T3=self.tension.T3, T5=self.tension.T5,
            lambda_1=self.lambda_1,
            phi_current=self.Phi, phi_prev=phi_prev, dt=dt,
        )
        return self


# ============================================================
# 介入ステップ
# ============================================================

@dataclass
class InterventionStep:
    t:                  float
    phi_before:         float
    phi_predicted:      float
    phi_threshold:      float
    needs_intervention: bool
    delta_phi_need:     float

    gradients:     Dict[str, float] = field(default_factory=dict)
    efficiencies:  Dict[str, float] = field(default_factory=dict)
    priorities:    Dict[str, float] = field(default_factory=dict)

    best_variable:   str   = ""
    delta_u_optimal: float = 0.0
    phi_after:       float = 0.0

    action_condition: Optional[ActionConditionDiagnosis] = None


def compute_intervention_step(
    state: StateSnapshot,
    phi_prev: float,
    dt: float, tau: float,
    phi_threshold: float,
    costs: Dict[str, float],
    lambda_count: float = 1.0,
) -> InterventionStep:
    """1 ステップの介入計算。"""
    phi_current   = state.Phi
    dPhi_dt       = (phi_current - phi_prev) / (dt + EPS)
    phi_predicted = phi_current + dPhi_dt * tau

    # ④ λ バグ修正：λ を介入発火の余裕量に反映
    effective_threshold = phi_threshold + lambda_count * 0.01
    needs = phi_predicted < effective_threshold
    delta_phi_need = max(0.0, effective_threshold - phi_predicted)

    # 勾配計算
    T2 = state.tension.T2 if state.tension else 0.0
    N  = state.N
    mu = state.mu

    grads = {
        'N':     -mu,                        # ∂Φ/∂N = -μ = -T₂
        'T2':    -N if T2 > EPS else 0.0,    # ∂Φ/∂T₂ = -N（μ=T₂ なので）
        'theta': 0.0,                        # 構造変化（数値微分は省略）
    }

    # F_geom の数値微分でθ勾配を推定
    if state.F_geom > EPS:
        grads['theta'] = state.F_geom * 0.1  # 仮の感度（Assumed）

    efficiencies = {
        k: abs(grads[k]) / (costs.get(k, 1.0) + EPS)
        for k in grads
    }

    # ④ バグ修正：λ による重み付けを分化させる
    priorities = {}
    for k in efficiencies:
        if k == 'N':
            priorities[k] = efficiencies[k]
        elif k == 'T2':
            priorities[k] = efficiencies[k] * (1.0 + 0.5 * lambda_count)
        else:
            priorities[k] = efficiencies[k] * (1.0 + lambda_count)

    best_var  = max(priorities, key=lambda k: priorities[k])
    best_sens = grads[best_var]
    delta_u   = delta_phi_need / (abs(best_sens) + EPS) if needs else 0.0
    phi_after = phi_current + best_sens * delta_u if needs else phi_current

    return InterventionStep(
        t=state.t, phi_before=phi_current,
        phi_predicted=phi_predicted, phi_threshold=phi_threshold,
        needs_intervention=needs, delta_phi_need=delta_phi_need,
        gradients=grads, efficiencies=efficiencies, priorities=priorities,
        best_variable=best_var, delta_u_optimal=delta_u, phi_after=phi_after,
        action_condition=state.action_condition,
    )


# ============================================================
# 最適化ループ
# ============================================================

@dataclass
class OptimizationResult:
    steps:               List[InterventionStep]
    total_interventions: int
    total_impact:        float
    total_cost:          float
    phi_trajectory:      List[float]
    intervention_times:  List[float]

    def summary(self) -> str:
        return "\n".join([
            "=" * 60,
            "介入最適化ループ 結果サマリー",
            "=" * 60,
            f"総ステップ数   : {len(self.steps)}",
            f"総介入回数     : {self.total_interventions}",
            f"総インパクト   : {self.total_impact:.4f}",
            f"Φ 初期値       : {self.phi_trajectory[0]:.4f}" if self.phi_trajectory else "",
            f"Φ 最終値       : {self.phi_trajectory[-1]:.4f}" if self.phi_trajectory else "",
            f"介入時刻       : {self.intervention_times}",
        ])


def run_intervention_loop(
    state_sequence: List[StateSnapshot],
    phi_threshold: float,
    tau: float = 3.0, dt: float = 1.0,
    costs: Optional[Dict[str, float]] = None,
    lambda_count: float = 1.0,
) -> OptimizationResult:
    if costs is None:
        costs = {'N': 1.0, 'T2': 1.0, 'theta': 10.0}

    steps: List[InterventionStep] = []
    phi_trajectory: List[float] = []
    intervention_times: List[float] = []
    total_impact = 0.0

    # ⑥ バグ修正：最初のステップは phi_prev = state.Phi
    phi_prev = state_sequence[0].Phi if state_sequence else 0.0

    for state in state_sequence:
        phi_trajectory.append(state.Phi)
        step = compute_intervention_step(
            state=state, phi_prev=phi_prev, dt=dt, tau=tau,
            phi_threshold=phi_threshold, costs=costs, lambda_count=lambda_count,
        )
        steps.append(step)

        if step.needs_intervention:
            intervention_times.append(state.t)
            total_impact += step.delta_u_optimal ** 2

        phi_prev = state.Phi

    total_cost = total_impact + lambda_count * len(intervention_times)
    return OptimizationResult(
        steps=steps, total_interventions=len(intervention_times),
        total_impact=total_impact, total_cost=total_cost,
        phi_trajectory=phi_trajectory, intervention_times=intervention_times,
    )


# ============================================================
# λ スキャン・最適 λ 選択
# ============================================================

@dataclass
class TradeoffPoint:
    lambda_val:          float
    total_interventions: int
    total_impact:        float
    total_cost:          float


def scan_lambda_tradeoff(
    state_sequence: List[StateSnapshot],
    phi_threshold: float,
    lambda_values: Optional[List[float]] = None,
    tau: float = 3.0, dt: float = 1.0,
    costs: Optional[Dict[str, float]] = None,
) -> List[TradeoffPoint]:
    if lambda_values is None:
        lambda_values = [10 ** x for x in np.linspace(-2, 2, 10)]
    results = []
    for lam in lambda_values:
        opt = run_intervention_loop(
            state_sequence, phi_threshold,
            tau=tau, dt=dt, costs=costs, lambda_count=lam,
        )
        results.append(TradeoffPoint(
            lambda_val=lam, total_interventions=opt.total_interventions,
            total_impact=opt.total_impact, total_cost=opt.total_cost,
        ))
    return results


def find_optimal_lambda(tradeoff_points: List[TradeoffPoint]) -> Tuple[float, Optional[TradeoffPoint]]:
    if not tradeoff_points:
        return 1.0, None
    impacts = [p.total_impact for p in tradeoff_points]
    counts  = [p.total_interventions for p in tradeoff_points]
    max_i   = max(impacts) if max(impacts) > EPS else 1.0
    max_c   = max(counts)  if max(counts)  > 0   else 1.0
    balanced = [
        (p.total_impact / max_i) + (p.total_interventions / max_c)
        for p in tradeoff_points
    ]
    best = int(np.argmin(balanced))
    return tradeoff_points[best].lambda_val, tradeoff_points[best]


# ============================================================
# Φ₀ と τ の自動設定
# ============================================================

def auto_phi_threshold(
    state: StateSnapshot,
    method: str = "ratio",
) -> Tuple[float, str]:
    """
    ⑤ バグ修正：符号不変式。
    Φ₀ = Phi - margin × (|Phi| + 1)
    """
    if method == "ratio":
        margin = 0.05
        Phi_0 = state.Phi - margin * (abs(state.Phi) + 1.0)
        desc = f"Φ から {margin*100:.0f}%·(|Φ|+1) 下方に設定"
    else:
        margin = 0.05
        Phi_0 = state.Phi - margin * (abs(state.Phi) + 1.0)
        desc = "デフォルト設定"
    return float(Phi_0), desc


def auto_tau(lambda_1: float) -> Tuple[float, str]:
    """
    予測ホライズン τ をリアプノフ指数から自動設定する。

    λ₁ < -0.1 → 安定、長期予測可能 → τ = 5
    λ₁ < 0    → やや安定           → τ = 3
    λ₁ < 0.05 → 中立              → τ = 2
    λ₁ ≥ 0.05 → カオス的          → τ = 1
    """
    if lambda_1 < -0.1:
        return 5.0, f"λ₁={lambda_1:.3f}（安定）→ τ=5"
    if lambda_1 < 0.0:
        return 3.0, f"λ₁={lambda_1:.3f}（やや安定）→ τ=3"
    if lambda_1 < 0.05:
        return 2.0, f"λ₁={lambda_1:.3f}（中立）→ τ=2"
    return 1.0, f"λ₁={lambda_1:.3f}（カオス的）→ τ=1"


# ============================================================
# 状態シーケンス生成
# ============================================================

def generate_state_sequence(
    embeddings_A: List[np.ndarray],
    embeddings_B: List[np.ndarray],
    N_series: List[float],
) -> List[StateSnapshot]:
    """埋め込み済み点群のシーケンスから StateSnapshot 列を生成する。"""
    snapshots = []
    phi_prev = None
    for t, (emb_A, emb_B, N) in enumerate(zip(embeddings_A, embeddings_B, N_series)):
        snap = StateSnapshot(t=float(t), embedded_A=emb_A, embedded_B=emb_B, N=N)
        snap.compute(phi_prev=phi_prev, dt=1.0)
        phi_prev = snap.Phi
        snapshots.append(snap)
    return snapshots
