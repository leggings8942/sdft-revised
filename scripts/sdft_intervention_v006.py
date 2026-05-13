"""
sdft_intervention_v006.py
SDFT revised v0.0.6 — 介入最適化ループと行動条件診断

このモジュールは SDFT の Layer 3 を提供する：
  ① StateSnapshot：時刻 t における全 SDFT 変数のスナップショット
  ② 行動条件診断（§4-2）：既存変数の読み替えによる6項目評価
  ③ 介入最適化ループ：Φ ≥ Φ₀ を維持する最小インパクト・最小回数の介入計画

【行動条件診断の6項目（既存変数の読み替えのみ）】
  T₂    → 到達可能性（移動コスト μ = 2√T₂ の低さ）
  T₅    → 伝達効率（情報損失の低さ）
  T₃    → 方向整合性（非対称差の小ささ）
  T₁    → スケール適合（レジーム間のスケール差の小ささ）
  dΦ/dt → 変化の好機（Φ が減少中＝変化の機が熟している）
  H     → 文脈の持続性（Hurst 指数 > 0.5 でトレンドが維持）

【介入最適化の定式化】
  min  Σ‖u_k‖² + λK
  s.t. Φ(t) ≥ Φ₀  ∀t

  操作変数：
    u₁ = N  （人数）       ∂Φ/∂N  = -2√T₂
    u₂ = T₂ （距離）       ∂Φ/∂T₂ = -N/√T₂
    u₃ = θ  （構造）       ∂Φ/∂θ  ≈ 数値微分

  アルゴリズム：
    Step 1: dΦ/dt を推定
    Step 2: Φ̂(t+τ) を予測
    Step 3: Φ̂(t+τ) < Φ₀ なら介入
    Step 4: priority_i = |∂Φ/∂u_i| / c_i を計算
    Step 5: 最高効率の操作変数 u* と最小操作量 Δu* を決定

設計原則：
  1. アフォーダンス概念を使わない
  2. 新しい変数を一切追加しない
  3. 全項目を Derived として計算する
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Optional
import numpy as np

from sdft_v006_core import (
    K_B, LN2, EPS,
    Evidence,
    calc_H_bit, calc_S, calc_D, calc_H,
    calc_fisher_matrix, calc_T_temperature, calc_U_internal_energy,
    calc_F_free_energy, calc_V_volume,
    calc_regime_tension_vector,
    calc_mu, calc_grand_potential,
    classify_phase, detect_phase_transition,
    RegimeTensionVector,
)


# ============================================================
# 行動条件診断
# ============================================================

@dataclass
class ActionConditionDiagnosis:
    """
    §4-2 行動条件診断。
    新しい変数を作らず、既存の𝓣・Φ・H の読み替えで構成する。

    6つの診断項目それぞれについて：
      変数値・ステータス（○/△/×）・解釈文を返す。
    """

    # 各診断項目の値（既存変数から直接取得）
    T1:      float  # 対数分配関数差（スケール適合）
    T2:      float  # Bhattacharyya距離（到達可能性）
    T3:      float  # α-ダイバージェンス（方向整合性）
    T5:      float  # 通信路損失率（伝達効率）
    mu:      float  # 化学ポテンシャル = 2√T₂（移動コスト）
    H:       float  # Hurst指数（文脈の持続性）
    dPhi_dt: float  # dΦ/dt（変化の好機）

    # 各診断項目のステータス（Derived）
    status_reachability:  str = ""  # 到達可能性
    status_transmission:  str = ""  # 伝達効率
    status_alignment:     str = ""  # 方向整合性
    status_scale:         str = ""  # スケール適合
    status_timing:        str = ""  # 変化の好機
    status_persistence:   str = ""  # 文脈の持続性

    # 解釈文
    interp_reachability:  str = ""
    interp_transmission:  str = ""
    interp_alignment:     str = ""
    interp_scale:         str = ""
    interp_timing:        str = ""
    interp_persistence:   str = ""

    # 総合判定
    overall_status:       str = ""
    overall_message:      str = ""

    def __post_init__(self):
        self._diagnose()

    def _diagnose(self):
        # ① 到達可能性（T₂・μ）
        if self.mu < 0.5:
            self.status_reachability = "good"
            self.interp_reachability = f"移動コスト μ={self.mu:.3f} — 低い。顧客は動きやすい状態。"
        elif self.mu < 1.5:
            self.status_reachability = "warn"
            self.interp_reachability = f"移動コスト μ={self.mu:.3f} — 中程度。一定の障壁がある。"
        else:
            self.status_reachability = "alert"
            self.interp_reachability = f"移動コスト μ={self.mu:.3f} — 高い。レジーム間の距離が大きく行動しにくい。"

        # ② 伝達効率（T₅）
        if self.T5 < 0.3:
            self.status_transmission = "good"
            self.interp_transmission = f"通信路損失 T₅={self.T5:.3f} — 低い。情報が届いており気づきやすい。"
        elif self.T5 < 0.6:
            self.status_transmission = "warn"
            self.interp_transmission = f"通信路損失 T₅={self.T5:.3f} — 中程度。一部の顧客に届いていない可能性。"
        else:
            self.status_transmission = "alert"
            self.interp_transmission = f"通信路損失 T₅={self.T5:.3f} — 高い。情報が大幅に損失しており、顧客が気づきにくい。"

        # ③ 方向整合性（T₃）
        if abs(self.T3) < 0.3:
            self.status_alignment = "good"
            self.interp_alignment = f"方向差 T₃={self.T3:.3f} — 小さい。双方向に自然な移動。"
        elif abs(self.T3) < 1.0:
            self.status_alignment = "warn"
            self.interp_alignment = f"方向差 T₃={self.T3:.3f} — 中程度。一方向への傾きがある。"
        else:
            self.status_alignment = "alert"
            self.interp_alignment = f"方向差 T₃={self.T3:.3f} — 大きい。非対称な差が強く、片方向の移動に偏っている。"

        # ④ スケール適合（T₁）
        if self.T1 < 0.3:
            self.status_scale = "good"
            self.interp_scale = f"スケール差 T₁={self.T1:.3f} — 小さい。2レジームの文脈が近い。"
        elif self.T1 < 1.0:
            self.status_scale = "warn"
            self.interp_scale = f"スケール差 T₁={self.T1:.3f} — 中程度。スケールのズレが一部行動を阻害。"
        else:
            self.status_scale = "alert"
            self.interp_scale = f"スケール差 T₁={self.T1:.3f} — 大きい。レジーム間のスケール差が行動の文脈を壊している。"

        # ⑤ 変化の好機（dΦ/dt）
        if self.dPhi_dt < -0.01:
            self.status_timing = "good"
            self.interp_timing = f"dΦ/dt={self.dPhi_dt:.4f} — Φ が減少中。変化の機が熟している。"
        elif self.dPhi_dt < 0.0:
            self.status_timing = "warn"
            self.interp_timing = f"dΦ/dt={self.dPhi_dt:.4f} — わずかに減少。好機の兆し。"
        else:
            self.status_timing = "mid"
            self.interp_timing = f"dΦ/dt={self.dPhi_dt:.4f} — Φ が安定または上昇中。変化の動機が弱い。"

        # ⑥ 文脈の持続性（H）
        if self.H >= 0.55:
            self.status_persistence = "good"
            self.interp_persistence = f"H={self.H:.3f} — 持続傾向（H>0.5）。行動の文脈が維持されやすい。"
        elif self.H >= 0.45:
            self.status_persistence = "warn"
            self.interp_persistence = f"H={self.H:.3f} — ランダムに近い。行動の文脈が不安定。"
        else:
            self.status_persistence = "alert"
            self.interp_persistence = f"H={self.H:.3f} — 反持続傾向（H<0.5）。行動が起きても継続しにくい。"

        # 総合判定
        statuses = [
            self.status_reachability,
            self.status_transmission,
            self.status_alignment,
            self.status_scale,
            self.status_timing,
            self.status_persistence,
        ]
        alert_count = statuses.count("alert")
        good_count  = statuses.count("good")

        if alert_count == 0 and good_count >= 4:
            self.overall_status  = "good"
            self.overall_message = "行動が起きやすい状態。全指標が良好または許容範囲内。"
        elif alert_count >= 3:
            self.overall_status  = "alert"
            self.overall_message = "行動が起きにくい状態。複数の障壁が同時に存在している。"
        elif alert_count >= 1:
            self.overall_status  = "warn"
            # 最も深刻な項目を特定
            alert_items = []
            if self.status_reachability == "alert":
                alert_items.append("移動コストが高い（T₂）")
            if self.status_transmission == "alert":
                alert_items.append("情報損失が大きい（T₅）")
            if self.status_alignment == "alert":
                alert_items.append("方向の非対称差が強い（T₃）")
            if self.status_scale == "alert":
                alert_items.append("スケール差が大きい（T₁）")
            self.overall_message = "一部に障壁あり。優先課題：" + "・".join(alert_items)
        else:
            self.overall_status  = "mid"
            self.overall_message = "行動条件は中程度。特定の介入で改善の余地がある。"

    def table_rows(self) -> List[Dict]:
        """レポートテンプレート用のテーブル行データを返す。"""
        return [
            {
                "item": "到達可能性",
                "variable": f"μ = 2√T₂ = {self.mu:.4f}",
                "status": self.status_reachability,
                "interpretation": self.interp_reachability,
            },
            {
                "item": "伝達効率",
                "variable": f"T₅ = {self.T5:.4f}",
                "status": self.status_transmission,
                "interpretation": self.interp_transmission,
            },
            {
                "item": "方向整合性",
                "variable": f"T₃ = {self.T3:.4f}",
                "status": self.status_alignment,
                "interpretation": self.interp_alignment,
            },
            {
                "item": "スケール適合",
                "variable": f"T₁ = {self.T1:.4f}",
                "status": self.status_scale,
                "interpretation": self.interp_scale,
            },
            {
                "item": "変化の好機",
                "variable": f"dΦ/dt = {self.dPhi_dt:.4f}",
                "status": self.status_timing,
                "interpretation": self.interp_timing,
            },
            {
                "item": "文脈の持続性",
                "variable": f"H = {self.H:.4f}",
                "status": self.status_persistence,
                "interpretation": self.interp_persistence,
            },
        ]

    def summary(self) -> str:
        lines = [
            "=== 行動条件診断 ===",
            f"総合：{self.overall_status.upper()}  {self.overall_message}",
            "",
        ]
        for row in self.table_rows():
            mark = {"good": "○", "warn": "△", "alert": "×", "mid": "－"}.get(
                row["status"], "－"
            )
            lines.append(f"  {mark} {row['item']:10s}  {row['variable']:20s}")
            lines.append(f"    {row['interpretation']}")
        return "\n".join(lines)


def calc_action_condition(
    T1:      float,
    T2:      float,
    T3:      float,
    T5:      float,
    H:       float,
    phi_current: float,
    phi_prev:    float,
    dt:          float = 1.0,
) -> ActionConditionDiagnosis:
    """
    行動条件診断を計算する。

    Parameters
    ----------
    T1, T2, T3, T5 : 𝓣 ベクトルの各成分
    H              : Hurst 指数
    phi_current    : 現在の Φ
    phi_prev       : 前時点の Φ
    dt             : 時間ステップ

    Returns
    -------
    ActionConditionDiagnosis
    """
    mu      = calc_mu(T2)
    dPhi_dt = (phi_current - phi_prev) / (dt + EPS)

    return ActionConditionDiagnosis(
        T1=T1, T2=T2, T3=T3, T5=T5,
        mu=mu, H=H, dPhi_dt=dPhi_dt,
    )


# ============================================================
# StateSnapshot
# ============================================================

@dataclass
class StateSnapshot:
    """
    時刻 t における SDFT 変数の完全なスナップショット。
    全ての構造場変数・𝓣・Φ・行動条件診断を一つのオブジェクトに保持する。
    """
    t:    float
    x:    List[float]   # レジームAの時系列データ
    x_B:  List[float]   # レジームBの時系列データ

    # 構造場変数
    H_bit: float = 0.0
    S:     float = 0.0
    D:     float = 0.0
    H:     float = 0.5

    # Fisher 行列由来
    G:    Optional[np.ndarray] = None
    T:    float = 0.0
    U:    float = 0.0
    F:    float = 0.0
    V:    float = 0.0

    # 𝓣 ベクトル
    tension: Optional[RegimeTensionVector] = None

    # グランドポテンシャル
    N:   float = 0.0
    mu:  float = 0.0
    Phi: float = 0.0

    # 行動条件診断（計算後に設定）
    action_condition: Optional[ActionConditionDiagnosis] = None

    def compute(
        self,
        phi_prev: float = 0.0,
        dt: float = 1.0,
        alpha_T3: float = 1.0,
    ) -> "StateSnapshot":
        """全変数を計算して自己を返す。"""
        # 構造場変数
        self.H_bit = calc_H_bit(self.x)
        self.S     = calc_S(self.x)
        self.D     = calc_D(self.x)
        self.H     = calc_H(self.x)

        # Fisher 行列由来
        G = calc_fisher_matrix(self.x)
        self.G = G
        self.T = calc_T_temperature(G)
        self.U = calc_U_internal_energy(G)
        self.F = calc_F_free_energy(G, self.H_bit)
        self.V = calc_V_volume(G)

        # 𝓣 ベクトル
        self.tension = calc_regime_tension_vector(
            self.x, self.x_B, alpha_T3=alpha_T3
        )

        # グランドポテンシャル
        self.mu  = calc_mu(self.tension.T2)
        self.Phi = calc_grand_potential(self.F, self.mu, self.N)

        # 行動条件診断（アフォーダンスの代替）
        self.action_condition = calc_action_condition(
            T1=self.tension.T1,
            T2=self.tension.T2,
            T3=self.tension.T3,
            T5=self.tension.T5,
            H=self.H,
            phi_current=self.Phi,
            phi_prev=phi_prev,
            dt=dt,
        )
        return self


# ============================================================
# 介入ステップ（アフォーダンス → 行動条件診断に置き換え）
# ============================================================

@dataclass
class InterventionStep:
    """1ステップの介入計算結果。"""
    t:                  float
    phi_before:         float
    phi_predicted:      float
    phi_threshold:      float
    needs_intervention: bool
    delta_phi_need:     float

    gradients:       Dict[str, float] = field(default_factory=dict)
    efficiencies:    Dict[str, float] = field(default_factory=dict)
    priorities:      Dict[str, float] = field(default_factory=dict)

    best_variable:   str   = ""
    delta_u_optimal: float = 0.0
    phi_after:       float = 0.0

    # 行動条件診断（アフォーダンスの代替）
    action_condition: Optional[ActionConditionDiagnosis] = None

    def summary(self) -> str:
        lines = [
            f"t={self.t:.1f}  Φ={self.phi_before:.4f}  "
            f"Φ̂={self.phi_predicted:.4f}  Φ₀={self.phi_threshold:.4f}",
            f"  介入: {'必要' if self.needs_intervention else '不要'}  "
            f"ΔΦ_need={self.delta_phi_need:.4f}",
        ]
        if self.needs_intervention:
            lines.append(
                f"  最優先: {self.best_variable}  "
                f"Δu*={self.delta_u_optimal:.4f}  "
                f"Φ_after={self.phi_after:.4f}"
            )
            lines.append(
                "  priority: "
                + "  ".join(f"{k}={v:.4f}" for k, v in self.priorities.items())
            )
        if self.action_condition:
            lines.append(
                f"  行動条件: {self.action_condition.overall_status.upper()}  "
                f"{self.action_condition.overall_message}"
            )
        return "\n".join(lines)


# ============================================================
# 勾配計算・介入ステップ計算
# ============================================================

def numerical_grad_Phi(
    state: StateSnapshot,
    epsilon: float = 1e-4,
) -> Dict[str, float]:
    """Φ の各操作変数に対する勾配を計算する。"""
    T2 = state.tension.T2 if state.tension else 0.0
    N  = state.N
    mu = calc_mu(T2)

    grad_N  = -mu
    grad_T2 = -N / (math.sqrt(max(T2, EPS)))

    if state.G is not None:
        G = state.G
        d = G.shape[0]
        try:
            G_inv          = np.linalg.inv(G + np.eye(d) * EPS)
            G_plus         = G * (1 + epsilon)
            G_minus        = G * (1 - epsilon)
            F_plus         = calc_F_free_energy(G_plus,  state.H_bit)
            F_minus        = calc_F_free_energy(G_minus, state.H_bit)
            Phi_plus       = calc_grand_potential(F_plus,  mu, N)
            Phi_minus      = calc_grand_potential(F_minus, mu, N)
            grad_theta     = (Phi_plus - Phi_minus) / (2 * epsilon)
        except Exception:
            grad_theta = 0.0
    else:
        grad_theta = 0.0

    return {
        'N':     float(grad_N),
        'T2':    float(grad_T2),
        'theta': float(grad_theta),
    }


def compute_intervention_step(
    state:          StateSnapshot,
    phi_prev:       float,
    dt:             float,
    tau:            float,
    phi_threshold:  float,
    costs:          Dict[str, float],
    lambda_count:   float = 1.0,
) -> InterventionStep:
    """1ステップの介入計算。"""
    phi_current   = state.Phi
    dPhi_dt       = (phi_current - phi_prev) / (dt + EPS)
    phi_predicted = phi_current + dPhi_dt * tau
    needs         = phi_predicted < phi_threshold
    delta_phi_need = max(0.0, phi_threshold - phi_predicted)

    grads = numerical_grad_Phi(state)
    efficiencies = {
        k: abs(grads[k]) / (costs.get(k, 1.0) + EPS)
        for k in grads
    }
    priorities = {
        k: efficiencies[k] * (1.0 + lambda_count)
        for k in efficiencies
    }

    best_var  = max(priorities, key=lambda k: priorities[k])
    best_sens = grads[best_var]
    delta_u   = delta_phi_need / best_sens if (needs and abs(best_sens) > EPS) else 0.0
    phi_after = phi_current + best_sens * delta_u if needs else phi_current

    return InterventionStep(
        t=state.t,
        phi_before=phi_current,
        phi_predicted=phi_predicted,
        phi_threshold=phi_threshold,
        needs_intervention=needs,
        delta_phi_need=delta_phi_need,
        gradients=grads,
        efficiencies=efficiencies,
        priorities=priorities,
        best_variable=best_var,
        delta_u_optimal=delta_u,
        phi_after=phi_after,
        action_condition=state.action_condition,
    )


# ============================================================
# 最適化ループ
# ============================================================

@dataclass
class OptimizationResult:
    """介入最適化ループ全体の結果。"""
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
            f"総コスト       : {self.total_cost:.4f}",
            f"Φ 初期値       : {self.phi_trajectory[0]:.4f}",
            f"Φ 最終値       : {self.phi_trajectory[-1]:.4f}",
            f"介入時刻       : {self.intervention_times}",
        ])


def run_intervention_loop(
    state_sequence: List[StateSnapshot],
    phi_threshold:  float,
    tau:            float = 3.0,
    dt:             float = 1.0,
    costs:          Dict[str, float] = None,
    lambda_count:   float = 1.0,
    max_iter:       int   = 100,
    convergence_tol: float = 1e-6,
) -> OptimizationResult:
    """介入最適化ループを実行する。"""
    if costs is None:
        costs = {'N': 1.0, 'T2': 1.0, 'theta': 10.0}

    steps:              List[InterventionStep] = []
    phi_trajectory:     List[float] = []
    intervention_times: List[float] = []
    total_impact:       float = 0.0
    phi_prev = state_sequence[0].Phi if state_sequence else 0.0

    for state in state_sequence:
        phi_trajectory.append(state.Phi)
        step = compute_intervention_step(
            state=state,
            phi_prev=phi_prev,
            dt=dt,
            tau=tau,
            phi_threshold=phi_threshold,
            costs=costs,
            lambda_count=lambda_count,
        )
        steps.append(step)

        if step.needs_intervention:
            intervention_times.append(state.t)
            total_impact += step.delta_u_optimal ** 2
            phi_est   = step.phi_after
            best_sens = step.gradients.get(step.best_variable, 0.0)
            for _ in range(max_iter):
                if phi_est >= phi_threshold:
                    break
                extra = (phi_threshold - phi_est) / (best_sens + EPS)
                phi_est    += best_sens * extra
                total_impact += extra ** 2
                if abs(extra) < convergence_tol:
                    break
            step.phi_after = phi_est

        phi_prev = state.Phi

    total_cost = total_impact + lambda_count * len(intervention_times)
    return OptimizationResult(
        steps=steps,
        total_interventions=len(intervention_times),
        total_impact=total_impact,
        total_cost=total_cost,
        phi_trajectory=phi_trajectory,
        intervention_times=intervention_times,
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
    phi_threshold:  float,
    lambda_values:  List[float] = None,
    tau:   float = 3.0,
    dt:    float = 1.0,
    costs: Dict[str, float] = None,
) -> List[TradeoffPoint]:
    if lambda_values is None:
        lambda_values = [10 ** x for x in np.linspace(-2, 2, 10)]
    results = []
    for lam in lambda_values:
        opt = run_intervention_loop(
            state_sequence=state_sequence,
            phi_threshold=phi_threshold,
            tau=tau, dt=dt, costs=costs, lambda_count=lam,
        )
        results.append(TradeoffPoint(
            lambda_val=lam,
            total_interventions=opt.total_interventions,
            total_impact=opt.total_impact,
            total_cost=opt.total_cost,
        ))
    return results


def find_optimal_lambda(
    tradeoff_points: List[TradeoffPoint],
) -> Tuple[float, Optional[TradeoffPoint]]:
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
# 状態シーケンス生成
# ============================================================

def generate_state_sequence(
    x_series:   List[List[float]],
    x_B_series: List[List[float]],
    N_series:   List[float],
    alpha_T3:   float = 1.0,
) -> List[StateSnapshot]:
    snapshots = []
    phi_prev  = 0.0
    for t, (x, x_B, N) in enumerate(zip(x_series, x_B_series, N_series)):
        snap = StateSnapshot(t=float(t), x=x, x_B=x_B, N=N)
        snap.compute(phi_prev=phi_prev, dt=1.0, alpha_T3=alpha_T3)
        phi_prev = snap.Phi
        snapshots.append(snap)
    return snapshots


# ============================================================
# メイン実行例
# ============================================================

if __name__ == "__main__":
    import random
    random.seed(42)

    print("=" * 60)
    print("SDFT revised v0.0.6 — 介入最適化・行動条件診断 動作確認")
    print("=" * 60)

    T_steps  = 15
    win_size = 200
    N_people = 100.0

    x_series, x_B_series, N_series = [], [], []
    for t in range(T_steps):
        noise = 50 + t * 8
        x_series.append([500 + random.gauss(0, noise) + i * 0.1 for i in range(win_size)])
        x_B_series.append([350 + random.gauss(0, noise * 1.2) for _ in range(win_size)])
        N_series.append(N_people)

    states = generate_state_sequence(x_series, x_B_series, N_series)
    print(f"\n生成した状態数: {len(states)}")

    print("\n--- 行動条件診断（t=0）---")
    print(states[0].action_condition.summary())

    print("\n--- λスキャン ---")
    phi_0    = states[0].Phi * 0.95
    costs    = {'N': 1.0, 'T2': 2.0, 'theta': 20.0}
    tradeoff = scan_lambda_tradeoff(
        states, phi_0,
        lambda_values=[0.01, 0.1, 0.5, 1.0, 5.0, 10.0],
        costs=costs,
    )
    print(f"{'λ':>8}  {'介入回数':>8}  {'インパクト':>12}  {'総コスト':>12}")
    print("-" * 50)
    for p in tradeoff:
        print(f"{p.lambda_val:>8.3f}  {p.total_interventions:>8d}  "
              f"{p.total_impact:>12.4f}  {p.total_cost:>12.4f}")

    opt_lambda, opt_point = find_optimal_lambda(tradeoff)
    print(f"\n最適 λ = {opt_lambda:.3f}  "
          f"介入回数: {opt_point.total_interventions}  "
          f"インパクト: {opt_point.total_impact:.4f}")

    result = run_intervention_loop(
        states, phi_0, tau=3.0, dt=1.0,
        costs=costs, lambda_count=opt_lambda,
    )
    print("\n" + result.summary())
    print("=" * 60)
    print("動作確認完了")
    print("=" * 60)
