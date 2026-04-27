"""
SDFT Unified v1.4-A Reference Implementation
Merged package of:
- sdft_unified_v1_4_reference.py
- sdft_affordance_extension_v0_1.py

Gauge-Invariant Semantic Field / Maxwellian Propagation / Hamiltonian Intervention
with Affordance Layer / Cognitive-Normative Mediation / Behavior-Oriented Intervention

This module preserves the original S, D, H, F, 𝓣 backbone of SDFT Unified v1.4.1
and extends it by adding A = Affordance as a second-layer field that mediates
between structural state and behavior.

This is a practical integrated reference implementation, not a proof of
mathematical completeness.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

try:
    import networkx as nx
except Exception:  # pragma: no cover
    nx = None

EPS = 1e-12


# ==========================================================
# Phase 1–4: Base utilities carried over from v1.4.1
# ==========================================================

def calc_entropy(data: Sequence[float], bins: int = 32) -> float:
    x = np.asarray(data, dtype=float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return 0.0
    counts, _ = np.histogram(x, bins=bins, density=False)
    p = counts[counts > 0] / np.sum(counts)
    return float(-np.sum(p * np.log(p + EPS)))


def higuchi_fd(data: Sequence[float], k_max: int = 10) -> float:
    x = np.asarray(data, dtype=float)
    n = len(x)
    if n < 8 or np.std(x) < EPS:
        return 1.0
    lk = []
    k_range = range(1, min(k_max, max(2, n // 2)) + 1)
    for k in k_range:
        lm = []
        for m in range(k):
            idx = np.arange(m, n, k)
            if len(idx) < 2:
                continue
            diff = np.abs(np.diff(x[idx])).sum()
            denom = ((n - m - 1) // k) * k
            if denom <= 0:
                continue
            norm = (n - 1) / denom
            lm.append((diff * norm) / k)
        if lm:
            lk.append((k, np.mean(lm)))
    if len(lk) < 2:
        return 1.0
    logs_k = np.log([1.0 / k for k, _ in lk])
    logs_l = np.log([l + EPS for _, l in lk])
    slope, _ = np.polyfit(logs_k, logs_l, 1)
    return float(max(1.0, min(2.5, slope)))


def hurst_exponent(data: Sequence[float], min_lag: int = 2, max_lag: int = 20) -> float:
    x = np.asarray(data, dtype=float)
    if len(x) < max_lag + 2 or np.std(x) < EPS:
        return 0.5
    lags = np.arange(min_lag, min(max_lag, len(x) // 2))
    if len(lags) < 2:
        return 0.5
    tau = []
    for lag in lags:
        diff = x[lag:] - x[:-lag]
        tau.append(np.sqrt(np.std(diff) + EPS))
    slope, _ = np.polyfit(np.log(lags), np.log(np.asarray(tau) + EPS), 1)
    return float(np.clip(2.0 * slope, 0.0, 1.5))


def sliding_windows(x: Sequence[float], window: int, step: int = 1):
    arr = np.asarray(x, dtype=float)
    for i in range(window, len(arr) + 1, step):
        yield i, arr[i - window:i]


def analyze_sdh_series(x: Sequence[float], window: int = 128, step: int = 8) -> pd.DataFrame:
    rows = []
    for t, w in sliding_windows(x, window=window, step=step):
        rows.append({
            't': int(t),
            'S': calc_entropy(w),
            'D': higuchi_fd(w),
            'H': hurst_exponent(w),
        })
    return pd.DataFrame(rows)


def phase_alert(df: pd.DataFrame, h_drop: float = -0.03, d_var_q: float = 0.75,
                s_rise: float = 0.0, lookback: int = 5) -> pd.DataFrame:
    out = df.copy()
    out['dH'] = out['H'].diff()
    out['dS'] = out['S'].diff()
    out['D_var'] = out['D'].rolling(lookback).var()
    thr = out['D_var'].quantile(d_var_q) if len(out) else np.nan
    out['alert'] = (out['dH'] < h_drop) & (out['D_var'] > thr) & (out['dS'] > s_rise)
    return out


# ==========================================================
# Quaternion layer
# ==========================================================

def q(w: float, x: float, y: float, z: float) -> np.ndarray:
    return np.asarray([w, x, y, z], dtype=float)


Q_ID = q(1.0, 0.0, 0.0, 0.0)


def q_norm(a: np.ndarray) -> float:
    return float(np.linalg.norm(a))


def q_unit(a: np.ndarray) -> np.ndarray:
    n = q_norm(a)
    if n < EPS:
        return Q_ID.copy()
    return np.asarray(a, dtype=float) / n


def q_conj(a: np.ndarray) -> np.ndarray:
    a = np.asarray(a, dtype=float)
    return np.asarray([a[0], -a[1], -a[2], -a[3]], dtype=float)


def q_inv(a: np.ndarray) -> np.ndarray:
    a = np.asarray(a, dtype=float)
    n2 = np.dot(a, a)
    if n2 < EPS:
        return Q_ID.copy()
    return q_conj(a) / n2


def q_mul(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    aw, ax, ay, az = a
    bw, bx, by, bz = b
    return np.asarray([
        aw * bw - ax * bx - ay * by - az * bz,
        aw * bx + ax * bw + ay * bz - az * by,
        aw * by - ax * bz + ay * bw + az * bx,
        aw * bz + ax * by - ay * bx + az * bw,
    ], dtype=float)


def q_distance(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.linalg.norm(np.asarray(a, dtype=float) - np.asarray(b, dtype=float)))


def q_log_unit(a: np.ndarray) -> np.ndarray:
    """Logarithm of a unit quaternion; returns a pure-imaginary quaternion."""
    u = q_unit(a)
    w = float(np.clip(u[0], -1.0, 1.0))
    v = u[1:]
    nv = np.linalg.norm(v)
    if nv < EPS:
        return np.asarray([0.0, 0.0, 0.0, 0.0], dtype=float)
    theta = math.acos(w)
    axis = v / nv
    return np.concatenate([[0.0], axis * theta])


def gauge_unit(g: np.ndarray) -> np.ndarray:
    return q_unit(g)


# ==========================================================
# Thermodynamic / regime layer
# ==========================================================

def free_energy(U: float, Theta: float, S: float) -> float:
    """Helmholtz-like free energy. Theta is temperature-like, not regime tension."""
    return float(U - Theta * S)


def regime_tension(LA: float, LB: float, delta_topo: float, delta_op: float,
                   kappa: float = 1.0, lam: float = 1.0,
                   gauge_dist: float = 0.0, eta: float = 1.0) -> float:
    ratio = max(float(LA), EPS) / max(float(LB), EPS)
    return float(abs(math.log(ratio)) + kappa * delta_topo + lam * delta_op + eta * gauge_dist)


# ==========================================================
# Gauge-invariant semantic field
# ==========================================================

def gauge_transform_node(q_i: np.ndarray, g_i: np.ndarray) -> np.ndarray:
    return q_mul(gauge_unit(g_i), np.asarray(q_i, dtype=float))


def gauge_transform_edge(W_ij: np.ndarray, g_i: np.ndarray, g_j: np.ndarray) -> np.ndarray:
    """W'_ij = g_j ⊗ W_ij ⊗ g_i^{-1}"""
    return q_mul(q_mul(gauge_unit(g_j), np.asarray(W_ij, dtype=float)), q_inv(gauge_unit(g_i)))


def apply_edge(W_ij: np.ndarray, q_i: np.ndarray) -> np.ndarray:
    return q_mul(np.asarray(W_ij, dtype=float), np.asarray(q_i, dtype=float))


def semantic_residual(q_j: np.ndarray, W_ij: np.ndarray, q_i: np.ndarray) -> float:
    return q_distance(np.asarray(q_j, dtype=float), apply_edge(W_ij, q_i))


def path_holonomy(W_path: Sequence[np.ndarray]) -> np.ndarray:
    acc = Q_ID.copy()
    for W in W_path:
        acc = q_mul(np.asarray(W, dtype=float), acc)
    return q_unit(acc)


def gauge_distortion(W_cycle: Sequence[np.ndarray]) -> float:
    hol = path_holonomy(W_cycle)
    return q_distance(hol, Q_ID)


def action_value(node_states: Mapping[str, np.ndarray], edge_ops: Mapping[Tuple[str, str], np.ndarray],
                 edges: Sequence[Tuple[str, str]], grad_phi_penalty: float = 0.0,
                 node_phi: Optional[Mapping[str, float]] = None,
                 loop_cycles: Optional[Sequence[Sequence[Tuple[str, str]]]] = None,
                 beta_hol: float = 1.0) -> float:
    total = 0.0
    for i, j in edges:
        total += semantic_residual(node_states[j], edge_ops[(i, j)], node_states[i]) ** 2
    if node_phi is not None and grad_phi_penalty > 0:
        for i, j in edges:
            total += grad_phi_penalty * (float(node_phi[j]) - float(node_phi[i])) ** 2
    if loop_cycles is not None:
        for cyc in loop_cycles:
            Wc = [edge_ops[e] for e in cyc]
            total += beta_hol * gauge_distortion(Wc) ** 2
    return float(total)


# ==========================================================
# Affordance layer (new in SDFT-A)
# ==========================================================

def sigmoid(x: float) -> float:
    x = float(np.clip(x, -60.0, 60.0))
    return float(1.0 / (1.0 + math.exp(-x)))


def affordance_logit(salience: float, feasibility: float, legitimacy: float,
                     reward: float, mutuality: float, cognitive_cost: float,
                     w_salience: float = 1.0, w_feasibility: float = 1.0,
                     w_legitimacy: float = 1.0, w_reward: float = 1.0,
                     w_mutuality: float = 1.0, w_cost: float = 1.0,
                     bias: float = 0.0) -> float:
    """Linear logit before the sigmoid transformation."""
    return float(
        bias
        + w_salience * salience
        + w_feasibility * feasibility
        + w_legitimacy * legitimacy
        + w_reward * reward
        + w_mutuality * mutuality
        - w_cost * cognitive_cost
    )


def affordance_score(salience: float, feasibility: float, legitimacy: float,
                     reward: float, mutuality: float, cognitive_cost: float,
                     **weights: float) -> float:
    return sigmoid(affordance_logit(
        salience, feasibility, legitimacy, reward, mutuality, cognitive_cost,
        **weights,
    ))


def affordance_potential(A: float, nu_A: float = 1.0) -> float:
    return float(nu_A * A)


def total_potential(phi_struct: float, phi_aff: float, gamma_aff: float = 1.0) -> float:
    return float(phi_struct + gamma_aff * phi_aff)


def behavior_probability(A: float, F: float, H: float, T_reg: float,
                         alpha_A: float = 2.0, alpha_F: float = 1.0,
                         alpha_H: float = 0.5, alpha_T: float = 1.25,
                         bias: float = 0.0) -> float:
    return sigmoid(bias + alpha_A * A + alpha_F * F + alpha_H * H - alpha_T * T_reg)


def edge_affordance_flow(A_src: float, affordance_gain: float = 0.0,
                         normative_barrier: float = 0.0, cognitive_barrier: float = 0.0,
                         semantic_penalty: float = 0.0, length: float = 1.0) -> float:
    """
    Expected transmitted action-possibility across an edge.
    Higher barriers and semantic mismatch reduce the transmitted affordance.
    """
    loss = normative_barrier + cognitive_barrier + semantic_penalty + max(length - 1.0, 0.0) * 0.1
    return float(np.clip(A_src + affordance_gain - loss, 0.0, 1.0))


def cognitive_barrier_proxy(cognitive_cost: float, lack_of_legitimacy: float,
                            uncertainty: float = 0.0) -> float:
    return float(max(0.0, cognitive_cost) + max(0.0, lack_of_legitimacy) + max(0.0, uncertainty))


# ==========================================================
# Maxwellian propagation layer
# ==========================================================

def scalar_potential(F: float, nu_F: float = 1.0) -> float:
    return float(nu_F * F)


def vector_potential(W_ij: np.ndarray) -> np.ndarray:
    return q_log_unit(W_ij)[1:]


def electric_field(phi_i: float, phi_j: float, dA_dt: float = 0.0, length: float = 1.0) -> float:
    ell = max(float(length), EPS)
    return float(-((phi_j - phi_i) / ell) - dA_dt)


def magnetic_field(W_cycle: Sequence[np.ndarray], loop_scale: float = 1.0) -> np.ndarray:
    hol = path_holonomy(W_cycle)
    return q_log_unit(hol)[1:] / max(float(loop_scale), EPS)


def material_params(S: float, T_reg: float,
                    eps0: float = 1.0, mu0: float = 1.0, sigma0: float = 0.05,
                    a_s: float = 1.0, b_t: float = 1.0) -> Tuple[float, float, float]:
    """
    Heuristic medium parameters.
    epsilon: receptivity / slack
    mu: inertia / institutional mass
    sigma: dissipation (entropy + friction)
    """
    S = max(float(S), 0.0)
    T_reg = max(float(T_reg), 0.0)
    eps = eps0 * (1.0 + 1.0 / (1.0 + T_reg + EPS))
    mu = mu0 * (1.0 + T_reg)
    sigma = sigma0 + a_s * S + b_t * T_reg
    return float(eps), float(mu), float(sigma)


def propagation_speed(eps: float, mu: float) -> float:
    return float(1.0 / math.sqrt(max(eps * mu, EPS)))


def impedance(eps: float, mu: float) -> float:
    return float(math.sqrt(max(mu / max(eps, EPS), EPS)))


def attenuation_alpha(eps: float, mu: float, sigma: float) -> float:
    return float(0.5 * sigma * math.sqrt(max(mu / max(eps, EPS), EPS)))


def transmission_coefficient(Z1: float, Z2: float) -> float:
    return float((4.0 * Z1 * Z2) / ((Z1 + Z2) ** 2 + EPS))


def arrival_time(path_lengths: Sequence[float], path_speeds: Sequence[float]) -> float:
    total = 0.0
    for ell, v in zip(path_lengths, path_speeds):
        total += float(ell) / max(float(v), EPS)
    return float(total)


def propagated_amplitude(A0: float, path_lengths: Sequence[float], alphas: Sequence[float]) -> float:
    loss = 0.0
    for ell, a in zip(path_lengths, alphas):
        loss += float(ell) * float(a)
    return float(A0 * math.exp(-loss))


# ==========================================================
# Hamiltonian / charged-node layer
# ==========================================================

def hamiltonian(p: Sequence[float], A: Sequence[float], charge: float, phi: float, mass: float,
                D: float = 0.0, T_reg: float = 0.0,
                alpha_D: float = 0.0, beta_T: float = 0.0,
                extra_potential: float = 0.0) -> float:
    p = np.asarray(p, dtype=float)
    A = np.asarray(A, dtype=float)
    eff = p - charge * A
    kinetic = float(np.dot(eff, eff) / (2.0 * max(float(mass), EPS)))
    potential = float(charge * phi + alpha_D * D ** 2 + beta_T * T_reg ** 2 + extra_potential)
    return kinetic + potential


def effective_hamiltonian(p: Sequence[float], A_vec: Sequence[float], charge: float,
                          phi_struct: float, phi_aff: float, mass: float,
                          D: float = 0.0, T_reg: float = 0.0,
                          gamma_A: float = 1.0, beta_C: float = 0.0,
                          beta_L: float = 0.0, cognitive_cost: float = 0.0,
                          legitimacy: float = 1.0,
                          alpha_D: float = 0.0, beta_T: float = 0.0,
                          extra_potential: float = 0.0) -> float:
    phi_eff = phi_struct + gamma_A * phi_aff
    aff_penalty = beta_C * cognitive_cost ** 2 + beta_L * (1.0 - legitimacy) ** 2
    return hamiltonian(
        p=np.asarray(p, dtype=float),
        A=np.asarray(A_vec, dtype=float),
        charge=charge,
        phi=phi_eff,
        mass=mass,
        D=D,
        T_reg=T_reg,
        alpha_D=alpha_D,
        beta_T=beta_T,
        extra_potential=extra_potential + aff_penalty,
    )


def lorentz_like_force(charge: float, E: Sequence[float], velocity: Sequence[float], B: Sequence[float]) -> np.ndarray:
    E = np.asarray(E, dtype=float)
    v = np.asarray(velocity, dtype=float)
    B = np.asarray(B, dtype=float)
    return float(charge) * (E + np.cross(v, B))


# ==========================================================
# Node / edge state definitions
# ==========================================================
@dataclass
class NodeState:
    node_id: str
    regime: str
    q_state: np.ndarray
    S: float
    D: float
    H: float
    U: float
    Theta: float
    salience: float = 0.5
    feasibility: float = 0.5
    legitimacy: float = 0.5
    reward: float = 0.5
    mutuality: float = 0.5
    cognitive_cost: float = 0.5
    affordance_bias: float = 0.0
    charge: float = 0.0
    mass: float = 1.0

    @property
    def F(self) -> float:
        return free_energy(self.U, self.Theta, self.S)

    @property
    def phi(self) -> float:
        return scalar_potential(self.F)

    @property
    def A_logit(self) -> float:
        return affordance_logit(
            salience=self.salience,
            feasibility=self.feasibility,
            legitimacy=self.legitimacy,
            reward=self.reward,
            mutuality=self.mutuality,
            cognitive_cost=self.cognitive_cost,
            bias=self.affordance_bias,
        )

    @property
    def A(self) -> float:
        return sigmoid(self.A_logit)

    @property
    def phi_A(self) -> float:
        return affordance_potential(self.A)

    @property
    def phi_total(self) -> float:
        return total_potential(self.phi, self.phi_A)

    @property
    def cognitive_barrier(self) -> float:
        return cognitive_barrier_proxy(
            cognitive_cost=self.cognitive_cost,
            lack_of_legitimacy=max(0.0, 1.0 - self.legitimacy),
            uncertainty=max(0.0, 0.5 - self.feasibility),
        )

    @property
    def behavior_p(self) -> float:
        return behavior_probability(self.A, self.F, self.H, self.cognitive_barrier)


@dataclass
class EdgeState:
    src: str
    dst: str
    W: np.ndarray
    length: float = 1.0
    LA: float = 1.0
    LB: float = 1.0
    delta_topo: float = 0.0
    delta_op: float = 0.0
    gauge_dist: float = 0.0
    kappa: float = 1.0
    lam: float = 1.0
    eta: float = 1.0
    normative_barrier: float = 0.0
    cognitive_barrier: float = 0.0
    affordance_gain: float = 0.0
    affordance_alignment: float = 1.0

    @property
    def T_reg_base(self) -> float:
        return regime_tension(
            self.LA, self.LB, self.delta_topo, self.delta_op,
            kappa=self.kappa, lam=self.lam,
            gauge_dist=self.gauge_dist, eta=self.eta,
        )

    @property
    def T_reg(self) -> float:
        social_term = self.normative_barrier + self.cognitive_barrier + max(0.0, 1.0 - self.affordance_alignment)
        return float(self.T_reg_base + social_term)


# ==========================================================
# DataFrame / graph level helpers
# ==========================================================

def nodes_from_dataframe(df: pd.DataFrame) -> Dict[str, NodeState]:
    out: Dict[str, NodeState] = {}
    for _, r in df.iterrows():
        out[str(r['node_id'])] = NodeState(
            node_id=str(r['node_id']),
            regime=str(r.get('regime', 'default')),
            q_state=q(float(r['q_w']), float(r['q_x']), float(r['q_y']), float(r['q_z'])),
            S=float(r['S']), D=float(r['D']), H=float(r.get('H', 0.5)),
            U=float(r['U']), Theta=float(r['Theta']),
            salience=float(r.get('salience', 0.5)),
            feasibility=float(r.get('feasibility', 0.5)),
            legitimacy=float(r.get('legitimacy', 0.5)),
            reward=float(r.get('reward', 0.5)),
            mutuality=float(r.get('mutuality', 0.5)),
            cognitive_cost=float(r.get('cognitive_cost', 0.5)),
            affordance_bias=float(r.get('affordance_bias', 0.0)),
            charge=float(r.get('charge', 0.0)), mass=float(r.get('mass', 1.0)),
        )
    return out


def edges_from_dataframe(df: pd.DataFrame) -> List[EdgeState]:
    out: List[EdgeState] = []
    for _, r in df.iterrows():
        out.append(EdgeState(
            src=str(r['src']), dst=str(r['dst']),
            W=q(float(r['W_w']), float(r['W_x']), float(r['W_y']), float(r['W_z'])),
            length=float(r.get('length', 1.0)),
            LA=float(r.get('L_A', 1.0)), LB=float(r.get('L_B', 1.0)),
            delta_topo=float(r.get('delta_topo', 0.0)),
            delta_op=float(r.get('delta_op', 0.0)),
            gauge_dist=float(r.get('gauge_dist', 0.0)),
            kappa=float(r.get('kappa', 1.0)), lam=float(r.get('lambda', 1.0)),
            eta=float(r.get('eta', 1.0)),
            normative_barrier=float(r.get('normative_barrier', 0.0)),
            cognitive_barrier=float(r.get('cognitive_barrier', 0.0)),
            affordance_gain=float(r.get('affordance_gain', 0.0)),
            affordance_alignment=float(r.get('affordance_alignment', 1.0)),
        ))
    return out


def build_graph(node_df: pd.DataFrame, edge_df: pd.DataFrame, gamma_aff: float = 1.0):
    if nx is None:
        raise ImportError('networkx is required for graph-level analysis in this reference implementation.')
    nodes = nodes_from_dataframe(node_df)
    edges = edges_from_dataframe(edge_df)
    G = nx.DiGraph()
    for node_id, st in nodes.items():
        G.add_node(
            node_id,
            state=st,
            F=st.F,
            phi=st.phi,
            A=st.A,
            phi_A=st.phi_A,
            phi_total=st.phi_total,
            behavior_p=st.behavior_p,
        )
    for e in edges:
        src_state = nodes[e.src]
        dst_state = nodes[e.dst]
        sem_res = semantic_residual(dst_state.q_state, e.W, src_state.q_state)
        eps, mu, sigma = material_params((src_state.S + dst_state.S) / 2.0, e.T_reg)
        v = propagation_speed(eps, mu)
        Z = impedance(eps, mu)
        alpha = attenuation_alpha(eps, mu, sigma)
        E_struct = electric_field(src_state.phi, dst_state.phi, dA_dt=0.0, length=e.length)
        E_aff = electric_field(src_state.phi_A, dst_state.phi_A, dA_dt=0.0, length=e.length)
        E_total = E_struct + gamma_aff * E_aff
        aff_flow = edge_affordance_flow(
            A_src=src_state.A,
            affordance_gain=e.affordance_gain,
            normative_barrier=e.normative_barrier,
            cognitive_barrier=e.cognitive_barrier,
            semantic_penalty=sem_res * 0.1,
            length=e.length,
        )
        G.add_edge(
            e.src, e.dst,
            edge_state=e,
            T_reg=e.T_reg,
            T_reg_base=e.T_reg_base,
            A_vec=vector_potential(e.W),
            E_struct=E_struct,
            E_aff=E_aff,
            E_total=E_total,
            epsilon=eps,
            mu=mu,
            sigma=sigma,
            speed=v,
            impedance=Z,
            alpha=alpha,
            semantic_residual=sem_res,
            affordance_flow=aff_flow,
        )
    return G


def cycle_magnetic_fields(G) -> List[Dict[str, object]]:
    if nx is None:
        raise ImportError('networkx is required for cycle analysis.')
    out = []
    for cyc in nx.simple_cycles(G):
        if len(cyc) < 2:
            continue
        cyc_edges = []
        cyc_pairs = []
        for i in range(len(cyc)):
            src = cyc[i]
            dst = cyc[(i + 1) % len(cyc)]
            data = G.get_edge_data(src, dst)
            if data is None:
                cyc_edges = []
                break
            cyc_edges.append(data['edge_state'].W)
            cyc_pairs.append((src, dst))
        if cyc_edges:
            B = magnetic_field(cyc_edges, loop_scale=len(cyc_edges))
            out.append({
                'cycle': cyc_pairs,
                'B': B,
                'gauge_distortion': gauge_distortion(cyc_edges),
            })
    return out


def best_path_metrics(G, src: str, dst: str, amplitude0: float = 1.0) -> Dict[str, float]:
    if nx is None:
        raise ImportError('networkx is required for path analysis.')
    path = nx.shortest_path(G, src, dst, weight=lambda u, v, d: d['edge_state'].length / max(d['speed'], EPS))
    lengths, speeds, alphas = [], [], []
    edge_impedances = []
    aff_flows = []
    tensions = []
    for i in range(len(path) - 1):
        d = G[path[i]][path[i + 1]]
        lengths.append(d['edge_state'].length)
        speeds.append(d['speed'])
        alphas.append(d['alpha'])
        edge_impedances.append(d['impedance'])
        aff_flows.append(d['affordance_flow'])
        tensions.append(d['T_reg'])
    arrival = arrival_time(lengths, speeds)
    amp = propagated_amplitude(amplitude0, lengths, alphas)
    trans = 1.0
    for z1, z2 in zip(edge_impedances[:-1], edge_impedances[1:]):
        trans *= transmission_coefficient(z1, z2)
    return {
        'arrival_time': float(arrival),
        'amplitude': float(amp),
        'transmission': float(trans),
        'path_length': float(np.sum(lengths)),
        'mean_affordance_flow': float(np.mean(aff_flows)) if aff_flows else 0.0,
        'mean_tension': float(np.mean(tensions)) if tensions else 0.0,
    }


def node_intervention_diagnostics(G, node_id: str,
                                  velocity: Sequence[float] = (0.0, 0.0, 0.0),
                                  gamma_A: float = 1.0) -> Dict[str, float]:
    st: NodeState = G.nodes[node_id]['state']
    A_vec = np.zeros(3, dtype=float)
    E_sum = np.zeros(3, dtype=float)
    B_sum = np.zeros(3, dtype=float)
    T_local = []
    for _, _, data in G.out_edges(node_id, data=True):
        A_vec += np.asarray(data['A_vec'], dtype=float)
        E_sum += np.asarray([data['E_total'], 0.0, 0.0], dtype=float)
        T_local.append(float(data['T_reg']))
    for _, _, data in G.in_edges(node_id, data=True):
        A_vec += np.asarray(data['A_vec'], dtype=float)
        E_sum += np.asarray([data['E_total'], 0.0, 0.0], dtype=float)
        T_local.append(float(data['T_reg']))
    for cyc in cycle_magnetic_fields(G):
        if any(a == node_id or b == node_id for a, b in cyc['cycle']):
            B_sum += np.asarray(cyc['B'], dtype=float)
    force = lorentz_like_force(st.charge, E_sum, velocity, B_sum)
    Hm = effective_hamiltonian(
        p=np.asarray(velocity, dtype=float),
        A_vec=A_vec,
        charge=st.charge,
        phi_struct=st.phi,
        phi_aff=st.phi_A,
        mass=st.mass,
        D=st.D,
        T_reg=float(np.mean(T_local)) if T_local else 0.0,
        gamma_A=gamma_A,
        beta_C=0.2,
        beta_L=0.2,
        cognitive_cost=st.cognitive_cost,
        legitimacy=st.legitimacy,
        alpha_D=0.1,
        beta_T=0.1,
    )
    return {
        'Fx': float(force[0]),
        'Fy': float(force[1]),
        'Fz': float(force[2]),
        'Hamiltonian': float(Hm),
        'A': float(st.A),
        'phi_A': float(st.phi_A),
        'behavior_p': float(st.behavior_p),
    }


def rank_interventions(node_table: pd.DataFrame, edge_table: pd.DataFrame,
                       top_n: int = 10) -> pd.DataFrame:
    """
    Heuristic ranking of intervention points.
    Lower affordance, higher tension, higher cognitive burden => higher priority.
    """
    edge_burden = edge_table.groupby('dst', as_index=False)['T_reg'].mean().rename(columns={'T_reg': 'incoming_T_reg'})
    merged = node_table.merge(edge_burden, how='left', left_on='node_id', right_on='dst')
    merged['incoming_T_reg'] = merged['incoming_T_reg'].fillna(0.0)
    merged['priority_score'] = (
        (1.0 - merged['A']) * 2.0
        + merged['incoming_T_reg'] * 0.7
        + merged['cognitive_cost'] * 1.0
        + (1.0 - merged['legitimacy']) * 1.2
    )
    cols = [
        'node_id', 'A', 'behavior_p', 'cognitive_cost', 'legitimacy',
        'incoming_T_reg', 'priority_score'
    ]
    return merged.sort_values('priority_score', ascending=False)[cols].head(top_n).reset_index(drop=True)


def analyze_sdhat_series(df: pd.DataFrame, a_drop: float = -0.05, h_drop: float = -0.03,
                         d_var_q: float = 0.75, lookback: int = 5) -> pd.DataFrame:
    """
    Optional trajectory monitor that includes A.
    Input columns: t, S, D, H, A
    """
    out = df.copy()
    out['dA'] = out['A'].diff()
    out['dH'] = out['H'].diff()
    out['D_var'] = out['D'].rolling(lookback).var()
    thr = out['D_var'].quantile(d_var_q) if len(out) else np.nan
    out['alert_bad_transition'] = (out['dA'] < a_drop) & (out['dH'] < h_drop) & (out['D_var'] > thr)
    out['alert_good_transition'] = (out['dA'] > abs(a_drop)) & (out['H'].diff() >= 0.0)
    return out


def analyze_semantic_field(node_df: pd.DataFrame, edge_df: pd.DataFrame,
                           src: Optional[str] = None, dst: Optional[str] = None) -> Dict[str, object]:
    G = build_graph(node_df, edge_df)
    node_table = []
    for node_id, data in G.nodes(data=True):
        st: NodeState = data['state']
        node_table.append({
            'node_id': node_id,
            'regime': st.regime,
            'S': st.S,
            'D': st.D,
            'H': st.H,
            'U': st.U,
            'Theta': st.Theta,
            'F': st.F,
            'phi': st.phi,
            'A': st.A,
            'phi_A': st.phi_A,
            'phi_total': st.phi_total,
            'salience': st.salience,
            'feasibility': st.feasibility,
            'legitimacy': st.legitimacy,
            'reward': st.reward,
            'mutuality': st.mutuality,
            'cognitive_cost': st.cognitive_cost,
            'behavior_p': st.behavior_p,
            'charge': st.charge,
            'mass': st.mass,
        })
    edge_table = []
    for i, j, data in G.edges(data=True):
        edge_table.append({
            'src': i,
            'dst': j,
            'T_reg_base': data['T_reg_base'],
            'T_reg': data['T_reg'],
            'E_struct': data['E_struct'],
            'E_aff': data['E_aff'],
            'E_total': data['E_total'],
            'speed': data['speed'],
            'impedance': data['impedance'],
            'alpha': data['alpha'],
            'affordance_flow': data['affordance_flow'],
            'semantic_residual': data['semantic_residual'],
            'normative_barrier': data['edge_state'].normative_barrier,
            'cognitive_barrier': data['edge_state'].cognitive_barrier,
            'affordance_gain': data['edge_state'].affordance_gain,
        })
    nodes_df = pd.DataFrame(node_table)
    edges_df = pd.DataFrame(edge_table)
    result = {
        'nodes': nodes_df,
        'edges': edges_df,
        'cycles': cycle_magnetic_fields(G),
        'intervention_priority': rank_interventions(nodes_df, edges_df),
        'node_diagnostics': {
            node_id: node_intervention_diagnostics(G, node_id)
            for node_id in G.nodes
        },
    }
    if src is not None and dst is not None:
        result['best_path'] = best_path_metrics(G, src, dst)
    return result


# ==========================================================
# Example usage
# ==========================================================
if __name__ == '__main__':
    node_df = pd.DataFrame([
        {
            'node_id': 'HQ', 'regime': 'A', 'q_w': 1, 'q_x': 0.10, 'q_y': 0.00, 'q_z': 0.00,
            'S': 0.40, 'D': 1.25, 'H': 0.62, 'U': 1.20, 'Theta': 0.70,
            'salience': 0.62, 'feasibility': 0.68, 'legitimacy': 0.72,
            'reward': 0.74, 'mutuality': 0.61, 'cognitive_cost': 0.28,
            'charge': 0.8, 'mass': 1.0,
        },
        {
            'node_id': 'Store', 'regime': 'B', 'q_w': 1, 'q_x': 0.00, 'q_y': 0.20, 'q_z': 0.00,
            'S': 0.75, 'D': 1.55, 'H': 0.45, 'U': 0.95, 'Theta': 0.80,
            'salience': 0.58, 'feasibility': 0.44, 'legitimacy': 0.36,
            'reward': 0.79, 'mutuality': 0.42, 'cognitive_cost': 0.66,
            'charge': 0.2, 'mass': 1.2,
        },
        {
            'node_id': 'Ops', 'regime': 'B', 'q_w': 1, 'q_x': 0.05, 'q_y': 0.05, 'q_z': 0.10,
            'S': 0.58, 'D': 1.40, 'H': 0.50, 'U': 1.05, 'Theta': 0.72,
            'salience': 0.51, 'feasibility': 0.49, 'legitimacy': 0.40,
            'reward': 0.60, 'mutuality': 0.47, 'cognitive_cost': 0.52,
            'charge': -0.3, 'mass': 1.1,
        },
        {
            'node_id': 'FieldStaff', 'regime': 'C', 'q_w': 1, 'q_x': 0.10, 'q_y': 0.12, 'q_z': 0.18,
            'S': 0.62, 'D': 1.47, 'H': 0.55, 'U': 1.01, 'Theta': 0.73,
            'salience': 0.65, 'feasibility': 0.55, 'legitimacy': 0.48,
            'reward': 0.70, 'mutuality': 0.58, 'cognitive_cost': 0.46,
            'charge': 0.1, 'mass': 1.0,
        },
    ])

    edge_df = pd.DataFrame([
        {
            'src': 'HQ', 'dst': 'Store', 'W_w': 1, 'W_x': 0.10, 'W_y': 0.05, 'W_z': 0.00,
            'length': 1.2, 'L_A': 1.0, 'L_B': 1.2, 'delta_topo': 0.10, 'delta_op': 0.20,
            'gauge_dist': 0.04, 'normative_barrier': 0.15, 'cognitive_barrier': 0.10,
            'affordance_gain': 0.05, 'affordance_alignment': 0.80,
        },
        {
            'src': 'Store', 'dst': 'Ops', 'W_w': 1, 'W_x': 0.00, 'W_y': 0.10, 'W_z': 0.05,
            'length': 1.0, 'L_A': 1.2, 'L_B': 1.0, 'delta_topo': 0.15, 'delta_op': 0.25,
            'gauge_dist': 0.08, 'normative_barrier': 0.20, 'cognitive_barrier': 0.18,
            'affordance_gain': -0.04, 'affordance_alignment': 0.70,
        },
        {
            'src': 'Ops', 'dst': 'FieldStaff', 'W_w': 1, 'W_x': 0.03, 'W_y': 0.02, 'W_z': 0.08,
            'length': 1.1, 'L_A': 1.0, 'L_B': 0.9, 'delta_topo': 0.08, 'delta_op': 0.12,
            'gauge_dist': 0.05, 'normative_barrier': 0.12, 'cognitive_barrier': 0.14,
            'affordance_gain': 0.06, 'affordance_alignment': 0.84,
        },
        {
            'src': 'FieldStaff', 'dst': 'Store', 'W_w': 1, 'W_x': 0.02, 'W_y': 0.06, 'W_z': 0.10,
            'length': 1.3, 'L_A': 0.9, 'L_B': 1.2, 'delta_topo': 0.11, 'delta_op': 0.18,
            'gauge_dist': 0.06, 'normative_barrier': 0.10, 'cognitive_barrier': 0.16,
            'affordance_gain': 0.03, 'affordance_alignment': 0.76,
        },
    ])

    result = analyze_semantic_field(node_df, edge_df, src='HQ', dst='FieldStaff')

    print('\n=== Nodes ===')
    print(result['nodes'].round(4).to_string(index=False))

    print('\n=== Edges ===')
    print(result['edges'].round(4).to_string(index=False))

    print('\n=== Intervention Priority ===')
    print(result['intervention_priority'].round(4).to_string(index=False))

    print('\n=== Best Path ===')
    print({k: round(v, 4) for k, v in result['best_path'].items()})

    print('\n=== Node Diagnostics ===')
    for node_id, diag in result['node_diagnostics'].items():
        print(node_id, {k: round(v, 4) for k, v in diag.items()})
