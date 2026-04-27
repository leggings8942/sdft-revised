"""
SDFT Unified Complete v1.4 Reference Implementation
Gauge-Invariant Semantic Field / Maxwellian Propagation / Hamiltonian Intervention

This module is intentionally self-contained so that it can be embedded in an HTML
knowledge document and directly reused by an AI agent or an analyst.
It is a reference implementation, not a proof of mathematical completeness.
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
# Phase 1–4: Base utilities already used in Unified v1.3.1
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
            norm = (n - 1) / (((n - m - 1) // k) * k)
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
    sigma: dissipation (entropy + silo wall)
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


def lorentz_like_force(charge: float, E: Sequence[float], velocity: Sequence[float], B: Sequence[float]) -> np.ndarray:
    E = np.asarray(E, dtype=float)
    v = np.asarray(velocity, dtype=float)
    B = np.asarray(B, dtype=float)
    return float(charge) * (E + np.cross(v, B))


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
    charge: float = 0.0
    mass: float = 1.0

    @property
    def F(self) -> float:
        return free_energy(self.U, self.Theta, self.S)

    @property
    def phi(self) -> float:
        return scalar_potential(self.F)


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

    @property
    def T_reg(self) -> float:
        return regime_tension(self.LA, self.LB, self.delta_topo, self.delta_op,
                              kappa=self.kappa, lam=self.lam,
                              gauge_dist=self.gauge_dist, eta=self.eta)


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
        ))
    return out


def build_graph(node_df: pd.DataFrame, edge_df: pd.DataFrame):
    if nx is None:
        raise ImportError('networkx is required for graph-level analysis in this reference implementation.')
    nodes = nodes_from_dataframe(node_df)
    edges = edges_from_dataframe(edge_df)
    G = nx.DiGraph()
    for node_id, st in nodes.items():
        G.add_node(node_id, state=st, F=st.F, phi=st.phi)
    for e in edges:
        src_state = nodes[e.src]
        dst_state = nodes[e.dst]
        eps, mu, sigma = material_params((src_state.S + dst_state.S) / 2.0, e.T_reg)
        v = propagation_speed(eps, mu)
        Z = impedance(eps, mu)
        alpha = attenuation_alpha(eps, mu, sigma)
        E = electric_field(src_state.phi, dst_state.phi, dA_dt=0.0, length=e.length)
        G.add_edge(e.src, e.dst,
                   edge_state=e,
                   T_reg=e.T_reg,
                   A_vec=vector_potential(e.W),
                   E=E,
                   epsilon=eps,
                   mu=mu,
                   sigma=sigma,
                   speed=v,
                   impedance=Z,
                   alpha=alpha,
                   semantic_residual=semantic_residual(dst_state.q_state, e.W, src_state.q_state))
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
    for i in range(len(path) - 1):
        d = G[path[i]][path[i + 1]]
        lengths.append(d['edge_state'].length)
        speeds.append(d['speed'])
        alphas.append(d['alpha'])
        edge_impedances.append(d['impedance'])
    arrival = arrival_time(lengths, speeds)
    amp = propagated_amplitude(amplitude0, lengths, alphas)
    trans = 1.0
    for z1, z2 in zip(edge_impedances[:-1], edge_impedances[1:]):
        trans *= transmission_coefficient(z1, z2)
    return {
        'arrival_time': float(arrival),
        'amplitude': float(amp),
        'transmission': float(trans),
        'path_length': float(np.sum(lengths)) if lengths else 0.0,
    }


def charged_node_effect(G, node_id: str, velocity: Sequence[float] = (0.0, 0.0, 0.0)) -> Dict[str, float]:
    st: NodeState = G.nodes[node_id]['state']
    A = np.zeros(3, dtype=float)
    E_sum = np.zeros(3, dtype=float)
    B_sum = np.zeros(3, dtype=float)
    for _, dst, data in G.out_edges(node_id, data=True):
        A += np.asarray(data['A_vec'], dtype=float)
        E_sum += np.asarray([data['E'], 0.0, 0.0], dtype=float)
    for cyc in cycle_magnetic_fields(G):
        if any(a == node_id for a, _ in cyc['cycle']):
            B_sum += np.asarray(cyc['B'], dtype=float)
    force = lorentz_like_force(st.charge, E_sum, velocity, B_sum)
    Hm = hamiltonian(np.asarray(velocity, dtype=float), A, st.charge, st.phi, st.mass,
                     D=st.D, T_reg=0.0, alpha_D=0.1, beta_T=0.1)
    return {
        'Fx': float(force[0]), 'Fy': float(force[1]), 'Fz': float(force[2]),
        'Hamiltonian': float(Hm),
    }


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
            'charge': st.charge,
            'mass': st.mass,
        })
    edge_table = []
    for i, j, data in G.edges(data=True):
        edge_table.append({
            'src': i,
            'dst': j,
            'T_reg': data['T_reg'],
            'E': data['E'],
            'speed': data['speed'],
            'impedance': data['impedance'],
            'alpha': data['alpha'],
            'semantic_residual': data['semantic_residual'],
        })
    result = {
        'nodes': pd.DataFrame(node_table),
        'edges': pd.DataFrame(edge_table),
        'cycles': cycle_magnetic_fields(G),
    }
    if src is not None and dst is not None:
        result['best_path'] = best_path_metrics(G, src, dst)
    return result


# ==========================================================
# Example usage
# ==========================================================
if __name__ == '__main__':
    node_df = pd.DataFrame([
        {'node_id': 'HQ', 'regime': 'A', 'q_w': 1, 'q_x': 0.1, 'q_y': 0.0, 'q_z': 0.0,
         'S': 0.40, 'D': 1.25, 'H': 0.62, 'U': 1.20, 'Theta': 0.70, 'charge': 0.8, 'mass': 1.0},
        {'node_id': 'Store', 'regime': 'B', 'q_w': 1, 'q_x': 0.0, 'q_y': 0.2, 'q_z': 0.0,
         'S': 0.75, 'D': 1.55, 'H': 0.45, 'U': 0.95, 'Theta': 0.80, 'charge': 0.2, 'mass': 1.2},
        {'node_id': 'Ops', 'regime': 'B', 'q_w': 1, 'q_x': 0.05, 'q_y': 0.05, 'q_z': 0.1,
         'S': 0.58, 'D': 1.40, 'H': 0.50, 'U': 1.05, 'Theta': 0.72, 'charge': -0.3, 'mass': 1.1},
    ])

    edge_df = pd.DataFrame([
        {'src': 'HQ', 'dst': 'Store', 'W_w': 0.98, 'W_x': 0.05, 'W_y': 0.02, 'W_z': 0.00,
         'length': 2.0, 'L_A': 1.4, 'L_B': 1.0, 'delta_topo': 0.25, 'delta_op': 0.35},
        {'src': 'Store', 'dst': 'Ops', 'W_w': 0.97, 'W_x': 0.02, 'W_y': 0.06, 'W_z': 0.01,
         'length': 1.2, 'L_A': 1.0, 'L_B': 1.0, 'delta_topo': 0.10, 'delta_op': 0.18},
        {'src': 'Ops', 'dst': 'HQ', 'W_w': 0.96, 'W_x': 0.01, 'W_y': 0.03, 'W_z': 0.06,
         'length': 1.5, 'L_A': 1.0, 'L_B': 1.4, 'delta_topo': 0.22, 'delta_op': 0.28},
    ])

    result = analyze_semantic_field(node_df, edge_df, src='HQ', dst='Store')
    print('=== Nodes ===')
    print(result['nodes'])
    print('\n=== Edges ===')
    print(result['edges'])
    print('\n=== Cycles ===')
    for cyc in result['cycles']:
        print(cyc)
    print('\n=== Best Path ===')
    print(result['best_path'])
