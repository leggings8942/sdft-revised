"""
sdft_v02_geometry.py
SDFT revised v0.2 — 多様体幾何構造

k-NN グラフをベースとした離散リーマン幾何の計算を提供する。

【提供する幾何量】
  k-NN グラフ        軌道点群の近傍構造            Derived
  測地距離 d_g        Dijkstra 最短経路長           Derived
  Ollivier-Ricci κ   離散リッチ曲率の平均          Derived
  スペクトルギャップ λ₂  グラフ Laplacian Fiedler 値   Derived (= F_geom)

設計原則：
  1. 連続版のリーマン計量 g を陽に構成しない（ノイズに弱いため）
  2. グラフ上の組み合わせ的手法で曲率・距離を定義する
  3. Ollivier-Ricci 曲率は Wasserstein-1 距離で定義される
     κ(x,y) = 1 - W₁(m_x, m_y) / d(x,y)
  4. スペクトルギャップ λ₂ は正規化 Laplacian の第 2 最小固有値
"""

from __future__ import annotations

import math
from typing import List, Tuple, Optional, Dict

import numpy as np
from scipy.spatial import KDTree
from scipy.sparse import csr_matrix, diags
from scipy.sparse.linalg import eigsh
from scipy.sparse.csgraph import shortest_path
from scipy.optimize import linear_sum_assignment

EPS = 1e-12


# ============================================================
# k-NN グラフ構築
# ============================================================

def build_knn_graph(
    points: np.ndarray,
    k: int = 10,
) -> Tuple[np.ndarray, np.ndarray, KDTree]:
    """
    k-NN グラフを構築する。

    Returns
    -------
    distances : (N, k+1) ndarray  自分を含む距離
    indices   : (N, k+1) ndarray  自分を含むインデックス
    tree      : KDTree
    """
    N = len(points)
    k_actual = min(k + 1, N)
    tree = KDTree(points)
    distances, indices = tree.query(points, k=k_actual)
    return distances, indices, tree


def build_adjacency_matrix(
    points: np.ndarray,
    k: int = 10,
    sigma: Optional[float] = None,
) -> csr_matrix:
    """
    k-NN グラフの重み付き隣接行列を構築する。

    重み: w_ij = exp(-‖x_i - x_j‖² / (2σ²))
    σ: 指定なしの場合は中央値距離を使用。

    Returns: (N, N) sparse matrix
    """
    N = len(points)
    dists, indices, _ = build_knn_graph(points, k=k)

    if sigma is None:
        nonzero_dists = dists[:, 1:]
        sigma = float(np.median(nonzero_dists[nonzero_dists > 0])) + EPS

    rows, cols, vals = [], [], []
    for i in range(N):
        for j_idx in range(1, dists.shape[1]):  # skip self
            j = indices[i, j_idx]
            d = dists[i, j_idx]
            w = math.exp(-d ** 2 / (2 * sigma ** 2))
            rows.append(i)
            cols.append(j)
            vals.append(w)
            # 対称化
            rows.append(j)
            cols.append(i)
            vals.append(w)

    W = csr_matrix((vals, (rows, cols)), shape=(N, N))
    # 重複エントリの平均（対称化）
    W = (W + W.T) / 2.0
    return W


# ============================================================
# グラフ Laplacian
# ============================================================

def build_graph_laplacian(
    points: np.ndarray,
    k: int = 10,
    normalized: bool = True,
) -> Tuple[csr_matrix, csr_matrix]:
    """
    正規化グラフ Laplacian を構築する。

    L_norm = I - D^{-1/2} W D^{-1/2}

    Returns: (L, W)
    """
    W = build_adjacency_matrix(points, k=k)
    N = W.shape[0]
    degrees = np.array(W.sum(axis=1)).flatten()
    degrees = np.maximum(degrees, EPS)

    if normalized:
        D_inv_sqrt = diags(1.0 / np.sqrt(degrees))
        L = diags(np.ones(N)) - D_inv_sqrt @ W @ D_inv_sqrt
    else:
        D = diags(degrees)
        L = D - W

    return L, W


# ============================================================
# スペクトルギャップ λ₂（= F_geom）
# ============================================================

def calc_spectral_gap(
    points: np.ndarray,
    k: int = 10,
) -> float:
    """
    正規化グラフ Laplacian の第 2 最小固有値（Fiedler 値）を計算する。

    λ₂ = F_geom（構造エネルギー）

    λ₂ 大：系の構造が密で安定（結合が強い）
    λ₂ 小：構造が疎で崩壊寸前（クラスタが分裂しかけ）

    Cheeger 不等式：λ₂/2 ≤ h(G) ≤ √(2λ₂)
      h(G) はグラフの Cheeger 定数（多様体の等周定数の離散版）

    証拠分類: Derived
    """
    N = len(points)
    if N < 4:
        return 0.0

    L, _ = build_graph_laplacian(points, k=k)

    n_eigs = min(3, N - 1)
    try:
        eigenvalues, _ = eigsh(L, k=n_eigs, which='SM')
        eigenvalues = np.sort(np.real(eigenvalues))
        # λ₀ ≈ 0（連結成分）、λ₁ = Fiedler 値
        if len(eigenvalues) >= 2:
            return float(max(eigenvalues[1], 0.0))
    except Exception:
        pass

    return 0.0


# ============================================================
# 測地距離（Dijkstra on k-NN graph）
# ============================================================

def calc_geodesic_distance_matrix(
    points: np.ndarray,
    k: int = 10,
) -> np.ndarray:
    """
    k-NN グラフ上の最短経路距離行列を計算する。

    Isomap (Tenenbaum et al., 2000) と同じ手法。

    Returns: (N, N) ndarray
    """
    N = len(points)
    dists, indices, _ = build_knn_graph(points, k=k)

    # 疎距離行列を構築
    rows, cols, vals = [], [], []
    for i in range(N):
        for j_idx in range(1, dists.shape[1]):
            j = indices[i, j_idx]
            d = dists[i, j_idx]
            rows.append(i)
            cols.append(j)
            vals.append(d)

    graph = csr_matrix((vals, (rows, cols)), shape=(N, N))
    # 対称化
    graph = graph.minimum(graph.T)
    mask = graph > 0
    graph[mask] = graph[mask]
    graph = (graph + graph.T) / 2

    dist_matrix = shortest_path(graph, directed=False)
    return dist_matrix


def calc_geodesic_distance(
    points_A: np.ndarray,
    points_B: np.ndarray,
    k: int = 10,
) -> float:
    """
    2つの点群のセントロイド間の測地距離を計算する。

    手順：
      1. 2点群を結合して k-NN グラフを作る
      2. 各セントロイドの最近傍を探す
      3. Dijkstra で最短経路距離を求める

    証拠分類: Derived
    """
    combined = np.vstack([points_A, points_B])
    N_A = len(points_A)

    centroid_A = points_A.mean(axis=0)
    centroid_B = points_B.mean(axis=0)

    tree = KDTree(combined)
    _, idx_A = tree.query(centroid_A, k=1)
    _, idx_B = tree.query(centroid_B, k=1)

    idx_A = int(np.atleast_1d(idx_A)[0])
    idx_B = int(np.atleast_1d(idx_B)[0])

    dists, indices, _ = build_knn_graph(combined, k=k)
    N = len(combined)
    rows, cols, vals = [], [], []
    for i in range(N):
        for j_idx in range(1, dists.shape[1]):
            j = indices[i, j_idx]
            d = dists[i, j_idx]
            rows.append(i)
            cols.append(j)
            vals.append(d)

    graph = csr_matrix((vals, (rows, cols)), shape=(N, N))
    graph = (graph + graph.T) / 2

    dist_mat = shortest_path(graph, directed=False,
                             indices=[idx_A])
    geodesic = float(dist_mat[0, idx_B])
    if np.isinf(geodesic):
        geodesic = float(np.linalg.norm(centroid_A - centroid_B))

    return geodesic


# ============================================================
# Ollivier-Ricci 曲率
# ============================================================

def _wasserstein1_uniform(
    neighbors_x: np.ndarray,
    neighbors_y: np.ndarray,
    all_points: np.ndarray,
) -> float:
    """
    2つの k-近傍集合の一様分布間の Wasserstein-1 距離を
    最適割当問題として計算する。

    W₁(m_x, m_y) = (1/k) min_{σ} Σ_i d(n_x(i), n_y(σ(i)))
    """
    k_x = len(neighbors_x)
    k_y = len(neighbors_y)
    k = min(k_x, k_y)
    if k == 0:
        return 0.0

    # コスト行列
    pts_x = all_points[neighbors_x[:k]]
    pts_y = all_points[neighbors_y[:k]]
    cost = np.linalg.norm(pts_x[:, None, :] - pts_y[None, :, :], axis=2)

    row_ind, col_ind = linear_sum_assignment(cost)
    return float(cost[row_ind, col_ind].sum()) / k


def calc_ollivier_ricci(
    points: np.ndarray,
    k: int = 10,
    n_edges: int = 500,
) -> float:
    """
    Ollivier-Ricci 曲率の平均値を計算する。

    κ(x, y) = 1 - W₁(m_x, m_y) / d(x, y)

    κ > 0：正曲率（近傍が収束・クラスタ的）
    κ = 0：平坦
    κ < 0：負曲率（近傍が発散・ハイパーボリック）

    References: Ollivier (2009), Lin-Lu-Yau (2011)
    証拠分類: Derived
    """
    N = len(points)
    if N < k + 2:
        return 0.0

    dists, indices, _ = build_knn_graph(points, k=k)
    curvatures = []

    # ランダムにエッジを選んで曲率を計算
    edge_count = min(n_edges, N * k)
    sampled = set()
    attempts = 0

    while len(sampled) < edge_count and attempts < edge_count * 3:
        i = np.random.randint(0, N)
        j_idx = np.random.randint(1, dists.shape[1])
        j = indices[i, j_idx]
        edge = (min(i, j), max(i, j))
        if edge not in sampled:
            sampled.add(edge)

            d_ij = float(np.linalg.norm(points[i] - points[j]))
            if d_ij < EPS:
                continue

            # x の k-近傍（自分を除く）
            nbrs_x = indices[i, 1:]
            nbrs_y = indices[j, 1:]

            w1 = _wasserstein1_uniform(nbrs_x, nbrs_y, points)
            kappa = 1.0 - w1 / d_ij
            curvatures.append(kappa)
        attempts += 1

    if not curvatures:
        return 0.0

    return float(np.mean(curvatures))


# ============================================================
# 統合解析関数
# ============================================================

def analyze_geometry(
    embedded: np.ndarray,
    k: int = 10,
) -> dict:
    """
    埋め込み点群の幾何構造を計算する。

    Returns
    -------
    dict: spectral_gap (F_geom), curvature, geodesic_dist_matrix
    """
    F_geom   = calc_spectral_gap(embedded, k=k)
    kappa    = calc_ollivier_ricci(embedded, k=k)

    return {
        'F_geom':           F_geom,
        'curvature':        kappa,
        'F_geom_evidence':  'Derived',
        'curvature_evidence': 'Derived',
    }


# ============================================================
# 動作確認
# ============================================================

if __name__ == "__main__":
    np.random.seed(42)

    print("=" * 60)
    print("SDFT v0.2 Geometry — 動作確認")
    print("=" * 60)

    # ① 2D ガウス雲（正曲率的）
    pts_gauss = np.random.randn(300, 3) * 0.5
    print("\n--- 3D ガウス雲 ---")
    sg = calc_spectral_gap(pts_gauss, k=10)
    kr = calc_ollivier_ricci(pts_gauss, k=10, n_edges=200)
    print(f"  λ₂ (spectral gap) = {sg:.4f}")
    print(f"  κ  (Ollivier-Ricci) = {kr:.4f}  (期待: > 0, 球的)")

    # ② 2クラスタ（低 λ₂）
    cluster_1 = np.random.randn(150, 3) * 0.3 + [3, 0, 0]
    cluster_2 = np.random.randn(150, 3) * 0.3 + [-3, 0, 0]
    pts_2clust = np.vstack([cluster_1, cluster_2])
    print("\n--- 2クラスタ ---")
    sg2 = calc_spectral_gap(pts_2clust, k=10)
    kr2 = calc_ollivier_ricci(pts_2clust, k=10, n_edges=200)
    print(f"  λ₂ = {sg2:.4f}  (期待: ガウス雲より小さい)")
    print(f"  κ  = {kr2:.4f}")

    # ③ 測地距離
    d_geo = calc_geodesic_distance(cluster_1, cluster_2, k=10)
    d_euc = float(np.linalg.norm(cluster_1.mean(0) - cluster_2.mean(0)))
    print(f"\n--- 測地距離 ---")
    print(f"  d_geo = {d_geo:.4f}")
    print(f"  d_euc = {d_euc:.4f}  (d_geo ≥ d_euc)")

    print("\n" + "=" * 60)
    print("動作確認完了")
