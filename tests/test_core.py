"""
SDFT v0.2 テストスイート

数学的不変量と回帰テストを実行する。
"""

import sys
import os
import numpy as np

# スクリプトディレクトリをパスに追加
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'scripts'))


def test_takens_embedding():
    """Takens 埋め込みの基本テスト"""
    from sdft_v02_embedding import takens_embed, estimate_delay, estimate_embedding_dim

    # 正弦波の埋め込み
    t = np.linspace(0, 10 * np.pi, 500)
    x = np.sin(t)

    tau = estimate_delay(x.tolist())
    assert tau >= 1, f"τ should be ≥ 1, got {tau}"

    d = estimate_embedding_dim(x.tolist(), tau)
    assert 1 <= d <= 10, f"d should be in [1, 10], got {d}"

    emb = takens_embed(x.tolist(), d=3, tau=5)
    expected_N = len(x) - (3 - 1) * 5
    assert emb.shape == (expected_N, 3), f"Shape mismatch: {emb.shape}"

    print("  ✅ test_takens_embedding passed")


def test_volume_entropy():
    """S_g がスケールに応答するか"""
    from sdft_v02_embedding import calc_volume_entropy

    np.random.seed(42)
    compact = np.random.randn(200, 3) * 0.1
    spread  = np.random.randn(200, 3) * 5.0

    S_compact = calc_volume_entropy(compact)
    S_spread  = calc_volume_entropy(spread)

    assert S_spread > S_compact, (
        f"S_g(spread)={S_spread:.4f} should > S_g(compact)={S_compact:.4f}"
    )
    print(f"  ✅ test_volume_entropy passed (compact={S_compact:.2f}, spread={S_spread:.2f})")


def test_correlation_dim():
    """相関次元が理論値に近いか"""
    from sdft_v02_embedding import calc_correlation_dim

    np.random.seed(42)
    # 1 次元の曲線（期待 D ≈ 1）
    t = np.linspace(0, 10, 500)
    curve = np.column_stack([np.sin(t), np.cos(t), np.sin(2 * t)])
    D_curve = calc_correlation_dim(curve)
    assert 0.5 < D_curve < 2.0, f"Curve D should ≈ 1, got {D_curve:.4f}"

    # 3D ガウス雲（期待 D → 3 に近い）
    cloud = np.random.randn(500, 3)
    D_cloud = calc_correlation_dim(cloud)
    assert D_cloud > D_curve, (
        f"Cloud D={D_cloud:.4f} should > Curve D={D_curve:.4f}"
    )

    print(f"  ✅ test_correlation_dim passed (curve={D_curve:.2f}, cloud={D_cloud:.2f})")


def test_lyapunov_sign():
    """リアプノフ指数の符号テスト"""
    from sdft_v02_embedding import calc_lyapunov

    np.random.seed(42)

    # 安定な軌道（λ₁ < 0 を期待）
    t = np.linspace(0, 20 * np.pi, 800)
    stable = np.column_stack([
        np.sin(t) + 0.01 * np.random.randn(len(t)),
        np.cos(t) + 0.01 * np.random.randn(len(t)),
        np.sin(2 * t) + 0.01 * np.random.randn(len(t)),
    ])
    lam_stable = calc_lyapunov(stable)

    # ランダムウォーク（λ₁ > 0 を期待）
    rw = np.cumsum(np.random.randn(800, 3), axis=0)
    lam_rw = calc_lyapunov(rw)

    # ランダムウォークの方がリアプノフ指数が高い
    assert lam_rw > lam_stable, (
        f"Random walk λ₁={lam_rw:.4f} should > stable λ₁={lam_stable:.4f}"
    )

    print(f"  ✅ test_lyapunov_sign passed (stable={lam_stable:.4f}, rw={lam_rw:.4f})")


def test_spectral_gap():
    """スペクトルギャップの比較テスト"""
    from sdft_v02_geometry import calc_spectral_gap

    np.random.seed(42)

    # 単一クラスタ（高 λ₂）
    single = np.random.randn(200, 3) * 0.5
    sg_single = calc_spectral_gap(single, k=10)

    # 2 クラスタ（低 λ₂）
    c1 = np.random.randn(100, 3) * 0.3 + [5, 0, 0]
    c2 = np.random.randn(100, 3) * 0.3 + [-5, 0, 0]
    two_clusters = np.vstack([c1, c2])
    sg_two = calc_spectral_gap(two_clusters, k=10)

    assert sg_single > sg_two, (
        f"Single cluster λ₂={sg_single:.4f} should > "
        f"two clusters λ₂={sg_two:.4f}"
    )

    print(f"  ✅ test_spectral_gap passed (single={sg_single:.4f}, two={sg_two:.4f})")


def test_geodesic_distance():
    """測地距離が Euclidean 距離以上か"""
    from sdft_v02_geometry import calc_geodesic_distance

    np.random.seed(42)
    A = np.random.randn(100, 3) * 0.3 + [3, 0, 0]
    B = np.random.randn(100, 3) * 0.3 + [-3, 0, 0]

    d_geo = calc_geodesic_distance(A, B, k=10)
    d_euc = float(np.linalg.norm(A.mean(0) - B.mean(0)))

    assert d_geo >= d_euc * 0.9, (
        f"d_geo={d_geo:.4f} should ≥ d_euc={d_euc:.4f} (approx)"
    )

    print(f"  ✅ test_geodesic_distance passed (geo={d_geo:.2f}, euc={d_euc:.2f})")


def test_tension_same_distribution():
    """同一分布で 𝓣 ≈ 0"""
    from sdft_v02_tension import calc_regime_tension_vector

    np.random.seed(42)
    A = np.random.randn(200, 3) * 0.5
    B = np.random.randn(200, 3) * 0.5

    t = calc_regime_tension_vector(A, B, k=10)
    assert t.T1 < 0.5, f"T₁={t.T1:.4f} should be small for same dist"
    assert t.T4 < 0.8, f"T₄={t.T4:.4f} should be moderate for same dist"

    print(f"  ✅ test_tension_same_distribution passed (T₁={t.T1:.3f}, T₂={t.T2:.3f})")


def test_tension_different_distribution():
    """異なる分布で 𝓣 が大きい"""
    from sdft_v02_tension import calc_regime_tension_vector

    np.random.seed(42)
    A = np.random.randn(200, 3) * 0.3 + [5, 0, 0]
    B = np.random.randn(200, 3) * 2.0 + [-5, 0, 0]

    t_same = calc_regime_tension_vector(
        np.random.randn(200, 3), np.random.randn(200, 3), k=10
    )
    t_diff = calc_regime_tension_vector(A, B, k=10)

    assert t_diff.T1 > t_same.T1, "T₁ should be larger for different dists"
    assert t_diff.T2 > t_same.T2, "T₂ should be larger for different dists"

    print(f"  ✅ test_tension_different_distribution passed")


def test_grand_potential():
    """Φ = F_geom - μN の基本テスト"""
    from sdft_v02_potential import calc_grand_potential, calc_mu

    F = 0.15
    T2 = 0.8
    N = 100.0

    mu = calc_mu(T2)
    Phi = calc_grand_potential(F, mu, N)

    expected = F - T2 * N  # μ = T₂
    assert abs(Phi - expected) < 1e-10, f"Φ={Phi}, expected={expected}"

    # Φ は F_geom 増加で増加、μN 増加で減少
    Phi2 = calc_grand_potential(F + 0.1, mu, N)
    assert Phi2 > Phi, "Increasing F_geom should increase Φ"

    Phi3 = calc_grand_potential(F, mu, N + 10)
    assert Phi3 < Phi, "Increasing N should decrease Φ"

    print(f"  ✅ test_grand_potential passed (Φ={Phi:.4f})")


def test_auto_phi_threshold_sign():
    """⑤ バグ修正確認：負の Φ で threshold が Φ より低いか"""
    from sdft_v02_intervention import auto_phi_threshold, StateSnapshot

    # 正の Φ
    snap_pos = StateSnapshot(t=0.0, embedded_A=np.zeros((10, 3)),
                              embedded_B=np.zeros((10, 3)))
    snap_pos.Phi = 10.0
    th_pos, _ = auto_phi_threshold(snap_pos)
    assert th_pos < snap_pos.Phi, f"Positive Φ: threshold={th_pos} should < Φ={snap_pos.Phi}"

    # 負の Φ
    snap_neg = StateSnapshot(t=0.0, embedded_A=np.zeros((10, 3)),
                              embedded_B=np.zeros((10, 3)))
    snap_neg.Phi = -50.0
    th_neg, _ = auto_phi_threshold(snap_neg)
    assert th_neg < snap_neg.Phi, (
        f"Negative Φ: threshold={th_neg} should < Φ={snap_neg.Phi}"
    )

    # ゼロの Φ
    snap_zero = StateSnapshot(t=0.0, embedded_A=np.zeros((10, 3)),
                               embedded_B=np.zeros((10, 3)))
    snap_zero.Phi = 0.0
    th_zero, _ = auto_phi_threshold(snap_zero)
    assert th_zero < snap_zero.Phi, f"Zero Φ: threshold={th_zero} should < 0"

    print(f"  ✅ test_auto_phi_threshold_sign passed (pos_th={th_pos:.2f}, neg_th={th_neg:.2f})")


def test_phase_classification_order():
    """Frozen Order が Living System に shadow されないか"""
    from sdft_v02_potential import classify_phase

    # Frozen Order に分類されるべきケース
    # S_norm ≤ 0.3, D_norm ≤ 0.3, λ₁ < -0.1
    label, _, _ = classify_phase(S_g=-4.0, D=0.5, lambda_1=-0.2, d_embed=3)
    assert label == "Frozen Order", f"Expected Frozen Order, got {label}"

    # Living System に分類されるべきケース
    label2, _, _ = classify_phase(S_g=-1.0, D=1.5, lambda_1=-0.05, d_embed=3)
    assert label2 == "Living System", f"Expected Living System, got {label2}"

    # Collapse
    label3, _, _ = classify_phase(S_g=4.0, D=2.5, lambda_1=0.2, d_embed=3)
    assert label3 == "Collapse", f"Expected Collapse, got {label3}"

    print(f"  ✅ test_phase_classification_order passed")


def test_auto_tau():
    """τ がリアプノフ指数に応じて変わるか"""
    from sdft_v02_intervention import auto_tau

    tau_stable, _ = auto_tau(-0.2)
    tau_chaotic, _ = auto_tau(0.1)

    assert tau_stable > tau_chaotic, (
        f"Stable τ={tau_stable} should > chaotic τ={tau_chaotic}"
    )
    print(f"  ✅ test_auto_tau passed (stable={tau_stable}, chaotic={tau_chaotic})")


# ============================================================
# メイン
# ============================================================

if __name__ == "__main__":
    print("=" * 60)
    print("SDFT v0.2 Test Suite")
    print("=" * 60)

    tests = [
        test_takens_embedding,
        test_volume_entropy,
        test_correlation_dim,
        test_lyapunov_sign,
        test_spectral_gap,
        test_geodesic_distance,
        test_tension_same_distribution,
        test_tension_different_distribution,
        test_grand_potential,
        test_auto_phi_threshold_sign,
        test_phase_classification_order,
        test_auto_tau,
    ]

    passed = 0
    failed = 0
    for test in tests:
        try:
            test()
            passed += 1
        except Exception as e:
            print(f"  ❌ {test.__name__} FAILED: {e}")
            failed += 1

    print(f"\n{'=' * 60}")
    print(f"Results: {passed} passed, {failed} failed out of {len(tests)}")
    print("=" * 60)

    if failed > 0:
        sys.exit(1)
