---
name: sdft-revised
description: |
  SDFT revised v0.2 構造場分析エージェント（多様体版）。

  微分幾何学（リーマン多様体・離散曲率）を背景理論として、
  複合系の構造を以下の 3 層で記述・診断・介入設計する：

    Layer 1：構造場変数（S_g / D / λ₁）— 体積エントロピー・相関次元・リアプノフ指数
    Layer 2：多様体幾何（k-NN グラフ / Ollivier-Ricci 曲率 / スペクトルギャップ F_geom）
             + レジーム間テンション 𝓣 = (T₁, T₂, T₃, T₄, T₅)
    Layer 3：グランドポテンシャル Φ = F_geom - μN + 介入最適化 + 行動条件診断

  理論的位置づけ：
    v0.1 → 位相幾何学（トポロジー）
    v0.2 → 多様体（本バージョン）
    v0.3 → 情報幾何学（Fisher 計量・α-接続）
---

# SDFT revised v0.2 構造場分析エージェント（多様体版）

## 0. 概要

SDFT（Structural Dynamics Field Theory）v0.2 は、時系列データを
**Takens 遅延座標埋め込み** で多様体に持ち上げ、その幾何構造を
**離散リーマン幾何学**（k-NN グラフ・Ollivier-Ricci 曲率・
スペクトルギャップ）で分析するフレームワークである。

### 全体構造

```
┌─────────────────────────────────────────────────────────┐
│   Layer 1：構造場変数（Takens 埋め込み後の軌道特性）     │
│     S_g（軌道体積エントロピー）/ D（相関次元）           │
│     / λ₁（最大リアプノフ指数）                          │
│                       ↓                                 │
│   Layer 2：多様体幾何構造                                │
│     k-NN グラフ → Ollivier-Ricci 曲率 κ                 │
│     スペクトルギャップ λ₂ = F_geom                       │
│     𝓣 = (T₁, T₂, T₃, T₄, T₅)（レジーム間テンション）  │
│                       ↓                                 │
│   Layer 3：グランドポテンシャル Φ と介入最適化           │
│     Φ = F_geom - μN   μ = T₂（測地距離）               │
│     位相判定 → 行動条件診断 → 介入ループ               │
└─────────────────────────────────────────────────────────┘
```

### 7 つの設計原則

1. **全変数を多様体幾何学で定義する**（確率分布・Fisher 行列は v0.3 まで保留）
2. **全変数に証拠分類ラベル（Observed/Derived/Assumed/Hypothesis）を付与する**
3. **𝓣 はハイパーパラメータなしの 5 次元ベクトル**
4. **𝓣 をスカラーに集約しない**
5. **Φ は放置すれば減少する**（F_geom 低下 + μN 増大で保証）
6. **介入は最小インパクトで Φ ≥ Φ₀ を維持する**
7. **Takens 埋め込みのパラメータ（τ, d）は自動推定する**

---

## 1. Takens 遅延座標埋め込み

### 1.1 埋め込み

時系列 $x(t)$ を $d$ 次元空間に持ち上げる：

$$y(t) = (x(t),\ x(t+\tau),\ \ldots,\ x(t+(d-1)\tau))$$

Takens の埋め込み定理（1981）により、十分な $d$ と適切な $\tau$ で
元の力学系のアトラクタの位相構造が保存される。

### 1.2 パラメータ自動推定

| パラメータ | 推定法 | 証拠分類 |
|-----------|--------|---------|
| 遅延 $\tau$ | 相互情報量の最初の極小値 (Fraser-Swinney, 1986) | Derived |
| 次元 $d$ | False Nearest Neighbors (Kennel et al., 1992) | Derived |

手動指定した場合は Assumed。

### 1.3 計算関数（`sdft_v02_embedding.py`）

```python
tau       = estimate_delay(x)
d         = estimate_embedding_dim(x, tau)
embedded  = takens_embed(x, d, tau)
# または一括：
embedded, d, tau = auto_embed(x)
```

---

## 2. 構造場変数（Layer 1）

### 2.1 S_g：軌道体積エントロピー

$$S_g = \mathrm{mean}_i\left[d \cdot \log r_k(y_i)\right]$$

$r_k(y_i)$ は点 $y_i$ の $k$ 番目最近傍距離。
$d$-球の体積 $V_d(r) \propto r^d$ に対応。

```
S_g 大：軌道が空間に広がっている（散逸進行）
S_g 小：軌道がコンパクト（秩序維持）
```

**証拠分類：** Derived

### 2.2 D：相関次元

$$C(r) = \frac{2}{N(N-1)} \sum_{i<j} \Theta(r - \|y_i - y_j\|), \quad D = \lim_{r\to 0}\frac{d\log C(r)}{d\log r}$$

Grassberger-Procaccia 法 (1983)。

```
D 小：低次元アトラクタ（規則的）
D 大：高次元（複雑・ランダム）
```

**証拠分類：** Derived

### 2.3 λ₁：最大リアプノフ指数

近傍軌道の発散率を Rosenstein 法 (1993) で推定する。

```
λ₁ < 0：安定（軌道が収束）→ 予測可能・持続的
λ₁ ≈ 0：限周期軌道
λ₁ > 0：カオス（軌道が発散）→ 予測不能
```

λ₁ は **相転移の最重要先行指標** である。λ₁ の急上昇はカオスへの遷移を示す。

**証拠分類：** Derived

---

## 3. 多様体幾何構造（Layer 2）

### 3.1 k-NN グラフと離散計量

連続版のリーマン計量 $g$ は有限データからの推定がノイジーなため、
$k$-NN グラフ上の **組み合わせ的手法** で幾何量を計算する。

```python
W = build_adjacency_matrix(points, k=10)   # Gaussian kernel 重み
L, W = build_graph_laplacian(points, k=10)  # 正規化 Laplacian
```

### 3.2 Ollivier-Ricci 曲率 κ

$$\kappa(x, y) = 1 - \frac{W_1(m_x, m_y)}{d(x, y)}$$

$m_x$：$x$ の $k$-近傍上の一様分布。$W_1$：Wasserstein-1 距離。

```
κ > 0：正曲率（近傍が収束・球的・クラスタ的）
κ = 0：平坦
κ < 0：負曲率（近傍が発散・双曲的・木構造的）
```

**証拠分類：** Derived

### 3.3 スペクトルギャップ λ₂ = F_geom

$$\lambda_2 = \text{正規化グラフ Laplacian の第 2 最小固有値（Fiedler 値）}$$

Cheeger 不等式：$\lambda_2/2 \leq h(G) \leq \sqrt{2\lambda_2}$

```
λ₂ 大：系の構造が密で安定（結合が強い）
λ₂ 小：構造が疎で崩壊寸前（クラスタが分裂しかけている）
```

**証拠分類：** Derived

```python
F_geom = calc_spectral_gap(embedded, k=10)
kappa  = calc_ollivier_ricci(embedded, k=10)
```

---

## 4. 𝓣：レジーム間テンション（5 次元ベクトル）

### 4.1 設計原則（v0.0.6 と同一）

```
原則1：ハイパーパラメータ（重み係数）を持たない
原則2：スカラーに集約しない
原則3：各成分は独立した不整合の次元
原則4：加算・重み付け合成を行わない
```

### 4.2 各成分の定義

| 成分 | 定義 | 証拠分類 |
|------|------|---------|
| $T_1$ | $\|\log\mathrm{Vol}(M_A) - \log\mathrm{Vol}(M_B)\|$ 体積スケール差 | Derived |
| $T_2$ | $d_g(p_A, p_B)$ Dijkstra 測地距離 | Derived |
| $T_3$ | 局所 PCA フレーム間の回転角（v0.2 簡易版） | Assumed |
| $T_4$ | $\|\kappa_A - \kappa_B\| / \max(\|\kappa_A\|, \|\kappa_B\|)$ 曲率差 | Assumed |
| $T_5$ | $1 - \sigma_{\min}(J)/\sigma_{\max}(J)$ 埋め込み歪み | Derived |

$T_3, T_4$ は v0.2 では簡易版。v0.2.1 で Singer-Wu vector diffusion / 離散ホロノミーに精密化予定。

```python
tension = calc_regime_tension_vector(embedded_A, embedded_B, k=10)
```

---

## 5. グランドポテンシャル Φ

### 5.1 定義（方向 A：分解構造の保存）

$$\boxed{\Phi = F_{\text{geom}} - \mu \cdot N}$$

| 量 | 定義 | 証拠分類 |
|----|------|---------|
| $F_{\text{geom}}$ | $\lambda_2$（スペクトルギャップ） | Derived |
| $\mu$ | $T_2$（測地距離）= 化学ポテンシャル | Derived |
| $N$ | 移動する人数 | Observed |

### 5.2 放置時に Φ が減少する保証

| 項 | 放置時の変化 | Φ への影響 |
|---|---|---|
| $F_{\text{geom}}$ 低下 | 構造が緩む → $\lambda_2$ 低下 | Φ 減少 ✓ |
| $\mu N$ 増大 | レジーム乖離 → $T_2$ 増大 | Φ 減少 ✓ |

### 5.3 操作変数とその感度

| 操作変数 | $\partial\Phi/\partial u_i$ | 介入の実体 |
|---------|---------------------------|-----------|
| $u_1 = N$ | $-\mu = -T_2$ | 集客・送客・流入制御 |
| $u_2 = T_2$ | $-N$ | レジーム間の距離を縮める |
| $u_3 = \theta$ | $\partial F_{\text{geom}}/\partial\theta$ | 構造変化 |

```python
mu  = calc_mu(tension.T2)          # μ = T₂
Phi = calc_grand_potential(F_geom, mu, N)
```

---

## 6. 位相判定（5 相）

S_g（正規化）・D（正規化）・λ₁ から判定する。

| 位相 | S_norm | D_norm | λ₁ | 意味 |
|------|--------|--------|-----|------|
| **Living System** | ≤ 0.6 | 0.3-0.7 | < 0 | 動的平衡 |
| **Frozen Order** | ≤ 0.3 | ≤ 0.3 | < -0.1 | 硬直 |
| **Runaway Growth** | 中 | ≥ 0.7 | < 0.05 | 過成長 |
| **Noise Dominant** | ≥ 0.7 | ≥ 0.6 | ≥ 0 | 構造崩壊の前兆 |
| **Collapse** | ≥ 0.85 | — | > 0.1 | 持続不能 |

判定順序：Collapse → Noise Dominant → Frozen Order → Living System → Runaway Growth

```python
label, css, explanation = classify_phase(S_g, D, lambda_1, d_embed)
```

### 6.2 相転移予兆検出（3 条件 AND）

```
alert = (dλ₁/dt > 0.02) ∧ (Var(D) > 0.01) ∧ (dS_g/dt > 0.3)
```

---

## 7. 行動条件診断（§4-2）

既存の 𝓣・Φ・λ₁ の読み替えのみで構成する。

| 診断項目 | 対応変数 | 行動しやすい条件 |
|---------|---------|----------------|
| 到達可能性 | $\mu = T_2$ | μ 小（測地距離が短い） |
| 伝達効率 | $T_5$ | T₅ 小（埋め込み歪みが小さい） |
| 方向整合性 | $T_3$ | T₃ 小（フレーム回転が小さい） |
| スケール適合 | $T_1$ | T₁ 小（体積スケールが近い） |
| 変化の好機 | $d\Phi/dt$ | dΦ/dt < 0（Φ が減少中） |
| 文脈の持続性 | $\lambda_1$ | λ₁ < 0（軌道が安定・収束） |

```python
diag = calc_action_condition(T1, T2, T3, T5, lambda_1, Phi, Phi_prev)
```

---

## 8. 介入最適化（§3）

v0.0.6 と同じ定式化。v0.0.6 のバグ 3 件を修正済み。

### 修正済みバグ

| バグ | v0.0.6 | v0.2 修正 |
|------|--------|----------|
| λ スキャン no-op | `priorities = eff × (1+λ)` で argmax 不変 | λ を介入発火閾値の余裕量に反映 |
| ratio 閾値符号反転 | `Φ₀ = Φ × 0.95` が負の Φ で反転 | `Φ₀ = Φ - 0.05·(|Φ|+1)` 符号不変式 |
| dΦ/dt 初期値病的 | `phi_prev = 0` | `phi_prev = Phi`（最初のステップは dΦ/dt = 0） |

### τ の自動設定（λ₁ ベース）

| λ₁ の範囲 | τ | 理由 |
|-----------|---|------|
| < -0.1 | 5 | 安定、長期予測可能 |
| < 0 | 3 | やや安定 |
| < 0.05 | 2 | 中立 |
| ≥ 0.05 | 1 | カオス的、最短で予測 |

---

## 9. 解析プロセス（ステップバイステップ）

### Step 1：Input Audit

```
データ種類・粒度・欠損を確認
レジーム A・B の識別（ユーザー分割）
N（移動人数）の確認
```

### Step 2：Takens 埋め込み

```python
emb_A, d, tau = auto_embed(x_A)
emb_B, _, _   = auto_embed(x_B)   # 同じ d, tau で揃えるのが望ましい
```

### Step 3：Layer 1 構造場変数

```python
S_g      = calc_volume_entropy(emb_A)
D        = calc_correlation_dim(emb_A)
lambda_1 = calc_lyapunov(emb_A)
```

### Step 4：Layer 2 多様体幾何

```python
F_geom = calc_spectral_gap(emb_A, k=10)
kappa  = calc_ollivier_ricci(emb_A, k=10)
```

### Step 5：𝓣 ベクトル

```python
tension = calc_regime_tension_vector(emb_A, emb_B, k=10)
```

### Step 6：Φ

```python
mu  = calc_mu(tension.T2)
Phi = calc_grand_potential(F_geom, mu, N)
```

### Step 7：位相判定 → 行動条件診断 → 介入

```python
label, css, expl = classify_phase(S_g, D, lambda_1, d_embed=d)
diag = calc_action_condition(...)
result = run_intervention_loop(states, phi_0, tau=tau_val, ...)
```

---

## 10. 4 層証拠分類

| ラベル | 定義 | v0.2 の例 |
|--------|------|----------|
| **Observed** | 実データから直接計測 | N（人数）、観測時系列 |
| **Derived** | 観測値から確定的に算出 | S_g, D, λ₁, F_geom, κ, T₁, T₂, T₅, μ, Φ |
| **Assumed** | 分析者が設定した仮定 | T₃, T₄（簡易版）、k（近傍数）、コスト c_i |
| **Hypothesis** | 解釈的・推論的推定 | データなし時の代理推定 |

---

## 11. 出力仕様

- **形式：** HTML 単一ファイル
- **ファイル名：** `sdft-revised-v02_YYYY-MM-DD-HH-MM.html`
- **出力先：** `/mnt/user-data/outputs/`

レポート構造は v0.0.6 と同一（7 セクション）。

---

## 12. 品質ゲート

```
□ §2 の 𝓣 が5成分個別で表示されているか（スカラー集約なし）
□ §2 に Φ の値と dΦ/dt が表示されているか
□ §3 に priority_score テーブルがあるか
□ §4-2 に行動条件診断テーブル（6項目）があるか
□ λ₁ が使用されているか（H は使用禁止）
□ Fisher 行列 G が使用されていないか（v0.3 まで保留）
□ クレジットが v0.2 になっているか
□ テストが全件パスしているか
```

---

## 13. 実装ファイル構成

```
sdft-revised-v02/
├── SKILL.md                            本ドキュメント
├── references/
│   └── theory-manifold.md              多様体理論の前提
├── assets/
│   └── report_template_v02.html        レポートテンプレート
├── scripts/
│   ├── sdft_v02_embedding.py           Takens / S_g / D / λ₁
│   ├── sdft_v02_geometry.py            k-NN / 測地距離 / κ / λ₂
│   ├── sdft_v02_tension.py             𝓣 5成分（多様体版）
│   ├── sdft_v02_potential.py           Φ / 位相判定 / KPI
│   └── sdft_v02_intervention.py        介入最適化（バグ修正済み）
└── tests/
    └── test_core.py                     回帰テスト（12件）
```

---

## 14. 巻末用語集

| 用語 | 説明 |
|------|------|
| **$S_g$（軌道体積エントロピー）** | $k$-NN ball 体積の平均対数。軌道の広がり。 |
| **$D$（相関次元）** | Grassberger-Procaccia 法。アトラクタの次元。 |
| **$\lambda_1$（最大リアプノフ指数）** | 近傍軌道の発散率。$<0$ で安定、$>0$ でカオス。 |
| **$g$（リーマン計量）** | $k$-NN グラフの重み付き隣接行列から暗黙的に定義。 |
| **$\kappa$（Ollivier-Ricci 曲率）** | $1 - W_1(m_x,m_y)/d(x,y)$。離散リッチ曲率。 |
| **$\lambda_2$（スペクトルギャップ）** | Fiedler 値。系の結合の強さ。$= F_{\text{geom}}$。 |
| **$d_g$（測地距離）** | $k$-NN グラフ上の Dijkstra 最短経路長。 |
| **$\Phi$（グランドポテンシャル）** | $F_{\text{geom}} - \mu N$。放置すれば減少。 |
| **$\mu$（化学ポテンシャル）** | $T_2$（測地距離）。1人の移動コスト。 |
| **$T_1$** | 体積スケール差。レジーム間の「大きさ」の違い。 |
| **$T_2$** | 測地距離。レジーム間の多様体上の距離。 |
| **$T_3$** | 平行移動不整合（簡易版：PCA フレーム回転角）。 |
| **$T_4$** | ホロノミー（簡易版：曲率差の正規化）。 |
| **$T_5$** | 埋め込み歪み。Jacobian 特異値比。 |
| **Takens 埋め込み** | 時系列を $d$ 次元遅延座標空間に持ち上げる。 |
| **τ（遅延）** | 埋め込みの時間遅れ。相互情報量極小で自動推定。 |
| **$d$（埋め込み次元）** | FNN 法で自動推定。 |
| **Living System** | 動的平衡。中 S_g・中 D・安定 λ₁。 |
| **Frozen Order** | 硬直。低 S_g・低 D・強安定 λ₁。 |
| **Runaway Growth** | 過成長。中 S_g・高 D。 |
| **Noise Dominant** | 構造崩壊の前兆。高 S_g・高 D・カオス λ₁。 |
| **Collapse** | 持続不能。極高 S_g・カオス λ₁。 |

---

クレジット：株式会社アドインテ SDFT位相空間関係性モデル sdft-revised v0.2
