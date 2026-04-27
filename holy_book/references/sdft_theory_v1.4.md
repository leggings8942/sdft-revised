# SDFT Unified v1.4.2 理論リファレンス

> 完全なHTMLドキュメントは次のパスに格納されています：
> '/Users/arakawa/Library/Application Support/Claude/local-agent-mode-sessions/skills-plugin/37f04196-6007-4008-b0ff-2b60ba7e8dea/dbbe5837-adb7-4108-8bdc-44cd306119c3/skills/sdft-v142-analysy/SDFT_Unified_Complete_v1.4.2.html'
>
> Pythonリファレンス実装：
> '/Users/arakawa/Library/Application Support/Claude/local-agent-mode-sessions/skills-plugin/37f04196-6007-4008-b0ff-2b60ba7e8dea/dbbe5837-adb7-4108-8bdc-44cd306119c3/skills/sdft-v142-analysy/scripts/sdft_unified_v1_4_A_reference.py'

---

## v1.4.2 コア数式

### 状態ベクトル（v1.4.2）

```
X(t)     = (S, D, H)              # 最小構成
X_ext(t) = (S, D, H, Θ, F, A, q) # 拡張（v1.4.2：アフォーダンス・四元数含む）
```

### 各指標の計算式

```python
S = calc_entropy(window, bins=32)       # シャノンエントロピー
D = higuchi_fd(window, k_max=10)        # ヒグチフラクタル次元
H = hurst_exponent(window)              # Hurst指数（R/S近似）
F = U - Theta * S                       # 自由エネルギー（Θ はTemperature-like Parameter）
```

### レジーム間テンション（v1.4.2）

```python
def regime_tension(LA, LB, delta_topo, delta_op,
                   kappa=1.0, lam=1.0,
                   gauge_dist=0.0, eta=1.0):
    ratio = max(LA, EPS) / max(LB, EPS)
    return abs(math.log(ratio)) + kappa * delta_topo + lam * delta_op + eta * gauge_dist

# v1.4.2 追加：エッジ上のレジームテンション（バリア項を含む）
T_reg(edge) = T_reg_base + normative_barrier + cognitive_barrier + max(0, 1 - affordance_alignment)
```

**記号整理（v1.4.2）：**
- Θ（Theta）：Temperature-like parameter（自由エネルギー計算用。𝓣 とは別物）
- 𝓣（Regime Tension）：レジーム間テンション（構造差・演算子差・ゲージ歪み + 規範的・認知的バリアの総和）

### ゲージ不変性・四元数層

```python
# 四元数状態
q(w, x, y, z) -> np.ndarray

# ゲージ変換（ノード）
q'_i = g_i ⊗ q_i

# ゲージ変換（エッジ）
W'_ij = g_j ⊗ W_ij ⊗ g_i^{-1}

# ゲージ歪み（ループホロノミー）
gauge_distortion = ||holonomy(W_cycle) - I||
holonomy = W_1 ⊗ W_2 ⊗ ... ⊗ W_n（ループ上の演算子積）
```

---

## アフォーダンス層（v1.4.2 新規追加）

v1.4.2 の最大の拡張。行為可能性を定量化する第二層フィールド。

### アフォーダンス算出

```python
# アフォーダンスlogit（6成分の線形結合）
A_logit = (w_salience * salience
         + w_feasibility * feasibility
         + w_legitimacy * legitimacy
         + w_reward * reward
         + w_mutuality * mutuality
         - w_cost * cognitive_cost
         + bias)

# アフォーダンススコア（sigmoid変換、0〜1）
A = sigmoid(A_logit)  # sigmoid(x) = 1 / (1 + exp(-x))

# アフォーダンスポテンシャル
phi_A = nu_A * A

# エッジ上のアフォーダンス流
aff_flow = clip(A_src + affordance_gain - normative_barrier - cognitive_barrier - semantic_penalty, 0, 1)

# 行動確率
P_behavior = sigmoid(alpha_A * A + alpha_F * F + alpha_H * H - alpha_T * T_reg + bias)
# デフォルト係数：alpha_A=2.0, alpha_F=1.0, alpha_H=0.5, alpha_T=1.25
```

### アフォーダンス6成分の定義

| フィールド | 意味 |
|-----------|------|
| **salience** | 顕在性：対象・機会がどれだけ知覚されやすいか |
| **feasibility** | 実行可能性：実際に行動できると感じる度合い |
| **legitimacy** | 正当性：社会的・規範的に許容されるか |
| **reward** | 報酬期待：行動から得られると期待する価値 |
| **mutuality** | 相互性：周囲との協調・共同利益感 |
| **cognitive_cost** | 認知的コスト：行動に伴う精神的・認知的負荷（マイナス項） |

### エッジ上の新規フィールド（v1.4.2）

| フィールド | 意味 |
|-----------|------|
| **normative_barrier** | 規範的障壁：制度・慣行・文化的抵抗 |
| **cognitive_barrier** | 認知的障壁：理解困難性・不確実性 |
| **affordance_gain** | エッジ通過によるアフォーダンス増減 |
| **affordance_alignment** | 送受信ノード間のアフォーダンス方向一致度 |

### アフォーダンスアラート条件（v1.4.2）

```python
alert_bad_transition  = (dA < -0.05) & (dH < -0.03) & (D_var > threshold)
alert_good_transition = (dA > 0.05)  & (H.diff() >= 0.0)
```

---

## マクスウェル伝播パラメータ（v1.4.2）

```python
# 媒質パラメータ
epsilon = eps0 * (1 + 1/(1 + T_reg))   # 受容性（スラック）
mu      = mu0 * (1 + T_reg)             # 慣性（制度的質量）
sigma   = sigma0 + a_s * S + b_t * T_reg  # 散逸（エントロピー＋サイロ壁）

# 伝播特性
v     = 1 / sqrt(epsilon * mu)          # 伝播速度
Z     = sqrt(mu / epsilon)              # インピーダンス
alpha = 0.5 * sigma * sqrt(mu / epsilon)  # 減衰係数
tau   = 4*Z1*Z2 / (Z1+Z2)^2            # 透過率（境界面）

# v1.4.2 追加：電場のアフォーダンス合算
E_struct = -(phi_j - phi_i) / length - dA_dt   # 構造電場
E_aff    = -(phi_A_j - phi_A_i) / length        # アフォーダンス電場
E_total  = E_struct + gamma_A * E_aff            # 合算電場（v1.4.2）
```

---

## ハミルトニアン・ローレンツ力（v1.4.2）

```python
# v1.4.2 有効ポテンシャル（構造 + アフォーダンス）
phi_eff = phi_struct + gamma_A * phi_A

# v1.4.2 アフォーダンスペナルティ項
aff_penalty = beta_C * cognitive_cost^2 + beta_L * (1 - legitimacy)^2

# v1.4.2 有効ハミルトニアン（アフォーダンス統合版）
H_eff = (p - Q*A_vec)^2 / (2*mass) + Q*phi_eff + alpha_D * D^2 + beta_T * T_reg^2 + aff_penalty

# ローレンツ力（行動変容シミュレーション）
F_Lorentz = Q * (E_total + v × B)
# E_total: 合算電場（構造電場 + アフォーダンス電場、v1.4.2）
# B: 磁場（ループホロノミーから）
# v × B: サイドフォース（側面的変容圧力）
```

---

## 位相ラベルと相転移条件

| Phase | S | D | H | A（v1.4.2補完） | 意味 |
|-------|---|---|---|----------------|------|
| Frozen Order | Low | Low | High | 高くても行動起きにくい | 安定だが変化乏しい |
| Living System | Mid | Mid | Mid-High | 中〜高：健全な動的平衡 | 動的平衡・自己修復 |
| Runaway Growth | Mid | High | High | 高：過剰・バブル注意 | 強い増幅相 |
| Noise Dominant | High | High | Low-Mid | 低下：認知的阻害が加速 | ノイズ優勢 |
| Collapse | High | Unstable | Low | 極低：行動不能 | 持続不能 |

**相転移前兆（v1.4.2 最重要条件）：**
```
(dH/dt < 0) かつ (Var(D) 増大) かつ (dA/dt < 0)  →  最強の相転移前兆

# 数値閾値付き版（スクリプト実装）：
(dH < -0.03) かつ (D_var > Q75) かつ (dS > 0)  →  alert = True
```

---

## 介入優先順位スコア（v1.4.2）

```python
priority_score = (
    (1 - A) * 2.0              # アフォーダンス不足（行為可能性の低さ）
    + incoming_T_reg * 0.7     # 受信テンション
    + cognitive_cost * 1.0     # 認知的コスト
    + (1 - legitimacy) * 1.2   # 正当性の欠如
)
# スコアが大きいノード・エッジへの介入を優先する
```

---

## 知識優先順位（v1.4.2）

1. `SDFT_Unified_Complete_v1.4.2.html` ← **最優先（カノニカルソース）**
2. `sdft_unified_v1_4_A_reference.py` ← 参照実装（v1.4.2 実装正本・アフォーダンス拡張統合版）
3. `SDFT_Unified_Complete_v1.4.1.html`（v1.4.1）← フォールバック
4. 旧版はフォールバックのみ

**競合ルール：**
- Θ（温度様パラメータ）と 𝓣（レジーム間テンション）の表記を混同しないこと
- A（アフォーダンス）を S（エントロピー）の代替として扱わないこと（行為可能性の低さ ≠ 無秩序性の高さ）
- ハミルトニアン式中の `A_vec`（ベクトルポテンシャル）と A（アフォーダンス）を混同しないこと
