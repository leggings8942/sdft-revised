---
name: sdft-revised
description: |
  SDFT revised v0.0.4 構造場分析エージェント。
  構造場変数（S/D/H/F）・レジーム間テンション𝓣（5次元ベクトル）・
  グランドポテンシャルΦによる介入最適化設計方針をHTMLレポートで生成します。

  v0.0.4 の設計思想：
    - 𝓣 はハイパーパラメータを持たない5次元ベクトル (T₁,T₂,T₃,T₄,T₅)
    - グランドポテンシャル Φ = F − μN で介入優先順位を設計
    - 放置すれば Φ は減少する → 最小インパクトで Φ ≥ Φ₀ を維持
    - 全変数を情報幾何学・熱力学的アナロジーで定義
    - 全変数に証拠分類ラベル（Observed/Derived/Assumed/Hypothesis）を付与

  対応領域：マーケット・組織・施設・都市・複合施設など複合領域の構造診断。
  解析順序：Input Audit → SDFT Mapping → 変数計算 → 位相判定 → Φ計算 → レポート生成。

  以下を含む場合は必ずこのスキルを使用：sdft-revised、SDFT構造場分析、
  エントロピー分析、Hurst指数、フラクタル次元、相転移予兆、位相空間、
  構造場診断、レジーム間テンション、グランドポテンシャル、情報幾何学、
  Bhattacharyya距離、αダイバージェンス、統計的ホロノミー、通信路損失、
  化学ポテンシャル、介入最適化。
---

# SDFT revised v0.0.4 構造場分析エージェント

---

## バージョンと対応範囲

**このスキルは v0.0.4 です。**

| 機能 | 状態 | 備考 |
|------|------|------|
| 構造場分析（S/D/H） | ✅ 対応 | 時系列から自動計算 |
| 自由エネルギー F | ✅ 再定義 | Fisher 行列から Derived |
| 𝓣（5次元ベクトル） | ✅ 対応 | ハイパーパラメータなし |
| グランドポテンシャル Φ | ✅ 新規 | Fisher 行列・T₂・N から定義 |
| 介入最適化（設計方針） | ✅ 新規 | ∂Φ/∂u_i に基づく priority_score |
| 位相判定（5相） | ✅ 対応 | S/D/H から判定 |
| 証拠分類（4層） | ✅ 対応 | 全変数に付与 |
| holy_book（v1.4.2原典） | ✅ 動的参照可能 | 比較・補完目的のみ |
| アフォーダンス層（A/P） | ❌ 未対応 | v0.1以降 |
| 介入最適化（Python完全実装） | ❌ 未完全実装 | v0.1以降 |
| ゲージ/Maxwell/Hamiltonian | ❌ 削除済み | 理論的整合性なし |

---

## カノニカルソース（最優先）

```
{スキルディレクトリ}/references/SDFT_Unified_Complete_v0.0.4.html
```

---

## holy_book — 原典（v1.4.2）アーカイブ

```
{スキルディレクトリ}/holy_book/
```

参照専用。分析には使用しない。

| 参照すべき状況 | 参照先ファイル |
|--------------|--------------|
| 𝓣 の旧定義との比較 | `holy_book/SDFT_Unified_Complete_v1.4.2.html` |
| v1.4.2 の変数定義を引用 | `holy_book/references/sdft_theory_v1.4.md` |
| 削除済み変数の元の意図 | `holy_book/SDFT_Unified_Complete_v1.4.2.html` |
| 旧実装アルゴリズム参照 | `holy_book/scripts/sdft_unified_v1_4_A_reference.py` |

**⚠️ holy_book の式・変数・実装を分析に直接使用してはならない。**

---

## 削除済み変数（使用禁止）

```
q, W_ij, ε, μ_maxwell, σ, v, Z, α_maxwell, τ
E, B, φ_struct, φ_A, E_total, H_eff, F_lorentz
charge, mass, LA, LB, ΔTopo, ΔOp
A, P_behavior, aff_flow
salience, feasibility, legitimacy, reward, mutuality, cognitive_cost
normative_barrier, cognitive_barrier, affordance_gain, affordance_alignment
```

---

# 各指標の詳細解説

## 1. S（シャノンエントロピー → 熱力学的エントロピー）

$$H_{bit} = -\sum_i p_i \log_2 p_i \quad [\text{bit}]$$

$$S = k_B \ln 2 \cdot H_{bit} \quad [J/K], \quad k_B = 1.380649 \times 10^{-23} \text{ J/K}$$

**意味：** 系の「無秩序さ・散逸度」。S が高いほど無秩序で散逸が進んでいる。

**計算関数：**
```python
H_bit = calc_H_bit(x, bins=32)   # [bit]
S     = calc_S(x, bins=32)        # [J/K]
```

**証拠分類：** Derived

---

## 2. D（フラクタル次元）

$$D = \lim_{k \to 0} \frac{\log L(k)}{\log(1/k)} \in [1.0,\ 2.0]$$

**意味：** D→1.0 で単純（規則的）、D→2.0 で完全ランダム（不規則）。時系列の「複雑さの構造」を測る。

**計算関数：**
```python
D = calc_D(x, k_max=10)
```

**証拠分類：** Derived

---

## 3. H（Hurst 指数）

$$H \in [0.0,\ 1.0]$$

**意味：** H>0.5 で持続性（トレンドが続く）、H<0.5 で反持続性（反転しやすい）。**相転移の最も重要な先行指標**（崩壊の前に H が低下する）。

**計算関数：**
```python
H = calc_H(x, max_lag=20)
```

**証拠分類：** Derived

---

## 4. Fisher 情報行列 G

$$G_{ij}(\theta) = \mathbb{E}\!\left[\frac{\partial \log p}{\partial \theta_i}\frac{\partial \log p}{\partial \theta_j}\right]$$

正規分布 $\mathcal{N}(\mu, \sigma^2)$ の閉形式：

$$G = \begin{pmatrix} 1/\sigma^2 & 0 \\ 0 & 2/\sigma^4 \end{pmatrix}$$

**意味：** 統計多様体上の計量テンソル。G が大きい = 感度が高い = エネルギーランドスケープが急峻。Jaynes の最大エントロピー原理により、G は熱力学的温度と直接対応する。

**計算関数：**
```python
G = calc_fisher_matrix(x, d=2)   # 正規分布仮定・2×2
```

**証拠分類：** Derived（正規分布仮定のもとで）

---

## 5. T（情報的温度）

$$T = \frac{\mathrm{tr}(G^{-1}(\theta))}{d \cdot k_B} \quad [K]$$

**理論的背景：** Jaynes（1957）により、$k_BT \longleftrightarrow \mathrm{tr}(G^{-1})/d$ の対応が厳密に導かれる（系の揺らぎのエネルギースケール）。T が大きい = 散逸が進んでいる。

**計算関数：**
```python
T = calc_T_temperature(G)
```

**証拠分類：** Derived

---

## 6. U（内部エネルギー）

$$U = \frac{\mathrm{tr}(G(\theta))}{d} \quad [J \text{ 相当}]$$

**意味：** Fisher 行列のトレース = 情報的エネルギー密度。等分配則との対応から導かれる。

**計算関数：**
```python
U = calc_U_internal_energy(G)
```

**証拠分類：** Derived

---

## 7. F（自由エネルギー）

$$F = U - T \cdot S = \frac{\mathrm{tr}(G)}{d} - \frac{\ln 2 \cdot H_{bit} \cdot \mathrm{tr}(G^{-1})}{d} \quad [J]$$

**意味：** 系が外部に対して仕事をできる最大量。F<0 は外部から仕事を与えないと維持できない状態。

**v0.0.4 の重要な変更：** 旧来の $F = U - \Theta S$（Assumed）から Fisher 行列由来の Derived に変更。

**計算関数：**
```python
F = calc_F_free_energy(G, H_bit)
```

**証拠分類：** Derived

---

## 8. 𝓣（レジーム間テンション・5次元ベクトル）

$$\vec{\mathcal{T}}(A, B) = (T_1,\ T_2,\ T_3,\ T_4,\ T_5) \in \mathbb{R}^5$$

**設計原則：** ハイパーパラメータなし・スカラー集約なし・各成分は独立した次元。

### 原典（v1.4.2）との対比

| 項 | v1.4.2 | v0.0.4 | 理論基盤の変化 |
|----|--------|--------|--------------|
| T₁ | $\|\log(L_A/L_B)\|$ | $\|\psi_A-\psi_B\|$ | ゲージ → 対数分配関数 |
| T₂ | $\kappa\cdot\Delta\text{Topo}$ | $d_B(p_A,p_B)$ | 位相幾何 → Bhattacharyya |
| T₃ | $\lambda\cdot\Delta\text{Op}$ | $D^{(\alpha\neq0)}(p_A\|p_B)$ | 作用素論 → α-ダイバージェンス |
| T₄ | $\eta\cdot\|hol-I\|$ | $\|\text{Hol}^{(\alpha)}-I\|$ | 四元数 → 統計的ホロノミー |
| T₅ | $b_{norm}+b_{cog}+\max(...)$ | $\mathcal{L}(e)$ | 社会学 → 通信路損失 |

---

### T₁：対数分配関数差

$$T_1 = |\psi_A - \psi_B|, \quad \psi(\theta) = \log Z(\theta)$$

**意味：** 2つのレジームのスケール・自由度の乖離。**証拠分類：** Derived

```python
T1, ev1 = calc_T1(data_A, data_B)
```

---

### T₂：Bhattacharyya 距離

$$T_2 = d_B(p_A, p_B) = -\log \int\sqrt{p_A(x)\,p_B(x)}\,dx$$

**意味：** 統計多様体上の2点間の幾何学的距離（Fisher-Rao 近似）。化学ポテンシャル $\mu = 2\sqrt{T_2}$ に直結。⚠️ T₃ との独立性のため α≠0 を守ること。**証拠分類：** Derived

```python
T2, ev2 = calc_T2(data_A, data_B)
```

---

### T₃：α-ダイバージェンス（α ≠ 0）

$$T_3 = D^{(\alpha)}(p_A\|p_B), \quad \alpha \neq 0$$

**意味：** 意味変換の非対称な差。α=+1 で A の視点から B の「意外さ」、α=−1 で逆。⚠️ α=0 は T₂ と重複するため禁止。**証拠分類：** Derived

```python
T3, ev3 = calc_T3(data_A, data_B, alpha=1.0)   # α≠0 必須
```

---

### T₄：統計的ホロノミー

$$T_4 = \|\mathrm{Hol}^{(\alpha)}(\gamma) - I\|$$

**意味：** 複数レジームを巡る閉ループで情報幾何学的整合性が保たれるか。**証拠分類：** Assumed（近似推定）

```python
T4, ev4 = calc_T4(data_A, data_B, alpha=1.0)
```

---

### T₅：通信路損失率

$$T_5 = \mathcal{L}(e) = 1 - \frac{I(X_A;\,X_B)}{H(X_A)}$$

**意味：** エッジ通過時の情報損失率（0=完全伝達、1=全損）。旧 $b_{norm}+b_{cog}+\max(...)$ の統合。**証拠分類：** Derived（同時観測データがある場合）/ Assumed（推定の場合）

```python
T5, ev5 = calc_T5(data_A, data_B)
```

---

### 𝓣 の統合計算

```python
tension = calc_regime_tension_vector(data_A, data_B, alpha_T3=1.0)
# tension.T1, T2, T3, T4, T5 を個別に使用
# スカラー集約禁止
```

---

## 9. グランドポテンシャル Φ

$$\Phi = F - \mu N = \frac{\mathrm{tr}(G)}{d} - \frac{\ln 2\cdot H_{bit}\cdot\mathrm{tr}(G^{-1})}{d} - 2\sqrt{T_2}\cdot N \quad [J]$$

$$d\Phi = -S\,dT - P\,dV - N\,d\mu$$

**各変数の対応：**

| 変数 | 定義 | 次元 | 証拠分類 |
|------|------|------|---------|
| $S$ | $k_B\ln 2\cdot H_{bit}$ | J/K | Derived |
| $T$ | $\mathrm{tr}(G^{-1})/(d\cdot k_B)$ | K | Derived |
| $P$ | $\mathrm{sign}(d\Phi)\cdot\|d\Phi\|_G$ | J | Derived |
| $V$，$dV$ | $\int\sqrt{\det G}\,d^d\theta$，$\sqrt{\det G}\,d^d\theta$ | m³相当 | Derived |
| $N$ | 移動する人数 | 人 | Observed |
| $\mu = 2\sqrt{T_2}$ | 1人の移動コスト | J/人 | Derived |

**放置すれば Φ が減少する理由：**

| 項 | 放置時 | 符号 | Φ への影響 |
|---|---|---|---|
| $-S\,dT$ | 散逸→$T\uparrow$（$dT>0$） | $-(+)(+)<0$ | ✅ 減少 |
| $-P\,dV$ | $P<0$，$dV<0$ | $-(-)(-)<0$ | ✅ 減少 |
| $-N\,d\mu$ | $T_2\uparrow$→$\mu\uparrow$（$d\mu>0$） | $-(+)(+)<0$ | ✅ 減少 |

**計算関数：**
```python
state = analyze_single_regime(x, N=N_people, T2=T2_value)
# state.Phi がグランドポテンシャル [J]
```

**証拠分類：** Derived（F が Derived かつ N が Observed の場合）

---

# 解析プロセス（ステップバイステップ）

## Step 1：Input Audit

```
・データ種類・粒度・欠損を確認
・各データに証拠分類ラベルを付ける
・レジームA・B の識別
・N（レジーム間を移動する人数）の確認
・𝓣 計算に使えるデータ量の確認
```

データがない場合の代理推定（全て Hypothesis ラベル）：
```
S → 行動多様性・需要変動係数から推定
D → 導線複雑度・競合多層性から推定
H → トレンド継続期間から推定
```

---

## Step 2：SDFT Mapping

```
ノード   = 構成要素（部門・ゾーン・セグメント）
リンク   = 接続（強度・方向のみ）
場       = 外部環境（市場・競合・マクロ）
境界     = レジームの境界
レジーム = 同一オペモードを持つ領域
```

---

## Step 3：構造場変数の計算

```python
H_bit = calc_H_bit(x)     # [bit]    Derived
S     = calc_S(x)          # [J/K]   Derived
D     = calc_D(x)          # [1-2]   Derived
H_exp = calc_H(x)          # [0-1]   Derived

G = calc_fisher_matrix(x)
T = calc_T_temperature(G)              # [K]     Derived
U = calc_U_internal_energy(G)         # [J相当] Derived
F = calc_F_free_energy(G, H_bit)      # [J]     Derived
V = calc_V_volume(G)                   # [m³相当] Derived
```

**優先解釈順序：**
1. **dH/dt**（最重要先行予兆。崩壊の前に H が低下する）
2. **Var(D)**（構造不安定化の兆候）
3. **dS/dt**（散逸の進行。遅行指標）
4. **Φ の変化方向**（介入が必要かどうかの判断）

---

## Step 4：𝓣 ベクトルの計算

```python
tension = calc_regime_tension_vector(
    data_A=regime_A,
    data_B=regime_B,
    alpha_T3=1.0    # α≠0 必須
)
# T₁〜T₅ を個別に記録・表示（スカラー集約禁止）
```

---

## Step 5：グランドポテンシャル Φ の計算

```python
mu  = calc_mu(tension.T2)               # μ = 2√T₂  [J/人]
Phi = calc_grand_potential(F, mu, N)    # Φ = F - μN [J]
```

---

## Step 6：位相判定

```python
label, css, explanation = classify_phase(S_val, D, H_exp)
```

| 位相 | S | D | H | 意味 |
|------|---|---|---|------|
| Living System | ≤0.6 | 1.2〜1.7 | ≥0.45 | 動的平衡・健全 |
| Frozen Order | ≤0.4 | ≤1.4 | ≥0.55 | 安定・変化なし |
| Runaway Growth | 中 | ≥1.6 | ≥0.6 | 過成長・バブル |
| Noise Dominant | ≥0.65 | ≥1.5 | ≤0.45 | 構造崩壊の前兆 |
| Collapse | ≥0.8 | 不安定 | ≤0.3 | 持続不能 |

**相転移予兆（3条件が全て揃った場合のみ）：**

$$\text{alert} = \left(\frac{dH}{dt}<0\right) \land \left(\mathrm{Var}(D)\uparrow\right) \land \left(\frac{dS}{dt}>0\right)$$

---

## Step 7：介入最適化の設計方針（Φ に基づく）

### 目標

$$\min_{u_1,...,u_K,K} \sum_{k=1}^{K} \|u_k\|^2 + \lambda K \quad \text{s.t.} \quad \Phi(t) \geq \Phi_0 \quad \forall t$$

「最小のインパクト・最小の介入回数で $\Phi(t) \geq \Phi_0$ を維持・増加させる」

### 操作変数と感度

| 操作変数 | $\partial\Phi/\partial u_i$ | 優先度が高い条件 |
|---------|---------------------------|----------------|
| $u_1=N$（人数） | $-\mu = -2\sqrt{T_2}$ | T₂ が大きいとき |
| $u_2=T_2$（Bhattacharyya距離） | $-N/\sqrt{T_2}$ | N が大きいとき |
| $u_3=\theta$（状態パラメータ） | $-S\partial T/\partial\theta - P\partial V/\partial\theta$ | 構造変化が必要なとき |

### アルゴリズムの手順

```python
result = calc_intervention_priority(
    Phi_current=Phi_t,
    Phi_prev=Phi_prev,
    dt=1.0,
    tau=5.0,               # 予測ホライズン
    Phi_threshold=Phi_0,   # 維持すべき下限
    T2=tension.T2,
    N=N_people,
    G=G,
    H_bit=H_bit,
)
# result.best_variable  → 最優先操作変数
# result.delta_u_optimal → 必要操作量
# result.priority_N/T2/theta → 各優先スコア
```

```
Step 1: Φ(t) を計算
Step 2: dΦ/dt ≈ (Φ(t) - Φ(t-Δt))/Δt を推定
Step 3: Φ̂(t+τ) = Φ(t) + dΦ/dt·τ を予測
Step 4: Φ̂(t+τ) < Φ₀ なら介入、そうでなければスキップ
Step 5: ΔΦ_need = Φ₀ - Φ̂(t+τ) を計算
Step 6: priority_i = |∂Φ/∂u_i| / c_i を計算
Step 7: best u* = argmax_i priority_i を選択
Step 8: Δu* = ΔΦ_need / (∂Φ/∂u*) で最小インパクトを計算
Step 9: 介入実行
```

⚠️ **介入最適化の Python 完全実装は v0.1 以降。** v0.0.4 では設計方針と priority_score の計算のみ。

---

## 基本動作原則

1. **カノニカルソースを参照する：** `SDFT_Unified_Complete_v0.0.4.html` が最優先
2. **4層証拠分類を全変数に付与する：**
   - **Observed**：実データから直接計測（N など）
   - **Derived**：観測値から確定的に算出（S/D/H/T/U/F/Φ/T₁/T₂/T₃/T₅ など）
   - **Assumed**：分析者が設定した仮定（T₄ の近似、介入コスト $c_i$ など）
   - **Hypothesis**：解釈的・推論的推定（データ代理推定値など）
3. **削除済み変数を使用しない**
4. **𝓣 をスカラー化しない**（5成分を常に個別に保持・表示）
5. **F の計算に旧来の $F = U - \Theta S$（Assumed 版）を使わない**
6. **エグゼクティブファーストを維持する**（§1 は平易な言語で書く）
7. **holy_book は参照専用**（分析には使用しない）

---

## 出力仕様

- **形式：** HTML 単一ファイル
- **ファイル名：** `sdft-revised-v004_YYYY-MM-DD-HH-MM.html`
- **テンプレート：** `{スキルディレクトリ}/assets/report_template_v004.html`
- **出力先：** `/mnt/user-data/outputs/`
- **生成方法：** テンプレート読込 → プレースホルダー置換 Python スクリプト生成 → 実行 → 保存

HTML をレスポンスに直接出力することは禁止。

### レポート構造（7セクション）

```
§1. エグゼクティブサマリー
    技術用語禁止・状況/問題/アクション3件/期待効果/トップリスク

§2. KPI ダッシュボード
    フェーズバッジ（5相）+ 相転移アラート
    KPI カード：S / D / H / F / Φ
    𝓣 の5成分個別カード（T₁〜T₅・スカラー集約なし）

§3. 施策提案・介入優先順位
    Φ に基づく priority_score を表示
    priority_i = |∂Φ/∂u_i| / c_i
    ⚠️ Python 完全実装は v0.1 以降（実装待ちを明記）

§4. 変革実行環境の診断
    4-1. SDFT 写像
    4-2. アフォーダンス層  ← ⚠️ v0.0.4 未対応
    4-3. 障壁ヒートマップ  ← ⚠️ v0.0.4 未対応
    4-4. 証拠分類表（全変数：S/D/H/T/U/F/Φ/μ/N/T₁〜T₅）

§5. フェーズ判定と環境変化リスク
    S-D 位相空間図（SVG）/ H 時系列グラフ（SVG）/ シナリオ表

§6. 総論
    不確実性の高い領域 / 追加で収集すべき観測データ

§7. 巻末用語集 + クレジット（v0.0.4）
```

### §2 の KPI 表示仕様

```
S:  [値] [J/K]  [Derived]  状態
D:  [値] [1-2]  [Derived]  状態
H:  [値] [0-1]  [Derived]  状態
F:  [値] [J]    [Derived]  状態
Φ:  [値] [J]    [Derived]  dΦ/dt=[値]

𝓣（5成分個別・集約なし）:
  T₁=[値] [Derived]          スケール差
  T₂=[値] [Derived]          幾何学的距離
  T₃=[値] [Derived] α=[値]   非対称差
  T₄=[値] [Assumed]           曲率差
  T₅=[値] [Derived]          通信路損失
```

### §3 の施策提案表示仕様

```
【グランドポテンシャルΦに基づく介入設計方針】

現在 Φ=[値]  目標 Φ₀=[値]
予測（τ後）Φ̂=[値]  → 介入[必要/不要]
ΔΦ_need=[値]

操作変数別 priority_score:
  u₁=N   priority=[値]  感度 ∂Φ/∂N  = -2√T₂ = [値]
  u₂=T₂  priority=[値]  感度 ∂Φ/∂T₂ = -N/√T₂ = [値]
  u₃=θ   priority=[値]  感度 ∂Φ/∂θ  ≈ [値]

推奨介入変数: [u_i]
必要操作量: Δu*=[値]

⚠️ Python完全実装はv0.1以降。上記は設計方針の表示です。
```

---

## 品質ゲート（出力前チェックリスト）

```
□ §1 が平易な言語で書かれているか（技術用語なし）
□ §2 の 𝓣 が5成分個別で表示されているか
□ §2 の 𝓣 各成分に証拠分類が付いているか
□ §2 の 𝓣 にスカラー集約値がないか
□ §2 に Φ の値と dΦ/dt が表示されているか
□ T₃ に使用した α 値が明示されているか（α=0 禁止）
□ §3 に priority_score が ∂Φ/∂u_i に基づいて表示されているか
□ §3 に「v0.1 以降で完全実装」の旨が明記されているか
□ §4-2 / §4-3 が「v0.0.4 未対応」と表示されているか
□ §4-4 に全変数の証拠分類が含まれているか
□ §5 に S-D 位相図と H 時系列が含まれているか
□ F の計算に旧来の F=U-ΘS（Assumed版）を使っていないか
□ 削除済み変数が出力に含まれていないか
□ 𝓣 の加算・スカラー集約が行われていないか
□ 巻末用語集が含まれているか
□ クレジットが v0.0.4 になっているか
```

### ハード制約（絶対禁止）

- 削除済み変数（q/W/ε/μ_maxwell/H_eff 等）を使用する
- A / P_behavior を計算・表示する
- 𝓣 にハイパーパラメータ（κ/λ/η）を掛ける
- 𝓣 をスカラーに集約する
- T₃ の α に 0 を使う
- F = U − ΘS（旧来の Assumed 版）を使用する
- §3 の介入設計を定式化なしで自由に埋める
- Φ を計算せずに介入優先度を出力する
- HTML をレスポンスに直接出力する

---

## 巻末用語集（v0.0.4）

| 用語 | 説明 |
|------|------|
| **$H_{bit}$** | $-\sum p_i\log_2 p_i$ [bit]。シャノンエントロピー。 |
| **$S$** | $k_B\ln 2\cdot H_{bit}$ [J/K]。熱力学的エントロピー。散逸度。 |
| **$D$** | ヒグチ法。1.0〜2.0。時系列の複雑性。 |
| **$H$** | 0〜1。Hurst指数。持続性。相転移の先行指標。 |
| **$G$** | Fisher情報行列。統計多様体上の計量テンソル。 |
| **$T$** | $\mathrm{tr}(G^{-1})/(d\cdot k_B)$ [K]。情報的温度。 |
| **$U$** | $\mathrm{tr}(G)/d$ [J相当]。内部エネルギー。 |
| **$F$** | $U-T\cdot S$ [J]。自由エネルギー。系が使えるエネルギー。 |
| **$\Phi$** | $F-\mu N$ [J]。グランドポテンシャル。放置すれば減少する。 |
| **$P$** | $\mathrm{sign}(d\Phi)\cdot\|d\Phi\|_G$ [J]。情報的圧力（符号付き）。 |
| **$V$** | $\int\sqrt{\det G}\,d^d\theta$。統計多様体の体積。 |
| **$N$** | 移動する人数 [人]。Observed。 |
| **$\mu$** | $2\sqrt{T_2}$ [J/人]。化学ポテンシャル。1人の移動コスト。 |
| **$\vec{\mathcal{T}}$** | $(T_1,...,T_5)$。レジーム間テンション。5次元ベクトル。 |
| **$T_1$** | 対数分配関数差。スケール・自由度の乖離。 |
| **$T_2$** | Bhattacharyya距離。状態点の幾何学的ずれ。 |
| **$T_3$** | α-ダイバージェンス（α≠0）。非対称な意味的差。 |
| **$T_4$** | 統計的ホロノミー。多様体の曲率差。 |
| **$T_5$** | 通信路損失率 $\mathcal{L}(e)$。情報損失割合。 |
| **Living System** | 中S・中D・中〜高H。動的平衡・健全。 |
| **Frozen Order** | 低S・低D・高H。安定だが変化なし。 |
| **Runaway Growth** | 中S・高D・高H。過成長の兆候。 |
| **Noise Dominant** | 高S・高D・低H。構造崩壊の前兆。 |
| **Collapse** | 高S・不安定D・低H。持続不能。 |
| **Observed** | 実測値 |
| **Derived** | 観測値から確定的に算出した値 |
| **Assumed** | 分析者の仮定値 |
| **Hypothesis** | 解釈的・推論的推定 |
| **Jeffreys事前分布密度** | $\sqrt{\det G}$。パラメータ変換に不変な事前分布の密度。 |
| **介入効率 $\eta_i$** | $\|\partial\Phi/\partial u_i\|/c_i$。単位コストあたりの Φ 変化量。 |

クレジット：株式会社アドインテ SDFT位相空間関係性モデル sdft-revised v0.0.4
