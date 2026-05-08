---
name: sdft-revised
description: |
  SDFT revised v0.0.5 構造場分析エージェント。
  構造場変数（S/D/H/F）・レジーム間テンション𝓣（5次元ベクトル）・
  グランドポテンシャルΦによる介入最適化・行動条件診断をHTMLレポートで生成します。

  v0.0.5 の設計思想：
    - アフォーダンス（A / P_behavior）を完全廃止
    - §4-2 を「行動条件診断」として既存変数の読み替えで再構成
    - 新しい変数は一切追加しない
    - 𝓣・Φ・H の読み替えで行動条件を診断する
    - 全変数が Derived または Observed で計算可能

  対応領域：マーケット・組織・施設・都市・複合施設など複合領域の構造診断。
  解析順序：Input Audit → SDFT Mapping → 変数計算 → 位相判定 → Φ計算
            → 行動条件診断 → 介入最適化 → レポート生成。

  以下を含む場合は必ずこのスキルを使用：sdft-revised、SDFT構造場分析、
  エントロピー分析、Hurst指数、フラクタル次元、相転移予兆、位相空間、
  構造場診断、レジーム間テンション、グランドポテンシャル、情報幾何学、
  Bhattacharyya距離、αダイバージェンス、統計的ホロノミー、通信路損失、
  化学ポテンシャル、介入最適化、行動条件診断。
---

# SDFT revised v0.0.5 構造場分析エージェント

---

## バージョンと対応範囲

**このスキルは v0.0.5 です。**

| 機能 | 状態 | v0.0.4 からの変更 |
|------|------|-----------------|
| 構造場分析（S/D/H） | ✅ 対応 | 変更なし |
| 自由エネルギー F | ✅ 対応 | 変更なし |
| 𝓣（5次元ベクトル） | ✅ 対応 | 変更なし |
| グランドポテンシャル Φ | ✅ 対応 | 変更なし |
| 介入最適化ループ | ✅ 対応 | 変更なし |
| **行動条件診断（§4-2）** | ✅ **新規** | **アフォーダンスを既存変数の読み替えで再構成** |
| 位相判定（5相） | ✅ 対応 | 変更なし |
| 証拠分類（4層） | ✅ 対応 | 変更なし |
| holy_book（v1.4.2原典） | ✅ 動的参照可能 | 変更なし |
| **アフォーダンス（A / P_behavior）** | ❌ **廃止** | **v0.0.5 で削除** |
| **§4-3 エッジ障壁ヒートマップ** | ✅ **新規実装** | T₂・T₅ のエッジ単位集計 |
| ゲージ/Maxwell/Hamiltonian | ❌ 削除済み | 変更なし |

---

## カノニカルソース（最優先）

```
{スキルディレクトリ}/references/SDFT_Unified_Complete_v0.0.5.html
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

---

## 削除済み変数（使用禁止）

```
【アフォーダンス系（v0.0.5 で廃止）】
A, P_behavior, aff_flow
salience, feasibility, legitimacy, reward, mutuality, cognitive_cost
affordance_gain, affordance_alignment

【物理アナロジー系（v0.0.1 以降削除継続）】
q, W_ij, ε, μ_maxwell, σ, v, Z, α_maxwell, τ
E, B, φ_struct, H_eff, F_lorentz, charge, mass
LA, LB, ΔTopo, ΔOp
normative_barrier, cognitive_barrier
```

---

# 各指標の詳細解説

## 1. S（シャノンエントロピー → 熱力学的エントロピー）

$$S = k_B \ln 2 \cdot H_{bit}, \quad H_{bit} = -\sum_i p_i\log_2 p_i \quad [\text{bit}]$$

**意味：** 系の散逸度。S が高いほど無秩序で散逸が進んでいる。

```python
H_bit = calc_H_bit(x)   # [bit]
S     = calc_S(x)        # [J/K]
```

**証拠分類：** Derived

---

## 2. D（フラクタル次元）

$$D \in [1.0, 2.0] \quad \text{（ヒグチ法）}$$

**意味：** D→1.0 で規則的、D→2.0 で完全ランダム。時系列の複雑性。

```python
D = calc_D(x)
```

**証拠分類：** Derived

---

## 3. H（Hurst 指数）

$$H \in [0.0, 1.0]$$

**意味：** H>0.5 で持続傾向（トレンドが続く）、H<0.5 で反持続（反転しやすい）。
**相転移の最も重要な先行指標**（崩壊の前に H が低下する）。

```python
H = calc_H(x)
```

**証拠分類：** Derived

---

## 4. Fisher 情報行列 G

$$G_{ij}(\theta) = \mathbb{E}\!\left[\frac{\partial\log p}{\partial\theta_i}\frac{\partial\log p}{\partial\theta_j}\right]$$

**意味：** 統計多様体上の計量テンソル。T/U/F/V/Φ の計算の基盤。

```python
G = calc_fisher_matrix(x, d=2)   # 正規分布仮定・2×2
```

**証拠分類：** Derived（正規分布仮定のもとで）

---

## 5. T（情報的温度）

$$T = \frac{\mathrm{tr}(G^{-1})}{d \cdot k_B} \quad [K]$$

**意味：** 系の揺らぎのエネルギースケール。Jaynes の最大エントロピー原理から $k_BT \leftrightarrow \mathrm{tr}(G^{-1})/d$ の対応が導かれる。

**証拠分類：** Derived

---

## 6. U（内部エネルギー）と F（自由エネルギー）

$$U = \frac{\mathrm{tr}(G)}{d} \quad [J\text{相当}]$$

$$F = U - T \cdot S = \frac{\mathrm{tr}(G)}{d} - \frac{\ln 2 \cdot H_{bit} \cdot \mathrm{tr}(G^{-1})}{d} \quad [J]$$

**意味：** F は「系が外部に仕事をできる最大量」。v0.0.5 では旧来の $F=U-\Theta S$（Assumed）を使わない。

```python
U = calc_U_internal_energy(G)
F = calc_F_free_energy(G, H_bit)
```

**証拠分類：** Derived

---

## 7. 𝓣（レジーム間テンション・5次元ベクトル）

$$\vec{\mathcal{T}}(A,B) = (T_1, T_2, T_3, T_4, T_5) \in \mathbb{R}^5$$

**設計原則：** ハイパーパラメータなし・スカラー集約なし・各成分は独立した次元。

| 成分 | 定義 | 意味 | 証拠分類 |
|------|------|------|---------|
| $T_1$ | $\|\psi_A-\psi_B\|$ | スケール・自由度の乖離 | Derived |
| $T_2$ | $d_B(p_A,p_B)$ | 状態点の幾何学的ずれ | Derived |
| $T_3$ | $D^{(\alpha\neq0)}(p_A\|p_B)$ | 非対称な意味的差 | Derived |
| $T_4$ | $\|\text{Hol}^{(\alpha)}-I\|$ | 多様体の曲率差 | Assumed |
| $T_5$ | $\mathcal{L}(e) = 1-I(X_A;X_B)/H(X_A)$ | 通信路損失率 | Derived |

**⚠️ T₃ の α=0 は禁止**（T₂ との重複が発生するため）

```python
tension = calc_regime_tension_vector(data_A, data_B, alpha_T3=1.0)
```

---

## 8. グランドポテンシャル Φ

$$\Phi = F - \mu N, \quad \mu = 2\sqrt{T_2} \quad [J]$$

$$d\Phi = -S\,dT - P\,dV - N\,d\mu$$

**放置すれば Φ は減少する**（3項全てが減少方向に作用）。

| 変数 | 定義 | 次元 | 証拠分類 |
|------|------|------|---------|
| $S$ | $k_B\ln 2\cdot H_{bit}$ | J/K | Derived |
| $T$ | $\mathrm{tr}(G^{-1})/(d\cdot k_B)$ | K | Derived |
| $P$ | $\mathrm{sign}(d\Phi)\cdot\|d\Phi\|_G$ | J | Derived |
| $V, dV$ | $\int\sqrt{\det G}\,d^d\theta$，$\sqrt{\det G}\,d^d\theta$ | m³相当 | Derived |
| $N$ | 移動する人数 | 人 | Observed |
| $\mu = 2\sqrt{T_2}$ | 1人の移動コスト | J/人 | Derived |

```python
mu  = calc_mu(tension.T2)
Phi = calc_grand_potential(F, mu, N)
```

**証拠分類：** Derived（N が Observed）

---

## 9. 行動条件診断（§4-2・アフォーダンスの代替）

### 設計思想

アフォーダンス（A / P_behavior）は学術的に確立した計算手法を持たないため廃止。
既存の𝓣・Φ・H の読み替えのみで「顧客が行動しやすい状態かどうか」を診断する。
**新しい変数は一切追加しない。**

### 診断の6項目と対応変数

| 診断項目 | 対応変数 | 行動しやすい条件 | 行動しにくい条件 |
|---------|---------|----------------|----------------|
| **到達可能性** | $\mu = 2\sqrt{T_2}$ | $\mu$ 小（移動コスト低） | $\mu$ 大（障壁が大きい） |
| **伝達効率** | $T_5$ | $T_5$ 小（情報損失少） | $T_5$ 大（気づかれない） |
| **方向整合性** | $T_3$ | $T_3$ 小（双方向に自然） | $T_3$ 大（一方向に偏る） |
| **スケール適合** | $T_1$ | $T_1$ 小（文脈が近い） | $T_1$ 大（スケールズレ） |
| **変化の好機** | $d\Phi/dt$ | $d\Phi/dt < 0$（機が熟している） | $d\Phi/dt \approx 0$（変化の動機弱） |
| **文脈の持続性** | $H$ | $H > 0.5$（トレンド持続） | $H < 0.5$（文脈が不安定） |

### 総合判定ロジック

```
全項目○（alert なし、good ≥ 4）→ 「行動が起きやすい状態」
alert ≥ 3 項目              → 「行動が起きにくい状態」
alert 1〜2 項目             → 「一部に障壁あり。優先課題を特定」
その他                       → 「中程度。特定の介入で改善可能」
```

### 各項目の閾値

| 項目 | ○ 基準 | △ 基準 | × 基準 |
|------|--------|--------|--------|
| 到達可能性 | μ < 0.5 | 0.5 ≤ μ < 1.5 | μ ≥ 1.5 |
| 伝達効率 | T₅ < 0.3 | 0.3 ≤ T₅ < 0.6 | T₅ ≥ 0.6 |
| 方向整合性 | \|T₃\| < 0.3 | 0.3 ≤ \|T₃\| < 1.0 | \|T₃\| ≥ 1.0 |
| スケール適合 | T₁ < 0.3 | 0.3 ≤ T₁ < 1.0 | T₁ ≥ 1.0 |
| 変化の好機 | dΦ/dt < -0.01 | -0.01 ≤ dΦ/dt < 0 | dΦ/dt ≥ 0 |
| 文脈の持続性 | H ≥ 0.55 | 0.45 ≤ H < 0.55 | H < 0.45 |

### 計算関数（`sdft_intervention_v005.py`）

```python
diag = calc_action_condition(
    T1=tension.T1, T2=tension.T2, T3=tension.T3, T5=tension.T5,
    H=H_val, phi_current=Phi_now, phi_prev=Phi_prev,
)
# diag.overall_status    → "good" / "warn" / "alert" / "mid"
# diag.overall_message   → 総合判定メッセージ
# diag.table_rows()      → §4-2 テーブル用データ
```

**証拠分類：** Derived（全項目が既存変数の読み替えのため）

---


---

## 10. §4-3 エッジ障壁ヒートマップ（v0.0.5 新規実装）

### 設計思想

アフォーダンス系の変数（normative_barrier / cognitive_barrier 等）を使わず、
T₂・T₅ のエッジ単位集計のみで「各エッジの流れやすさ」を可視化する。

### flow_ease の定義

$$	ext{flow\_ease}_{e} = \left(1 - rac{\mu_e}{\mu_{max}}
ight) \cdot (1 - T_{5,e}) \in [0, 1]$$

| 要素 | 定義 | 証拠分類 |
|------|------|---------|
| $\mu_e = 2\sqrt{T_{2,e}}$ | エッジ $e$ の移動コスト | Derived |
| $T_{5,e}$ | エッジ $e$ の通信路損失率 | Derived |
| $\mu_{max}$ | μ の正規化用最大値 | Assumed |

**積の正当化：** 「動けるが信号が届かない」「信号は届くが動けない」いずれもフローは起きない（AND 条件）。

### 計算関数（`sdft_report_v005.py`）

```python
edge_flows = calc_edge_flows(
    regime_data={"A": data_A, "B": data_B},
    edges=[("A", "B"), ("B", "A")],
    mu_max=3.0,   # Assumed
)
# edge_flows[i].flow_ease, .status, .label
```

### ステータスの判定

| flow_ease | ステータス | 意味 |
|---|---|---|
| ≥ 0.6 | ○ easy | 流れやすい |
| 0.3〜0.6 | △ moderate | 中程度 |
| < 0.3 | × hard | 流れにくい |

---

## 11. Φ₀（介入閾値）の自動設定（v0.0.5 新規実装）

### 設計思想

Φ₀ を外部から手動設定する代わりに、Living System 相の境界条件から自動導出する。

### 自動設定のロジック

```python
phi_0, desc = auto_phi_threshold(state, method="living_system_boundary")
```

**method = "living_system_boundary"（推奨）：**

Living System 相の境界での S（目標値 0.55 付近）を使って：

$$F_{boundary} = U - T \cdot S_{target}$$
$$\Phi_0 = F_{boundary} - \mu \cdot N$$

**その他のオプション：**
- `method="ratio"`：Φ_current × 0.95（簡易設定）
- `method="phase_transition_boundary"`：Noise Dominant 手前の境界

### τ（予測ホライズン）の自動設定

```python
tau, desc = auto_tau(H)
```

| H の範囲 | τ | 理由 |
|---|---|---|
| H ≥ 0.6 | 5 ステップ | 持続傾向が強く長期予測が有効 |
| H ≥ 0.45 | 3 ステップ | 中程度の持続性 |
| H ≥ 0.3 | 2 ステップ | 持続性が低い |
| H < 0.3 | 1 ステップ | 反持続傾向・最短で予測 |

---

## 12. レポート生成パイプライン（v0.0.5 新規実装）

### 全体の流れ

```python
from sdft_report_v005 import ReportContext, generate_report

# StateSnapshot を計算
snap = StateSnapshot(t=0.0, x=data_A, x_B=data_B, N=100.0)
snap.compute()

# 介入最適化ループを実行
states   = generate_state_sequence(x_series, x_B_series, N_series)
phi_0, _ = auto_phi_threshold(snap)
tau, _   = auto_tau(snap.H)
tradeoff = scan_lambda_tradeoff(states, phi_0)
opt_lam, _ = find_optimal_lambda(tradeoff)
result   = run_intervention_loop(states, phi_0, tau=tau, lambda_count=opt_lam)

# ReportContext の生成（全プレースホルダーを自動埋め込み）
ctx = ReportContext.from_state(
    state=snap,
    regime_data={"A": data_A, "B": data_B},
    edges=[("A", "B"), ("B", "A")],
    intervention_result=result,
    opt_lambda=opt_lam,
    H_series=[s.H for s in states],
    title="分析タイトル",
    subject="分析対象",
)

# HTML レポートを生成・保存
html = generate_report(
    ctx=ctx,
    template_path="assets/report_template_v005.html",
    output_path="/mnt/user-data/outputs/report.html",
)
```

### 主要な自動埋め込み対象

```
構造場変数：S/D/H/F/T/U/V/Φ/μ/N/dΦ/dt/Φ₀ → 全て自動
𝓣 ベクトル：T₁〜T₅ + 証拠分類ラベル → 全て自動
位相判定：フェーズバッジ・説明文 → 全て自動
SVG 図：S-D位相図・H時系列 → 全て自動
行動条件診断：6項目テーブル・総合判定 → 全て自動
エッジヒートマップ：flow_ease テーブル + SVG → 全て自動
介入最適化：priority_score・Δu*・推奨介入変数 → 全て自動
証拠分類表：全変数の分類表 → 全て自動
```

## 解析プロセス（ステップバイステップ）

### Step 1：Input Audit
```
・データ種類・粒度・欠損を確認
・各データに証拠分類ラベルを付ける
・レジームA・B の識別
・N（レジーム間を移動する人数）の確認
```

### Step 2：SDFT Mapping
```
ノード = 構成要素（部門・ゾーン・セグメント）
リンク = 接続（強度・方向のみ）
場     = 外部環境（市場・競合・マクロ）
境界   = レジームの境界
レジーム = 同一オペモードを持つ領域
```

### Step 3：構造場変数の計算
```python
H_bit = calc_H_bit(x)
S     = calc_S(x)
D     = calc_D(x)
H_exp = calc_H(x)
G     = calc_fisher_matrix(x)
T     = calc_T_temperature(G)
U     = calc_U_internal_energy(G)
F     = calc_F_free_energy(G, H_bit)
V     = calc_V_volume(G)
```

### Step 4：𝓣 ベクトルの計算
```python
tension = calc_regime_tension_vector(data_A, data_B, alpha_T3=1.0)
# T₁〜T₅ を個別に記録・表示（スカラー集約禁止）
```

### Step 5：グランドポテンシャル Φ の計算
```python
mu  = calc_mu(tension.T2)
Phi = calc_grand_potential(F, mu, N)
```

### Step 6：位相判定
```python
label, css, explanation = classify_phase(S_val, D, H_exp)
```

**相転移予兆（3条件全て揃った場合のみ）：**
$$\text{alert} = \left(\frac{dH}{dt}<0\right) \land \left(\mathrm{Var}(D)\uparrow\right) \land \left(\frac{dS}{dt}>0\right)$$

### Step 7：行動条件診断
```python
diag = calc_action_condition(
    T1=tension.T1, T2=tension.T2,
    T3=tension.T3, T5=tension.T5,
    H=H_exp, phi_current=Phi, phi_prev=Phi_prev,
)
```

**優先解釈順序：**
1. **dH/dt**（最重要先行予兆）
2. **Var(D)**（構造不安定化）
3. **dS/dt**（散逸の進行・遅行指標）
4. **行動条件診断の総合判定**（介入の方向性）
5. **priority_score**（最優先操作変数）

### Step 8：介入最適化
```python
from sdft_intervention_v005 import (
    generate_state_sequence, scan_lambda_tradeoff,
    find_optimal_lambda, run_intervention_loop,
)
states     = generate_state_sequence(x_series, x_B_series, N_series)
tradeoff   = scan_lambda_tradeoff(states, phi_threshold=Phi_0)
opt_lambda, _ = find_optimal_lambda(tradeoff)
result     = run_intervention_loop(states, Phi_0, lambda_count=opt_lambda)
```

---

## 基本動作原則

1. **カノニカルソースを参照する：** `SDFT_Unified_Complete_v0.0.5.html` が最優先
2. **4層証拠分類を全変数に付与する**
3. **削除済み変数（アフォーダンス系含む）を使用しない**
4. **𝓣 をスカラー化しない**（5成分を個別に保持・表示）
5. **F の計算に旧来の $F=U-\Theta S$（Assumed版）を使わない**
6. **行動条件診断は既存変数の読み替えのみで行う**（新変数を追加しない）
7. **エグゼクティブファーストを維持する**（§1 は平易な言語）

---

## 出力仕様

- **形式：** HTML 単一ファイル
- **ファイル名：** `sdft-revised-v005_YYYY-MM-DD-HH-MM.html`
- **テンプレート：** `{スキルディレクトリ}/assets/report_template_v005.html`
- **出力先：** `/mnt/user-data/outputs/`

### レポート構造（7セクション）

```
§1. エグゼクティブサマリー
    技術用語禁止・状況/問題/アクション3件/期待効果/トップリスク

§2. KPI ダッシュボード
    フェーズバッジ（5相）+ 相転移アラート
    KPI カード：S / D / H / F / Φ
    𝓣 の5成分個別カード（T₁〜T₅・スカラー集約なし）

§3. 施策提案・介入優先順位
    Φ に基づく priority_score（∂Φ/∂u_i / c_i）
    λ 最適化結果
    推奨介入変数と必要操作量

§4. 変革実行環境の診断
    4-1. SDFT 写像
    4-2. 行動条件診断（6項目・既存変数の読み替え）← v0.0.5 で新規実装
    4-3. エッジ障壁ヒートマップ（flow_ease = (1-μ/μ_max)×(1-T₅)）
    4-4. 証拠分類表

§5. フェーズ判定と環境変化リスク
    S-D 位相空間図 / H 時系列 / シナリオ表

§6. 総論
    不確実性の高い領域 / 追加で収集すべき観測データ

§7. 巻末用語集 + クレジット（v0.0.5）
```

### §4-2 行動条件診断の表示仕様

```
【行動条件診断】
総合判定：[○ 行動が起きやすい / △ 一部に障壁あり / × 行動が起きにくい]
  [メッセージ]

┌──────────────┬──────────────────────┬──────────────────────────────────┐
│ 診断項目     │ 対応変数             │ 解釈                             │
├──────────────┼──────────────────────┼──────────────────────────────────┤
│ ○ 到達可能性  │ μ=2√T₂=[値]          │ 移動コスト低い。顧客は動きやすい。 │
│ × 伝達効率   │ T₅=[値]              │ 通信路損失大。気づかれにくい。     │
│ △ 方向整合性  │ T₃=[値]              │ やや非対称。一方向への傾き。        │
│ ○ スケール適合│ T₁=[値]              │ スケール差小。文脈が近い。          │
│ ○ 変化の好機  │ dΦ/dt=[値]           │ Φ 減少中。変化の機が熟している。   │
│ △ 文脈の持続性│ H=[値]               │ H≈0.5。文脈がやや不安定。          │
└──────────────┴──────────────────────┴──────────────────────────────────┘
```

---


## v0.0.5 で修正したバグ（v0.0.4 比）

| バグ | 内容 | 修正 |
|------|------|------|
| `calc_H()` の `/2.0` | std を使う場合 slope = H のため /2.0 は不要 → H が半分になっていた | `/2.0` を削除 |
| `calc_fisher_matrix()` 非標準化 | 元データで G を計算すると σ²≫1 により T≈10²⁹K になっていた | データを標準化してから G を計算 |
| `calc_T4()` 常に 0 | 標準化後 G_A=G_B となり差分が消えていた | 非標準化の分散比（対数行列式差）で計算するよう変更 |

## ペンディング項目（中優先・v0.1 以降）

| 項目 | 現状 | 理由 |
|------|------|------|
| Fisher 行列の一般的な経験的推定 | 正規分布仮定の閉形式（2×2） | 精度向上・任意分布への対応 |
| T₅ の KSG 推定量 | ガウス近似の相互情報量 | 非線形依存関係の捕捉 |
| H の精度向上 | ラグ別std法（短い系列で不安定） | DFA など高精度手法への切替 |
| 介入コスト c_i の設定ガイドライン | 全て Assumed | ドメイン固有の設定根拠 |
| V の数値積分 | √det G × θ_range^d（近似） | 厳密な積分 |
| T₄ の厳密なホロノミー計算 | 分散比の近似 | Fisher-Rao 接続による経路積分 |
| P の自己整合解法 | 単純差分による近似 | 最適化ループ |

## 品質ゲート（出力前チェックリスト）

```
□ §1 が平易な言語で書かれているか（技術用語なし）
□ §2 の 𝓣 が5成分個別で表示されているか（スカラー集約なし）
□ §2 の 𝓣 各成分に証拠分類が付いているか
□ §2 に Φ の値と dΦ/dt が表示されているか
□ T₃ に使用した α 値が明示されているか（α=0 禁止）
□ §3 に priority_score が ∂Φ/∂u_i に基づいて表示されているか
□ §4-2 が行動条件診断テーブルで表示されているか（未対応表示は禁止）
□ §4-2 の各項目に対応変数と解釈が記載されているか
□ §4-3 エッジ障壁ヒートマップに flow_ease テーブルと SVG が含まれているか
□ §4-4 に全変数の証拠分類が含まれているか
□ §5 に S-D 位相図と H 時系列が含まれているか
□ A / P_behavior が出力に含まれていないか（廃止済み）
□ アフォーダンス関連変数（salience 等）が出力に含まれていないか
□ F の計算に旧来の F=U-ΘS（Assumed版）を使っていないか
□ 行動条件診断に新しい変数が導入されていないか（既存変数の読み替えのみ）
□ 巻末用語集が含まれているか（アフォーダンス関連用語は削除）
□ クレジットが v0.0.5 になっているか
```

### ハード制約（絶対禁止）

- A / P_behavior を計算・表示する
- アフォーダンス関連変数（salience/feasibility/legitimacy/reward/mutuality/cognitive_cost）を使用する
- 削除済み変数（q/W/ε/μ_maxwell/H_eff 等）を使用する
- 𝓣 にハイパーパラメータ（κ/λ/η）を掛ける
- 𝓣 をスカラーに集約する
- T₃ の α に 0 を使う
- F = U − ΘS（旧来の Assumed 版）を使用する
- 行動条件診断に新しい変数を導入する
- HTML をレスポンスに直接出力する

---

## 巻末用語集（v0.0.5）

| 用語 | 説明 |
|------|------|
| **$H_{bit}$** | $-\sum p_i\log_2 p_i$ [bit]。シャノンエントロピー。 |
| **$S$** | $k_B\ln 2\cdot H_{bit}$ [J/K]。熱力学的エントロピー。 |
| **$D$** | ヒグチ法。1.0〜2.0。時系列の複雑性。 |
| **$H$** | Hurst指数。0〜1。相転移の先行指標。 |
| **$G$** | Fisher情報行列。統計多様体上の計量テンソル。 |
| **$T$** | $\mathrm{tr}(G^{-1})/(d\cdot k_B)$ [K]。情報的温度。 |
| **$U$** | $\mathrm{tr}(G)/d$ [J相当]。内部エネルギー。 |
| **$F$** | $U-T\cdot S$ [J]。自由エネルギー。 |
| **$\Phi$** | $F-\mu N$ [J]。グランドポテンシャル。放置すれば減少する。 |
| **$\mu$** | $2\sqrt{T_2}$ [J/人]。化学ポテンシャル。1人の移動コスト。 |
| **$\vec{\mathcal{T}}$** | $(T_1,...,T_5)\in\mathbb{R}^5$。レジーム間テンション5次元ベクトル。 |
| **$T_1$** | 対数分配関数差。スケール・自由度の乖離。 |
| **$T_2$** | Bhattacharyya距離。状態点の幾何学的ずれ。 |
| **$T_3$** | α-ダイバージェンス（α≠0）。非対称な意味的差。 |
| **$T_4$** | 統計的ホロノミー。多様体の曲率差。 |
| **$T_5$** | 通信路損失率 $1-I(X_A;X_B)/H(X_A)$。情報損失割合。 |
| **行動条件診断** | 𝓣・Φ・H の読み替えで行動しやすさを診断。新変数なし。 |
| **到達可能性** | $\mu=2\sqrt{T_2}$ の低さ。移動コストが低いほど行動しやすい。 |
| **伝達効率** | $1-T_5$。情報損失が少ないほど気づきやすい。 |
| **方向整合性** | $|T_3|$ の小ささ。双方向に自然な移動。 |
| **スケール適合** | $T_1$ の小ささ。2レジームの文脈が近い。 |
| **変化の好機** | $d\Phi/dt < 0$。Φ が減少中＝変化の機が熟している。 |
| **文脈の持続性** | $H > 0.5$。行動が起きても継続しやすい状態。 |
| **Living System** | 中S・中D・中〜高H。健全な動的平衡。 |
| **Frozen Order** | 低S・低D・高H。安定だが変化なし。 |
| **Runaway Growth** | 中S・高D・高H。過成長の兆候。 |
| **Noise Dominant** | 高S・高D・低H。構造崩壊の前兆。 |
| **Collapse** | 高S・不安定D・低H。持続不能。 |
| **Observed** | 実測値 |
| **Derived** | 観測値から確定的に算出した値 |
| **Assumed** | 分析者の仮定値 |
| **Hypothesis** | 解釈的・推論的推定 |

クレジット：株式会社アドインテ SDFT位相空間関係性モデル sdft-revised v0.0.5
