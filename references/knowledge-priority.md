# SDFT revised v0.0.6 — Knowledge Priority

## 1. 知識参照の優先順位

```
優先度1：references/SDFT_Unified_Complete_v0.0.6.html（カノニカルソース）
優先度2：本 SKILL.md
優先度3：ユーザーから直接提供されたデータファイル（実測値）
優先度4：scripts/ 配下の Python 実装（参照実装）
```

複数データが競合する場合：

```
観測値 > 集計導出値 > 古いデータ
カノニカルソース > SKILL.md
SKILL.md > Python 実装
```

不確実性を維持し、強制的な断定をしない。

## 2. 4層証拠分類

| ラベル | 定義 | 例 |
|--------|------|-----|
| **Observed** | 実データから直接計測した値 | N（移動する人数）、観測時系列 |
| **Derived** | 観測値から確定的に算出した値 | S, D, H, T, U, F, V, Φ, μ, T₁, T₂, T₃, T₅ |
| **Assumed** | 分析者が設定した仮定値 | T₄（近似推定）、介入コスト c_i、μ_max |
| **Hypothesis** | 解釈的・推論的推定 | データがない場合の代理推定 |

## 3. 変数の証拠分類デフォルト

| 変数 | デフォルト分類 | 計算根拠 |
|------|--------------|---------|
| H_bit | Derived | シャノンエントロピー（時系列から計算） |
| S | Derived | k_B·ln2·H_bit |
| D | Derived | ヒグチ法（時系列から計算） |
| H | Derived | ラグ別 std 法（時系列から計算） |
| G | Derived | Fisher 行列（標準化データの閉形式） |
| T | Derived | tr(G⁻¹)/(d·k_B) |
| U | Derived | tr(G)/d |
| F | Derived | U - T·S |
| V | Derived | √det(G)·θ_range^d（θ_range は Assumed） |
| N | Observed | 直接計測 |
| μ | Derived | 2√T₂ |
| Φ | Derived | F - μN |
| T₁ | Derived | 経験的 ψ 差 |
| T₂ | Derived | Bhattacharyya 距離 |
| T₃ | Derived | α-ダイバージェンス（α ≠ 0） |
| T₄ | Assumed | 非標準化分散の対数行列式差による近似 |
| T₅ | Derived（同時観測時）/ Assumed（推定時） | 相互情報量推定 |
| 介入コスト c_i | Assumed | ドメイン固有・外部設定 |
| μ_max | Assumed | 正規化用・既定 3.0 |

## 4. データ欠損時の対応

時系列データが存在しない場合、代理変数を **Hypothesis** ラベルで推定する：

```
S → 行動多様性・需要変動係数から推定
D → 導線複雑度・競合多層性から推定
H → トレンド継続期間から推定
```

代理変数を使う場合は必ずその旨をレポートの証拠分類表に明記すること。

## 5. 証拠の劣化を防ぐ原則

```
1. Derived 変数を Assumed として扱わない
2. Assumed 値を Derived として表示しない
3. Hypothesis（推論）と Observed（実測）を混同しない
4. 代理推定を行った場合は必ず Hypothesis ラベルを付与する
5. 介入コスト c_i のような Assumed 値は出力に明記する
```
