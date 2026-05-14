# SDFT v0.2 出力契約

## 1. 出力ファイル

| 項目 | 仕様 |
|---|---|
| 形式 | HTML 単一ファイル |
| ファイル名 | `sdft-revised-v02_YYYY-MM-DD-HH-MM.html` |
| 出力先 | `/mnt/user-data/outputs/` |
| 文字コード | UTF-8 |
| テンプレート | `assets/report_template_v02.html` |

## 2. 必須セクション

| § | タイトル | 必須プレースホルダー |
|---|---|---|
| 1 | エグゼクティブサマリー | EXEC_SITUATION, EXEC_PROBLEM, EXEC_ACTION_1-3, EXEC_IMPACT, EXEC_TOP_RISK |
| 2 | KPI ダッシュボード | SG_VALUE, D_VALUE, LAMBDA1_VALUE, FGEOM_VALUE, MU_VALUE, PHI_VALUE, T1-T5 |
| 3 | 施策提案・介入優先順位 | SENSITIVITY_*, PRIORITY_*, BEST_VARIABLE, PHI_ENERGY_GRAPH_SVG |
| 4 | 変革実行環境の診断 | ACTION_CONDITION_ROWS, EDGE_HEATMAP_*, EVIDENCE_ROWS |
| 5 | フェーズ判定と環境変化リスク | SD_X/Y, LAM_POLYLINE, SCENARIO_ROWS |
| 6 | 総論 | UNCERTAINTY_DESC, ADDITIONAL_OBS_ROWS |
| 7 | 巻末用語集 | （テンプレート固定） |

## 3. 品質ゲート

レポート出力前に以下を全て確認すること：

```
□ §2 の 𝓣 が5成分個別で表示されているか（スカラー集約なし）
□ §2 に Φ の値と dΦ/dt が表示されているか
□ §3 に priority_score テーブルがあるか
□ §4-2 に行動条件診断テーブル（6項目）があるか
□ λ₁（リアプノフ指数）が使用されているか（H / Hurst 指数は使用禁止）
□ Fisher 行列 G が使用されていないか（v0.3 まで保留）
□ 温度 T [K] が使用されていないか（v0.2 では廃止）
□ バージョン表記が v0.2 になっているか
□ テストが全件パスしているか（12件）
```

## 4. 変数の単位

v0.2 では熱力学的単位系（J, K, J/K）を**使用しない**。
全変数は無次元量または幾何学的量として扱う。

| 変数 | 単位 |
|---|---|
| S_g | 無次元（対数スケール） |
| D | 無次元 |
| λ₁ | 無次元（1/step） |
| λ₂ = F_geom | 無次元 |
| κ | 無次元 |
| T₁-T₅ | 無次元 |
| μ = T₂ | 距離（埋め込み空間の単位） |
| N | 人 |
| Φ | 無次元 |

## 5. 証拠分類ラベル

全変数に以下のいずれかを付与すること：

| ラベル | 定義 |
|---|---|
| Observed | 実データから直接計測 |
| Derived | 観測値から確定的に算出 |
| Assumed | 分析者が設定した仮定 |
| Hypothesis | 解釈的・推論的推定 |
