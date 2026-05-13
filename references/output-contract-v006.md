# SDFT revised v0.0.6 — Output Contract

## 1. レポート出力仕様

- **形式：** HTML 単一ファイル
- **ファイル名：** `sdft-revised-v006_YYYY-MM-DD-HH-MM.html`
- **テンプレート：** `assets/report_template_v006.html`
- **出力先：** `/mnt/user-data/outputs/`
- **生成スクリプト：** `scripts/sdft_report_v006.py`

## 2. 必須セクション

| セクション | 必須要素 |
|-----------|---------|
| §1 エグゼクティブサマリー | 状況・問題・アクション3件・期待効果・トップリスク |
| §2 KPI ダッシュボード | S/D/H/F + Φ + μ + 𝓣5成分個別 + フェーズバッジ + 相転移アラート |
| §3 施策提案・介入優先順位 | priority_score テーブル + Φエネルギー-時刻グラフ + 介入詳細テーブル |
| §4-1 SDFT 写像 | ノード/リンク/場/境界/レジーム |
| §4-2 行動条件診断 | 6項目テーブル + 総合判定 |
| §4-3 エッジ障壁ヒートマップ | flow_ease テーブル + SVG バーグラフ |
| §4-4 証拠分類表 | 全変数（H_bit/S/D/H/T/U/F/V/Φ/μ/N/T₁〜T₅） |
| §5 フェーズ判定 | S-D 位相図 + H 時系列 + 予測シナリオ表 |
| §6 総論 | 不確実性 + 追加観測の推奨 |
| §7 巻末用語集 + クレジット | 全用語 + クレジット v0.0.6 |

## 3. プレースホルダー一覧（87項目）

### ヘッダー（5項目）
- `{{REPORT_TITLE}}` / `{{ANALYSIS_DATE}}` / `{{SUBJECT}}` / `{{ANALYSIS_PURPOSE}}` / `{{PHASE_LABEL}}`

### §1 エグゼクティブサマリー（7項目）
- `{{EXEC_SITUATION}}` / `{{EXEC_PROBLEM}}` / `{{EXEC_IMPACT}}` / `{{EXEC_TOP_RISK}}`
- `{{EXEC_ACTION_1}}` / `{{EXEC_ACTION_2}}` / `{{EXEC_ACTION_3}}`

### §2 KPI ダッシュボード（25項目）
- `{{PHASE_CSS}}` / `{{PHASE_EXPLANATION}}` / `{{PHASE_TRANSITION_ALERT}}`
- KPI カード：`{{S_VALUE/S_STATUS_CSS/S_STATUS}}`、D / H / F 同様
- Φ 関連：`{{MU_VALUE}}` / `{{PHI_VALUE}}` / `{{PHI_DELTA}}` / `{{PHI_THRESHOLD}}` / `{{N_VALUE}}` / `{{PHI_PREDICTED}}`
- 𝓣 ベクトル：`{{T1_VALUE}}` 〜 `{{T5_VALUE}}` / `{{T1_EVIDENCE}}` 〜 `{{T5_EVIDENCE}}` / `{{T3_ALPHA}}`

### §3 介入優先順位（18項目）
- 判定：`{{INTERVENTION_NEEDED}}` / `{{DELTA_PHI_NEED}}` / `{{BEST_VARIABLE}}` / `{{DELTA_U_OPTIMAL}}`
- 感度：`{{SENSITIVITY_N}}` / `{{SENSITIVITY_T2}}` / `{{SENSITIVITY_THETA}}`
- priority：`{{PRIORITY_N}}` / `{{PRIORITY_T2}}` / `{{PRIORITY_THETA}}`
- バー幅：`{{PRIORITY_N_PCT}}` / `{{PRIORITY_T2_PCT}}` / `{{PRIORITY_THETA_PCT}}`
- バッジ：`{{BEST_N_BADGE}}` / `{{BEST_T2_BADGE}}` / `{{BEST_THETA_BADGE}}`
- コスト：`{{COST_N}}` / `{{COST_T2}}` / `{{COST_THETA}}`
- グラフ：`{{PHI_ENERGY_GRAPH_SVG}}` / `{{PHI_INTERVENTION_TABLE}}`

### §4 診断（10項目）
- §4-1：`{{MAPPING_NODES/LINKS/FIELD/BOUNDARY/REGIME}}`
- §4-2：`{{ACTION_CONDITION_ROWS}}` / `{{ACTION_OVERALL_CSS}}` / `{{ACTION_OVERALL_STATUS}}` / `{{ACTION_OVERALL_MESSAGE}}`
- §4-3：`{{EDGE_HEATMAP_TABLE}}` / `{{EDGE_HEATMAP_SVG}}`
- §4-4：`{{EVIDENCE_ROWS}}`

### §5 フェーズ判定（8項目）
- SVG 位相図：`{{SD_X}}` / `{{SD_Y}}` / `{{SD_Y_LABEL}}` / `{{H_POLYLINE}}` / `{{H_CURRENT_X}}` / `{{H_CURRENT_Y}}` / `{{H_CURRENT_Y_LABEL}}`
- `{{SCENARIO_ROWS}}`

### §6 総論（2項目）
- `{{UNCERTAINTY_DESC}}` / `{{ADDITIONAL_OBS_ROWS}}`

## 4. KPI ステータス判定基準

### 構造場変数
| 変数 | alert | warn | mid | good |
|------|-------|------|-----|------|
| S（正規化） | ≥0.8 | 0.6-0.8 | 0.4-0.6 | <0.4 |
| D | ≥1.8 | 1.6-1.8 | 1.4-1.6 | 1.0-1.4 |
| H | <0.3 | 0.3-0.45 | 0.45-0.6 | ≥0.6 |
| F | < 0（要警戒） | — | — | > 0（健全） |

### 𝓣 成分
| 成分 | alert | warn | good |
|------|-------|------|------|
| T₁ | ≥0.5 | ≥0.2 | <0.2 |
| T₂ | ≥1.0 | ≥0.3 | <0.3 |
| T₃ | ≥0.5 | ≥0.1 | <0.1 |
| T₄ | ≥0.5 | ≥0.2 | <0.2 |
| T₅ | ≥0.7 | ≥0.3 | <0.3 |

## 5. 出力品質ゲート

```
□ 全プレースホルダーが置換されている（残留ゼロ）
□ §1 が平易な言語で書かれている（技術用語禁止）
□ §2 の 𝓣 5成分が個別に表示されている（スカラー集約禁止）
□ §2 に Φ の値と dΦ/dt が表示されている
□ §3 に priority_score テーブルがある
□ §3 に Φエネルギー-時刻グラフが含まれている
□ §3 に介入詳細テーブルが含まれている
□ §4-1 に SDFT 写像（5項目）がある
□ §4-2 に行動条件診断（6項目）テーブルがある
□ §4-3 にエッジ障壁ヒートマップ（テーブル + SVG）がある
□ §4-4 に全変数の証拠分類が含まれている
□ §5 に S-D 位相図と H 時系列がある
□ T₃ の α 値が明示されている（α=0 禁止）
□ アフォーダンス系の変数が出力に含まれていない
□ 物理アナロジー系の変数が出力に含まれていない
□ クレジットが v0.0.6 になっている
```
