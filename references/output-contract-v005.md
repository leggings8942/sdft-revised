# output-contract v0.0.5

## 必須セクション（7セクション）

| セクション | 必須要素 | 状態 |
|-----------|---------|------|
| §1 エグゼクティブサマリー | 状況/問題/アクション3件/効果/リスク | 必須 |
| §2 KPI ダッシュボード | S/D/H/F/Φ + 𝓣5成分個別 + フェーズバッジ | 必須 |
| §3 施策提案 | priority_score（∂Φ/∂u_i）+ 実装待ち明記 | 必須 |
| §4-1 SDFT写像 | ノード/リンク/場/境界/レジーム | 必須 |
| §4-2 アフォーダンス層 | v0.0.5 未対応表示 | 必須 |
| §4-3 障壁ヒートマップ | v0.0.5 未対応表示 | 必須 |
| §4-4 証拠分類表 | 全変数（S/D/H/T/U/F/Φ/μ/N/T₁〜T₅） | 必須 |
| §5 フェーズ判定 | S-D位相図/H時系列/シナリオ表 | 必須 |
| §6 総論 | 不確実性/追加観測 | 必須 |
| §7 用語集 | 全用語 + クレジット v0.0.5 | 必須 |

## テンプレート

`assets/report_template_v004.html`

## 新規プレースホルダー（v0.0.5 追加分）

### Φ 関連
- `{{PHI_VALUE}}`：Φ の値 [J]
- `{{PHI_DELTA}}`：dΦ/dt の推定値
- `{{PHI_THRESHOLD}}`：介入閾値 Φ₀
- `{{PHI_PREDICTED}}`：予測 Φ（τ後）
- `{{N_VALUE}}`：人数 N [人]
- `{{MU_VALUE}}`：化学ポテンシャル μ = 2√T₂

### §3 介入設計
- `{{INTERVENTION_NEEDED}}`：「介入必要」または「介入不要」
- `{{DELTA_PHI_NEED}}`：ΔΦ_need
- `{{SENSITIVITY_N}}`：∂Φ/∂N
- `{{SENSITIVITY_T2}}`：∂Φ/∂T₂
- `{{SENSITIVITY_THETA}}`：∂Φ/∂θ
- `{{PRIORITY_N}}`：priority_score（N）
- `{{PRIORITY_T2}}`：priority_score（T₂）
- `{{PRIORITY_THETA}}`：priority_score（θ）
- `{{PRIORITY_N_PCT}}`：バー幅（0-100）
- `{{PRIORITY_T2_PCT}}`：バー幅（0-100）
- `{{PRIORITY_THETA_PCT}}`：バー幅（0-100）
- `{{BEST_VARIABLE}}`：最優先操作変数
- `{{DELTA_U_OPTIMAL}}`：必要操作量
- `{{BEST_N_BADGE}}`：最優先なら badge HTML
- `{{BEST_T2_BADGE}}`：最優先なら badge HTML
- `{{BEST_THETA_BADGE}}`：最優先なら badge HTML
- `{{COST_N}}`：介入コスト c_N（Assumed）
- `{{COST_T2}}`：介入コスト c_T₂（Assumed）
- `{{COST_THETA}}`：介入コスト c_θ（Assumed）

### 𝓣 成分
- `{{T1_VALUE}}`, `{{T1_EVIDENCE}}`
- `{{T2_VALUE}}`, `{{T2_EVIDENCE}}`
- `{{T3_VALUE}}`, `{{T3_EVIDENCE}}`, `{{T3_ALPHA}}`
- `{{T4_VALUE}}`, `{{T4_EVIDENCE}}`
- `{{T5_VALUE}}`, `{{T5_EVIDENCE}}`

## KPI STATUS_CSS 判定基準

| 変数 | alert | warn | good | mid |
|------|-------|------|------|-----|
| S | >0.8 | 0.6-0.8 | <0.4 | 0.4-0.6 |
| D | >1.8 | 1.6-1.8 | 1.0-1.4 | 1.4-1.6 |
| H | <0.3 | 0.3-0.45 | >0.6 | 0.45-0.6 |
| F | 文脈依存 | — | 文脈依存 | — |

## v0.0.5 の変更点（v0.0.4 からの差分）

### 廃止
- A（アフォーダンス）/ P_behavior：全箇所で使用禁止
- calc_affordance()：削除済み

### §4-2 の変更
v0.0.4：未対応表示
v0.0.5：行動条件診断テーブル（既存変数の読み替え）

### 新規プレースホルダー（§4-2）
- `{{ACTION_CONDITION_ROWS}}`：6行のテーブル行 HTML
- `{{ACTION_OVERALL_CSS}}`：総合判定の CSS クラス
- `{{ACTION_OVERALL_STATUS}}`：総合判定テキスト
- `{{ACTION_OVERALL_MESSAGE}}`：総合判定メッセージ
