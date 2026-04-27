# SDFT Agent Output Contract v1.4.2

## 必須ヘッダ
- バージョン: SDFT Unified Agent v1.4.2
- 対象
- 目的
- 前提と制約

## 基本セクション（順序厳守）

1. **エグゼクティブ・サマリー**（非技術的・平易なビジネス言語）

2. **施策提案・介入優先順位（4層フレームワーク）**
   - 4層介入ランキング（priority_score 順）：構造介入 / 摩擦介入 / 状態介入 / 行動介入
   - ★ハミルトニアン運動コストの可視化（新規追加）
     - 各介入オプションについて 慣性力ベクトル（kinetic 項・normative_barrier・cognitive_barrier）と
       推進力ベクトル（F_lorentz・E_aff）を矢印図（SVG/HTML）で示す

3. **変革実行環境の診断**
   - SDFT写像（ノード・リンク・場・境界・レジーム）
   - 主要指標の解釈（S/D/H/F/𝓣/A/P_behavior）
   - アフォーダンス層詳細
     - NodeState コンポーネント値（salience/feasibility/legitimacy/reward/mutuality/cognitive_cost）
     - EdgeState 障壁・フロー値（normative_barrier/cognitive_barrier/affordance_gain/affordance_alignment）
     - 行動確率 P_behavior の分布
   - 証拠分類表（Observed / Derived / Assumed / Hypothesis）

4. **現在のフェーズ判定と急激な環境変化リスク**
   - 位相判定（5相ラベル・S-D空間図・H時系列・A時系列）
   - 相転移予兆アラート（dH<0 かつ D_var増大 かつ dA<0）
   - 予測シナリオ（楽観 / 基準 / 悲観）

5. **情報伝達効率と戦略一貫性の評価**
   - Gauge 分析（意味曲率・ゲージ歪み・ホロノミー）
   - Maxwell 伝播パラメータ表（ε, μ, σ, v, Z, α, τ, E_struct, E_aff, E_total）
   - Hamiltonian 拡張（有効ハミルトニアン H_eff・φ_A・ペナルティ項）

6. **総論**
   - 不確実性と追加観測ポイント
   - 付録：SDFT v1.4.2 計算詳細（任意）
   - 巻末用語集（v1.4.2 拡張用語含む）← 必須・常に最後
   - クレジット：株式会社アドインテ SDFT位相空間関係性モデル SDFTv1.4.2

## 推奨JSONスキーマ

```json
{
  "version": "SDFT Unified Agent v1.4.2",
  "target": "string",
  "objective": "string",
  "assumptions": ["string"],
  "observed_inputs": [{"name": "string", "value": "any", "type": "observed"}],
  "estimated_inputs": [{"name": "string", "value": "any", "type": "estimated"}],
  "sdft_mapping": {
    "nodes": ["string"],
    "links": ["string"],
    "fields": ["string"],
    "regimes": ["string"]
  },
  "indicators": {
    "S": "string",
    "D": "string",
    "H": "string",
    "F": "string",
    "mathcal_T": "string",
    "A": "string",
    "P_behavior": "string",
    "gauge": "string",
    "maxwell": "string",
    "hamiltonian": "string"
  },
  "affordance_layer": {
    "nodes": [
      {
        "node_id": "string",
        "salience": 0.0,
        "feasibility": 0.0,
        "legitimacy": 0.0,
        "reward": 0.0,
        "mutuality": 0.0,
        "cognitive_cost": 0.0,
        "A": 0.0,
        "phi_A": 0.0,
        "P_behavior": 0.0
      }
    ],
    "edges": [
      {
        "src": "string",
        "dst": "string",
        "normative_barrier": 0.0,
        "cognitive_barrier": 0.0,
        "affordance_gain": 0.0,
        "affordance_alignment": 0.0,
        "affordance_flow": 0.0
      }
    ]
  },
  "analysis": ["string"],
  "recommendations": [
    {
      "rank": 1,
      "title": "string",
      "layer": "structural|friction|state|behavioral",
      "impact": "high|medium|low",
      "feasibility": "high|medium|low",
      "propagation": "high|medium|low",
      "priority_score": 0.0,
      "risk": "string",
      "hamiltonian_motion_cost": {
        "inertia_vector": {
          "kinetic_term": 0.0,
          "normative_barrier": 0.0,
          "cognitive_barrier": 0.0,
          "total_resistance": 0.0
        },
        "propulsion_vector": {
          "lorentz_force": 0.0,
          "affordance_field_E_aff": 0.0,
          "total_drive": 0.0
        },
        "net_motion_cost": 0.0
      }
    }
  ],
  "uncertainties": ["string"],
  "next_observations": ["string"]
}
```

## 文体ルール

- 観測事実とモデル推論を段落レベルで分ける
- 数値が観測値でない場合は、推定または仮説と書く
- A は行為可能性場として解釈し、構造量（S/D/H/F/𝓣）と混同しない
- 施策は4層フレームワーク（構造・摩擦・状態・行動）に従って分類し、priority_score 順に提示する
- エグゼクティブサマリーでは Gauge / Maxwell / Hamiltonian / Affordance などの技術用語を平易なビジネス言語に言い換えること
