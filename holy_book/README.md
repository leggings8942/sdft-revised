# holy_book — SDFT v1.4.2 原典アーカイブ

## 位置づけ

`holy_book` は SDFT の原典（v1.4.2）を保存したアーカイブです。
**sdft-revised（v0.0.x）のカノニカルソースではありません。**

sdft-revised のカノニカルソースは：
```
{スキルディレクトリ}/references/SDFT_Unified_Complete_v0.0.2.html
```

---

## holy_book を参照すべき状況

以下の状況では holy_book を動的に参照してください：

### 参照すべき状況 ✅

| 状況 | 参照先ファイル |
|------|--------------|
| 𝓣 の旧定義と新定義を比較・説明したい | `SDFT_Unified_Complete_v1.4.2.html` |
| v1.4.2 の変数定義を正確に引用したい | `references/sdft_theory_v1.4.md` |
| v1.4.2 の計算実装を参照したい | `scripts/sdft_unified_v1_4_A_reference.py` |
| 旧バージョンとの互換性を確認したい | `references/output-contract-v1.4.2.md` |
| 削除済み変数の元の定義を確認したい | `SDFT_Unified_Complete_v1.4.2.html` |
| データマッピングの参考にしたい | `references/data-mapping.md` |
| 品質チェックの基準を参考にしたい | `references/quality-checklist-v1.4.2.md` |

### 参照不要な状況 ❌

| 状況 | 理由 |
|------|------|
| 通常の分析・レポート生成時 | sdft-revised の定義で完結する |
| S/D/H/F の計算時 | v0.0.2 と v1.4.2 で定義が同一 |
| 位相判定時 | v0.0.2 と v1.4.2 で定義が同一 |
| 削除済み変数の使用 | 削除済み・参照しても使用禁止 |

---

## ⚠️ 重要：参照時の制約

```
holy_book の内容は「参照・比較・説明」のためのものです。
holy_book の式・変数・実装を sdft-revised の分析に直接使用してはなりません。

絶対禁止：
  - 削除済み変数（q/W/ε/μ/H_eff 等）の使用
  - 𝓣 の v1.4.2 定義式を分析に使用すること
  - アフォーダンス層（A/P_behavior）の使用
```

---

## ファイル構成

```
holy_book/
├── README.md                                    ← このファイル
├── SDFT_Unified_Complete_v1.4.2.html            ← 原典理論書（最重要）
├── references/
│   ├── sdft_theory_v1.4.md                      ← 理論サマリー
│   ├── knowledge-priority-v1.4.2.md             ← 知識優先順位
│   ├── output-contract-v1.4.2.md                ← 出力契約
│   ├── data-mapping.md                          ← データマッピング
│   ├── quality-checklist-v1.4.2.md              ← 品質チェックリスト
│   ├── executive-summary-template-v1.4.3.md     ← エグゼクティブサマリーテンプレート
│   ├── business-impact-template-v1.4.3.md       ← ビジネスインパクトテンプレート
│   └── parameter-provenance-template-v1.4.2.md  ← パラメータ来歴テンプレート
└── scripts/
    ├── sdft_unified_v1_4_A_reference.py          ← v1.4.2 参照実装（最重要）
    ├── sdft_compute_v1.4.py                      ← 計算ヘルパー
    └── validate_report.py                        ← バリデーション
```
