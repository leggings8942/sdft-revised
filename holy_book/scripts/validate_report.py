#!/usr/bin/env python3
"""SDFT v1.4.2 レポート品質ゲート検証スクリプト"""
import sys, pathlib, re

REQUIRED_SECTIONS = [
    "エグゼクティブ・サマリー",
    "施策提案",
    "変革実行環境の診断",
    "フェーズ判定",
    "情報伝達効率",
    "用語集",
    "証拠分類",
    "アフォーダンス",
    "P_behavior",
]

REQUIRED_KPI = ["S:", "D:", "H:", "F:", "𝓣", "A:", "P_behavior"]

REQUIRED_CSS_CLASSES = [
    "kpi-grid-7",
    "affordance-grid",
    "heatmap-table",
    "vector-diagram",
    "phase-badge",
    "glossary-grid",
]

FORBIDDEN_PLACEHOLDERS = re.compile(r'\{\{[A-Z_]+\}\}')


def main():
    if len(sys.argv) < 2:
        print("usage: validate_report.py <report.html>")
        raise SystemExit(2)

    p = pathlib.Path(sys.argv[1])
    if not p.exists():
        print(f"ERROR: File not found: {p}")
        raise SystemExit(1)

    text = p.read_text(encoding="utf-8")
    errors = []
    warnings = []

    # Check required sections
    for sec in REQUIRED_SECTIONS:
        if sec not in text:
            errors.append(f"必須セクション未発見: '{sec}'")

    # Check KPI indicators
    for kpi in REQUIRED_KPI:
        if kpi not in text:
            errors.append(f"KPI指標未発見: '{kpi}'")

    # Check CSS classes (template structure)
    for cls in REQUIRED_CSS_CLASSES:
        if cls not in text:
            warnings.append(f"テンプレートCSSクラス未使用: '{cls}' (レイアウト崩れの可能性)")

    # Check for unfilled placeholders
    remaining = FORBIDDEN_PLACEHOLDERS.findall(text)
    if remaining:
        unique = list(dict.fromkeys(remaining))[:10]
        errors.append(f"未置換プレースホルダー発見: {unique}")

    # Check print CSS
    if "@media print" not in text:
        warnings.append("@media print CSS が未埋め込み（印刷非対応）")

    # Check credit
    if "株式会社アドインテ" not in text or "SDFTv1.4.2" not in text:
        warnings.append("クレジット（株式会社アドインテ SDFTv1.4.2）が未記載")

    # Check v1.4.2 badge
    if "v1.4.2" not in text:
        warnings.append("v1.4.2 バッジが未表示")

    # Report
    if errors:
        print("=== ERRORS ===")
        for e in errors:
            print(f"  [ERROR] {e}")
    if warnings:
        print("=== WARNINGS ===")
        for w in warnings:
            print(f"  [WARN]  {w}")

    if not errors and not warnings:
        print("OK — 品質ゲートすべて通過")
    elif not errors:
        print(f"\nOK (警告 {len(warnings)}件) — エラーなし")
    else:
        print(f"\nFAIL — エラー {len(errors)}件 / 警告 {len(warnings)}件")
        raise SystemExit(1)


if __name__ == "__main__":
    main()
