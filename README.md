# ペロブスカイト超電導 × トリテトラ理論（TTT）

[![Theory](https://img.shields.io/badge/Theory-TriTetra%20TTT-7c6af7?style=for-the-badge)](https://github.com/kiki054-n/TTT)
[![License](https://img.shields.io/badge/License-MIT-2ecc71?style=for-the-badge)](LICENSE)
[![Status](https://img.shields.io/badge/Status-研究中-e74c3c?style=for-the-badge)]()

> **ベクトルの総和はゼロである — この一原理でペロブスカイト超電導の構造・電子・熱的安定性を俯瞰する**

本リポジトリは、[TTT（トリテトラ理論）](https://github.com/kiki054-n/TTT) の幾何学的枠組みをペロブスカイト超電導材料に適用し、構造安定性・超電導転移・Tc 予測を単一の原理で記述する試みです。

---

## 核心命題

```
R + I + J = 0                   （三軸ベクトルの総和）
R + I + J + V = 0               （正四面体の第四頂点）
∯ v dΩ = 0                      （OOπ：球面への完全展開）
```

| TTT 概念 | ペロブスカイト超電導への対応 |
|---|---|
| 双極のゼロ（00） | Cooper 対（↑↓スピンペア、運動量・スピン総和 = 0） |
| 9点正八面体系 | BO₆ 正八面体（ペロブスカイトの基本単位） |
| ベクトル総和ゼロ | 超電導状態（完全反磁性・零抵抗） |
| 閾値超過 → 系の破壊 | 脱対（pair-breaking）、Tc での相転移 |
| OOπ（球面展開） | 完全立方対称（Pm3̄m）、最高 Tc 状態 |

---

## リポジトリ構成

```
perovskite-ttt/
├── README.md                        # 本ファイル
├── theory/
│   ├── 01_ttt_principles.md         # TTT 基本原理とペロブスカイトへの適用
│   ├── 02_nine_point_octahedron.md  # 9点正八面体系の詳細
│   ├── 03_cooper_pair_as_00.md      # Cooper 対 = 双極のゼロの量子的実現
│   ├── 04_tc_threshold_theory.md    # Tc と TTT 閾値理論
│   └── 05_structural_instability.md # ヤーン・テラー効果と点の歪み
├── docs/
│   ├── correspondence_table.md      # TTT ↔ 物理学 対応表（詳細版）
│   └── materials_catalog.md        # 代表的ペロブスカイト超電導材料一覧
├── simulation/
│   ├── bo6_octahedron.py            # BO₆ 正八面体の TTT 9点系シミュレーション
│   ├── vector_sum_analyzer.py       # ベクトル総和ゼロの可視化・解析
│   └── tc_predictor.py             # Tc 予測モデル（TTT 閾値ベース）
├── assets/
│   └── ttt-perovskite-demo.html    # インタラクティブデモ
└── LICENSE
```

---

## クイックスタート

```bash
git clone https://github.com/<your-username>/perovskite-ttt
cd perovskite-ttt
pip install numpy matplotlib scipy
python simulation/bo6_octahedron.py
```

---

## 理論的背景

### なぜ TTT がペロブスカイト超電導に有効か

ペロブスカイト（ABO₃）の基本単位である **BO₆ 正八面体**は、TTT の **9点正八面体系**の直接的な実現です：

```
TTT 9点系                BO₆ 正八面体
───────────────────────────────────────
実数軸上の3点（R）   ←→  ±x, ±y, ±z の酸素ペア（3軸）
虚数軸上の4点（I,J） ←→  遷移金属 B + 軸外対称点
中心原点（00）       ←→  B サイトイオン（反転対称中心）
ベクトル総和 = 0     ←→  完全立方対称 Pm3̄m
```

### Cooper 対と双極のゼロ

BCS 理論の Cooper 対：
```
電子A (+k↑) + 電子B (−k↓) → Cooper 対
スピン総和 = 0、運動量総和 = 0   ← TTT「00」の量子的実現
```

超電導転移 = **全電子ペアが「ベクトル総和ゼロ」に収束する相転移**

> 「原点は出発点ではなく、到達点である」— TTT

### Tc と閾値理論

| TTT 安定条件 | 超電導安定条件 |
|---|---|
| `\|Δv\| < v_c`（臨界変位以下） | `k_BT < 2Δ`（熱エネルギー < ペアギャップ） |
| 閾値超過 → 系の崩壊 | `T > Tc` → Cooper 対の破壊 |

---

## 関連リポジトリ

- **TTT（トリテトラ理論）本体**: [github.com/kiki054-n/TTT](https://github.com/kiki054-n/TTT)

---

## ライセンス

MIT License — 詳細は [LICENSE](LICENSE) を参照

---

*ベクトルの総和はゼロになる — 宇宙の釣り合いの原理*  
**by kiki054-n**
