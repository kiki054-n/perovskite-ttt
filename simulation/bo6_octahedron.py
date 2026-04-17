"""
bo6_octahedron.py
==================
BO₆ 正八面体を TTT 9点系として解析するシミュレーション。
ベクトル総和ゼロへの収束と歪み指数（TTT-D）を計算・可視化する。

使い方：
  python bo6_octahedron.py

依存ライブラリ：
  numpy, matplotlib, scipy
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.proj3d import proj_transform
from matplotlib import rcParams

rcParams['font.family'] = 'sans-serif'
rcParams['axes.unicode_minus'] = False


# ────────────────────────────────────────────
# 1. TTT 9点系の定義
# ────────────────────────────────────────────

class BO6_TTT:
    """
    BO₆ 正八面体を TTT 9点正八面体系として表現するクラス。

    TTT 9点系：
      中心原点（00）：B サイト金属イオン + 対称双極点
      R 軸：±x 方向の酸素（O1, O2）
      I 軸：±y 方向の酸素（O3, O4）
      J 軸：±z 方向の酸素（O5, O6 = 頂点酸素）
    """

    def __init__(self, d_x=2.0, d_y=2.0, d_z=2.0, label="Ideal BO₆"):
        """
        Parameters
        ----------
        d_x : float  B-O 結合長（x 方向）[Å]
        d_y : float  B-O 結合長（y 方向）[Å]
        d_z : float  B-O 結合長（z 方向）[Å]
        label : str  材料名または状態名
        """
        self.d_x = d_x
        self.d_y = d_y
        self.d_z = d_z
        self.label = label

        # 9点系の座標定義（TTT マッピング）
        self._build_points()

    def _build_points(self):
        """9点（6酸素 + B + 双極点2個）の座標を構築"""
        # 中心 B サイト（00）
        self.B_center = np.array([0.0, 0.0, 0.0])

        # 6つの酸素イオン（ベクトル）
        self.oxygen = {
            "O_R+": np.array([+self.d_x, 0, 0]),   # R+ 軸
            "O_R-": np.array([-self.d_x, 0, 0]),   # R- 軸
            "O_I+": np.array([0, +self.d_y, 0]),   # I+ 軸
            "O_I-": np.array([0, -self.d_y, 0]),   # I- 軸
            "O_J+": np.array([0, 0, +self.d_z]),   # J+ 軸（頂点）
            "O_J-": np.array([0, 0, -self.d_z]),   # J- 軸（頂点）
        }

    def vector_sum(self):
        """全酸素ベクトルの総和（TTT: = 0 が理想）"""
        return sum(self.oxygen.values())

    def ttt_distortion_index(self):
        """
        TTT 歪み指数 D = |ベクトル総和| の大きさ
        D = 0 → 完全均衡（OOπ 状態）
        """
        return np.linalg.norm(self.vector_sum())

    def jt_parameter(self):
        """
        ヤーン・テラー歪みパラメータ Q_JT
        Q > 0 : z 伸長（oblate）
        Q < 0 : z 圧縮（prolate）
        """
        d0 = (self.d_x + self.d_y + self.d_z) / 3
        Q = (2*self.d_z - self.d_x - self.d_y) / np.sqrt(6)
        return Q

    def equilibrium_score(self):
        """
        均衡度 EQ（0〜1）
        EQ = 1 → 完全均衡（超電導 Tc 最大を期待）
        EQ = 0 → 完全崩壊
        """
        d_mean = (self.d_x + self.d_y + self.d_z) / 3
        D = self.ttt_distortion_index()
        return max(0.0, 1.0 - D / d_mean)

    def predict_tc(self, tc_max=92.0, beta=2.0, tolerance_factor=1.0):
        """
        TTT 理論に基づく Tc 予測（仮説モデル）

        Tc = Tc_max × EQ × exp(−β|t−1|)

        Parameters
        ----------
        tc_max : float  材料クラスの最大 Tc [K]
        beta : float    tolerance factor 感度
        tolerance_factor : float  ゴールドシュミット t 値
        """
        eq = self.equilibrium_score()
        f_t = np.exp(-beta * abs(tolerance_factor - 1.0))
        return tc_max * eq * f_t

    def report(self):
        """解析結果をテキストで出力"""
        vs = self.vector_sum()
        print(f"\n{'='*50}")
        print(f"  {self.label}")
        print(f"{'='*50}")
        print(f"  結合長: d_x={self.d_x:.3f}, d_y={self.d_y:.3f}, d_z={self.d_z:.3f} Å")
        print(f"  ベクトル総和: ({vs[0]:.4f}, {vs[1]:.4f}, {vs[2]:.4f})")
        print(f"  TTT-D（歪み指数）: {self.ttt_distortion_index():.4f}")
        print(f"  JT パラメータ Q: {self.jt_parameter():.4f} Å")
        print(f"  均衡度 EQ: {self.equilibrium_score():.4f}")
        print(f"  予測 Tc（仮説）: {self.predict_tc():.1f} K")


# ────────────────────────────────────────────
# 2. 可視化
# ────────────────────────────────────────────

def plot_bo6_3d(bo6: BO6_TTT, ax=None, show_vectors=True):
    """BO₆ 正八面体の 3D 可視化"""
    if ax is None:
        fig = plt.figure(figsize=(7, 7))
        ax = fig.add_subplot(111, projection='3d')

    colors = {
        "O_R+": "#E24B4A", "O_R-": "#E24B4A",
        "O_I+": "#378ADD", "O_I-": "#378ADD",
        "O_J+": "#1D9E75", "O_J-": "#1D9E75",
    }
    axis_labels = {
        "O_R+": "R+", "O_R-": "R−",
        "O_I+": "I+", "O_I-": "I−",
        "O_J+": "J+", "O_J-": "J−",
    }

    # B サイト（中心）
    ax.scatter([0], [0], [0], s=300, c='#534AB7', zorder=5, label='B（中心）')

    # 酸素イオン
    for name, pos in bo6.oxygen.items():
        c = colors[name]
        ax.scatter(*pos, s=200, c=c, zorder=5)
        ax.text(pos[0]*1.15, pos[1]*1.15, pos[2]*1.15,
                axis_labels[name], fontsize=9, ha='center', color=c)

        # B-O 結合線
        ax.plot([0, pos[0]], [0, pos[1]], [0, pos[2]],
                color=c, linewidth=1.5, alpha=0.7)

    # ベクトル総和
    vs = bo6.vector_sum()
    if np.linalg.norm(vs) > 1e-6 and show_vectors:
        ax.quiver(0, 0, 0, vs[0], vs[1], vs[2],
                  color='#D85A30', linewidth=2, arrow_length_ratio=0.3,
                  label=f'ベクトル総和 |D|={np.linalg.norm(vs):.3f}')

    d_max = max(bo6.d_x, bo6.d_y, bo6.d_z) * 1.3
    ax.set_xlim(-d_max, d_max)
    ax.set_ylim(-d_max, d_max)
    ax.set_zlim(-d_max, d_max)
    ax.set_xlabel('R 軸 (x)', labelpad=5)
    ax.set_ylabel('I 軸 (y)', labelpad=5)
    ax.set_zlabel('J 軸 (z)', labelpad=5)
    ax.set_title(f"{bo6.label}\nTTT-D={bo6.ttt_distortion_index():.4f}, EQ={bo6.equilibrium_score():.4f}",
                 fontsize=11)
    ax.legend(fontsize=8)
    return ax


def plot_distortion_vs_tc(materials):
    """TTT-D（歪み指数）vs 予測 Tc の相関プロット"""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    D_vals = [m.ttt_distortion_index() for m in materials]
    EQ_vals = [m.equilibrium_score() for m in materials]
    Tc_vals = [m.predict_tc() for m in materials]
    labels = [m.label for m in materials]

    # パネル1：TTT-D vs Tc
    ax1 = axes[0]
    colors = plt.cm.plasma(np.linspace(0.2, 0.9, len(materials)))
    for i, (d, tc, lbl) in enumerate(zip(D_vals, Tc_vals, labels)):
        ax1.scatter(d, tc, s=120, color=colors[i], zorder=3, label=lbl)
    ax1.set_xlabel('TTT 歪み指数 D', fontsize=11)
    ax1.set_ylabel('予測 Tc [K]', fontsize=11)
    ax1.set_title('TTT-D vs 予測 Tc', fontsize=12)
    ax1.legend(fontsize=8, loc='upper right')
    ax1.grid(True, alpha=0.3)

    # パネル2：均衡度 EQ vs Tc
    ax2 = axes[1]
    for i, (eq, tc, lbl) in enumerate(zip(EQ_vals, Tc_vals, labels)):
        ax2.scatter(eq, tc, s=120, color=colors[i], zorder=3, label=lbl)
    eq_line = np.linspace(0, 1, 100)
    ax2.plot(eq_line, 92.0 * eq_line, 'k--', alpha=0.4, linewidth=1,
             label='Tc_max × EQ')
    ax2.set_xlabel('均衡度 EQ', fontsize=11)
    ax2.set_ylabel('予測 Tc [K]', fontsize=11)
    ax2.set_title('均衡度 EQ vs 予測 Tc', fontsize=12)
    ax2.legend(fontsize=8, loc='upper left')
    ax2.grid(True, alpha=0.3)

    plt.suptitle('TTT理論によるペロブスカイト超電導安定性解析', fontsize=13)
    plt.tight_layout()
    plt.savefig('ttt_perovskite_analysis.png', dpi=150, bbox_inches='tight')
    print("\n  → ttt_perovskite_analysis.png を保存しました")
    plt.show()


# ────────────────────────────────────────────
# 3. メイン実行
# ────────────────────────────────────────────

def main():
    print("=" * 60)
    print("  BO₆ 正八面体 × TTT 9点系 シミュレーション")
    print("  ベクトルの総和はゼロになる")
    print("=" * 60)

    # 代表的な材料のモデル化
    materials = [
        # 完全正八面体（理想状態）
        BO6_TTT(d_x=2.00, d_y=2.00, d_z=2.00, label="理想 BO₆（OOπ）"),

        # YBCO 的 CuO₂ 面（弱いヤーン・テラー）
        BO6_TTT(d_x=1.94, d_y=1.94, d_z=2.30, label="YBCO 型 CuO₂"),

        # LSCO 最適ドープ
        BO6_TTT(d_x=1.90, d_y=1.90, d_z=2.42, label="LSCO（最適ドープ）"),

        # LaMnO₃（強いヤーン・テラー歪み）
        BO6_TTT(d_x=1.96, d_y=1.96, d_z=2.28, label="LaMnO₃（JT 歪み）"),

        # 強歪み（超電導なし領域）
        BO6_TTT(d_x=1.85, d_y=1.85, d_z=2.50, label="過剰 JT 歪み"),
    ]

    # 各材料の解析レポート
    for m in materials:
        m.report()

    # 可視化
    print("\n  → 3D 可視化を生成中...")

    # 代表的な 2 ケースを 3D プロット
    fig = plt.figure(figsize=(14, 6))
    ax1 = fig.add_subplot(121, projection='3d')
    ax2 = fig.add_subplot(122, projection='3d')
    plot_bo6_3d(materials[0], ax=ax1)   # 理想状態
    plot_bo6_3d(materials[2], ax=ax2)   # LSCO

    plt.suptitle('TTT 9点系：BO₆ 正八面体の構造比較', fontsize=13)
    plt.tight_layout()
    plt.savefig('bo6_3d_comparison.png', dpi=150, bbox_inches='tight')
    print("  → bo6_3d_comparison.png を保存しました")
    plt.show()

    # TTT-D vs Tc 相関プロット
    plot_distortion_vs_tc(materials)

    print("\n  シミュレーション完了。")
    print("  TTT 理論：ベクトルの総和はゼロになる")


if __name__ == "__main__":
    main()
