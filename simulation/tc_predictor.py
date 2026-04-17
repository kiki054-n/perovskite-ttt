"""
tc_predictor.py
================
TTT 閾値理論に基づく Tc 予測モデル。
BO₆ 歪みパラメータ・ドーピング・tolerance factor を入力として
Tc を推定する。

使い方：
  python tc_predictor.py
"""

import numpy as np
import matplotlib.pyplot as plt

# ────────────────────────────────────────────
# TTT Tc 予測モデル
# ────────────────────────────────────────────

def ttt_tc_model(
    eq: float,
    tolerance_factor: float,
    doping: float,
    tc_max: float = 92.0,
    beta: float = 2.0,
    x_opt: float = 0.15,
    sigma: float = 0.10,
) -> float:
    """
    TTT 理論に基づく Tc 予測式

    Tc = Tc_max × EQ × f(t) × g(x)

    Parameters
    ----------
    eq               : float  BO₆ 均衡度（0〜1、TTT 歪み指数から算出）
    tolerance_factor : float  ゴールドシュミット t 値
    doping           : float  キャリアドーピング量 x
    tc_max           : float  材料クラスの最大 Tc [K]
    beta             : float  tolerance factor 感度係数
    x_opt            : float  最適ドーピング量（Tc 最大点）
    sigma            : float  ドーピング依存性の幅

    Returns
    -------
    Tc : float  予測超電導転移温度 [K]
    """
    # tolerance factor 補正（OOπ からの逸脱ペナルティ）
    f_t = np.exp(-beta * abs(tolerance_factor - 1.0))

    # ドーピング依存性（ガウス型：最適ドープで最大）
    g_x = np.exp(-((doping - x_opt) ** 2) / (2 * sigma ** 2))

    return tc_max * eq * f_t * g_x


def plot_doping_phase_diagram():
    """
    ドーピング x に対する Tc の相図
    LSCO 系（La₁₋ₓSrₓCuO₄）を模擬
    """
    x_vals = np.linspace(0, 0.35, 300)

    # 異なる EQ（BO₆ 均衡度）での比較
    eq_scenarios = [
        (1.00, "EQ=1.00（理想正八面体）", "#1D9E75"),
        (0.92, "EQ=0.92（YBCO 型）", "#378ADD"),
        (0.82, "EQ=0.82（LSCO 最適）", "#7F77DD"),
        (0.65, "EQ=0.65（強歪み）", "#D85A30"),
    ]

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    ax1 = axes[0]
    for eq, label, color in eq_scenarios:
        tc_vals = [ttt_tc_model(eq=eq, tolerance_factor=0.90, doping=x,
                                tc_max=92.0) for x in x_vals]
        ax1.plot(x_vals, tc_vals, color=color, linewidth=2, label=label)

    ax1.set_xlabel('ドーピング量 x', fontsize=11)
    ax1.set_ylabel('予測 Tc [K]', fontsize=11)
    ax1.set_title('ドーピング依存性（均衡度 EQ 別）', fontsize=12)
    ax1.legend(fontsize=9)
    ax1.grid(True, alpha=0.3)
    ax1.axvline(0.15, color='gray', linestyle='--', alpha=0.5, linewidth=1,
                label='最適ドープ x=0.15')

    # tolerance factor 依存性
    ax2 = axes[1]
    t_vals = np.linspace(0.75, 1.15, 300)
    for eq, label, color in eq_scenarios:
        tc_vals = [ttt_tc_model(eq=eq, tolerance_factor=t, doping=0.15,
                                tc_max=92.0) for t in t_vals]
        ax2.plot(t_vals, tc_vals, color=color, linewidth=2, label=label)

    ax2.axvline(1.0, color='gray', linestyle='--', alpha=0.5, linewidth=1)
    ax2.text(1.01, 5, 'OOπ\n(t=1)', fontsize=8, color='gray')
    ax2.set_xlabel('Tolerance factor t', fontsize=11)
    ax2.set_ylabel('予測 Tc [K]', fontsize=11)
    ax2.set_title('Tolerance factor 依存性（EQ 別）', fontsize=12)
    ax2.legend(fontsize=9)
    ax2.grid(True, alpha=0.3)

    plt.suptitle('TTT 閾値理論による Tc 予測モデル', fontsize=13)
    plt.tight_layout()
    plt.savefig('tc_prediction_phase_diagram.png', dpi=150, bbox_inches='tight')
    print("  → tc_prediction_phase_diagram.png を保存しました")
    plt.show()


def plot_ttt_vector_convergence():
    """
    TTT ベクトル総和ゼロへの収束過程の可視化
    超電導転移（T → 0）に対応する
    """
    T = np.linspace(0, 150, 500)  # 温度 [K]
    Tc = 92.0                      # YBCO の Tc
    Delta0 = 30e-3                 # ペアギャップ（eV → kB 単位）

    # BCS ギャップ関数の近似
    Delta_T = Delta0 * np.where(T < Tc, np.sqrt(1 - (T / Tc) ** 3), 0)

    # TTT ベクトル収束度（正規化）
    # T → 0 で「00」に完全収束
    convergence = 1 - (T / Tc) ** 2
    convergence = np.clip(convergence, 0, 1)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 7), sharex=True)

    # ペアギャップ
    ax1.plot(T, Delta_T / Delta0, color='#378ADD', linewidth=2)
    ax1.fill_between(T, 0, Delta_T / Delta0, alpha=0.15, color='#378ADD')
    ax1.axvline(Tc, color='#E24B4A', linestyle='--', linewidth=1.5)
    ax1.text(Tc + 2, 0.8, f'Tc = {Tc} K\n（閾値）', color='#E24B4A', fontsize=9)
    ax1.set_ylabel('正規化ペアギャップ Δ(T)/Δ₀', fontsize=10)
    ax1.set_title('温度とTTT収束度の関係（YBCO 型）', fontsize=12)
    ax1.grid(True, alpha=0.3)
    ax1.set_ylim(-0.05, 1.15)

    # TTT ベクトル収束度
    ax2.plot(T, convergence * 100, color='#1D9E75', linewidth=2)
    ax2.fill_between(T, 0, convergence * 100, alpha=0.15, color='#1D9E75')
    ax2.axvline(Tc, color='#E24B4A', linestyle='--', linewidth=1.5)
    ax2.axhline(100, color='#534AB7', linestyle=':', linewidth=1, alpha=0.7)
    ax2.text(5, 95, '完全収束（OOπ）', color='#534AB7', fontsize=9)
    ax2.set_xlabel('温度 T [K]', fontsize=10)
    ax2.set_ylabel('TTT 収束度（Cooper 対凝縮率 %）', fontsize=10)
    ax2.grid(True, alpha=0.3)
    ax2.set_ylim(-5, 115)

    plt.tight_layout()
    plt.savefig('ttt_vector_convergence.png', dpi=150, bbox_inches='tight')
    print("  → ttt_vector_convergence.png を保存しました")
    plt.show()


# ────────────────────────────────────────────
# メイン実行
# ────────────────────────────────────────────

if __name__ == "__main__":
    print("=" * 60)
    print("  TTT 閾値理論 — Tc 予測モデル")
    print("=" * 60)

    # ドーピング相図
    print("\n  ドーピング相図を生成中...")
    plot_doping_phase_diagram()

    # ベクトル収束過程
    print("\n  TTT ベクトル収束過程を生成中...")
    plot_ttt_vector_convergence()

    # 数値例
    print("\n  代表的な条件での Tc 予測（数値）:")
    print(f"  {'材料':25s}  {'EQ':5s}  {'t':5s}  {'x':5s}  {'Tc[K]':7s}")
    print("  " + "-" * 57)

    cases = [
        ("理想 BO₆ + 最適ドープ", 1.00, 1.00, 0.15),
        ("YBCO 型",              0.92, 0.90, 0.15),
        ("LSCO（最適）",         0.82, 0.90, 0.15),
        ("アンドープ（x=0）",    0.82, 0.90, 0.00),
        ("過剰ドープ（x=0.30）", 0.82, 0.90, 0.30),
        ("強 JT 歪み",           0.65, 0.85, 0.15),
    ]
    for name, eq, t, x in cases:
        tc = ttt_tc_model(eq=eq, tolerance_factor=t, doping=x)
        print(f"  {name:25s}  {eq:.2f}   {t:.2f}   {x:.2f}   {tc:6.1f}")
