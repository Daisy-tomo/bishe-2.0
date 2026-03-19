"""
参数敏感性分析
=============
逐一测试法 (One-At-a-Time, OAT)：
  - 固定其他12个参数为真值
  - 依次让每个参数在其先验区间内变化（N_points个均匀采样点）
  - 运行前向模型，计算输出 R_ud 和 C_ud 的变化幅度
  - 以观测噪声水平（1%）为阈值，剔除不敏感参数
  - 保存每一个逐一测试的详细结果
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import json
from datetime import datetime

# ─────────────────────────────────────────────
# 1. 分段函数（与 MATLAB 一致）
# ─────────────────────────────────────────────

def piecewise_kT(T_g):
    """燃气绝热指数"""
    if T_g > 1600:
        return 1.25
    elif T_g > 1400:
        return 1.30
    elif T_g > 800:
        return 1.33
    else:
        return 1.33


def piecewise_RT(T_g):
    """燃气气体常数 [J/(kg·K)]"""
    if T_g > 1600:
        return 288.6
    elif T_g > 1400:
        return 288.0
    elif T_g > 800:
        return 287.6
    else:
        return 287.6


def delta_cooling(T_g):
    """涡轮冷却引气系数"""
    d = 0.02 + (T_g - 1200) / 100 * 0.02
    return max(0.0, min(d, 0.15))


# ─────────────────────────────────────────────
# 2. 前向模型（与 MATLAB engine_forward 完全一致）
# ─────────────────────────────────────────────

def engine_forward(theta, cond):
    """
    前向模型：给定13个参数和工况，返回 (R_ud, C_ud)。
    若计算失败返回 (nan, nan)。

    Parameters
    ----------
    theta : array-like, length 13
    cond  : dict with keys T_H, M_flight, m, pi_k, T_g

    Returns
    -------
    y   : np.ndarray [R_ud, C_ud]
    aux : dict  中间变量
    """
    (eta_k, eta_t, eta_m, eta_v, eta_tv,
     eta_c1, eta_c2,
     sigma_cc, sigma_kan, sigma_kask, sigma_ks,
     eta_T, lam) = theta

    T_H      = cond['T_H']
    M_flight = cond['M_flight']
    m        = cond['m']
    pi_k     = cond['pi_k']
    T_g      = cond['T_g']

    nan_y = np.array([np.nan, np.nan])

    try:
        # (1) 基本常数
        k_air = 1.4
        R_air = 287.3
        a = np.sqrt(k_air * R_air * T_H)
        if not (np.isfinite(a) and a > 0):
            return nan_y, {}
        V_flight = a * M_flight

        # (2) 分段函数值
        kT = piecewise_kT(T_g)
        RT = piecewise_RT(T_g)
        d  = delta_cooling(T_g)

        # (3) 进口总压比与压气机入口温度
        inner = 1 + V_flight**2 / (2 * (k_air / (k_air - 1)) * R_air * T_H)
        if inner <= 0:
            return nan_y, {}
        tau_v = inner ** (k_air / (k_air - 1))
        T_B   = T_H * (inner ** k_air)
        if not (np.isfinite(T_B) and T_B > 0):
            return nan_y, {}

        # (4) 压气机出口温度
        pi_k_ratio = pi_k ** ((k_air - 1) / k_air)
        if not (np.isfinite(pi_k_ratio) and pi_k_ratio >= 1):
            return nan_y, {}
        T_k = T_B * (1 + (pi_k_ratio - 1) / eta_k)
        if not (np.isfinite(T_k) and T_k > 0):
            return nan_y, {}

        # (5) 相对耗油量
        g_T = 3e-5 * T_g - 2.69e-5 * T_k - 0.003
        if not (np.isfinite(g_T) and g_T > 0):
            return nan_y, {}

        # (6) 热恢复系数 lambda_heat
        compress_work = (k_air / (k_air - 1)) * R_air * T_B * (pi_k_ratio - 1)
        gas_enthalpy  = (kT / (kT - 1)) * RT * T_g
        if abs(gas_enthalpy) < 1e-6:
            return nan_y, {}
        num_lambda = 1 - compress_work / (gas_enthalpy * eta_k)
        den_lambda = 1 - compress_work / (gas_enthalpy * eta_k * eta_t)
        if abs(den_lambda) < 1e-10:
            return nan_y, {}
        lambda_heat = num_lambda / den_lambda
        if not np.isfinite(lambda_heat):
            return nan_y, {}

        # (7) 进口总压恢复系数
        sigma_bx = sigma_cc * sigma_kan

        # (8) 单位自由能 L_sv
        exp_T = (kT - 1) / kT
        expansion_pr_denom = tau_v * sigma_bx * pi_k * sigma_kask * sigma_ks
        if expansion_pr_denom <= 0:
            return nan_y, {}
        expansion_term = (1.0 / expansion_pr_denom) ** exp_T
        if not np.isfinite(expansion_term):
            return nan_y, {}
        term1 = (kT / (kT - 1)) * RT * T_g * (1 - expansion_term)
        compress_work2 = (k_air / (k_air - 1)) * R_air * T_B * (pi_k ** ((k_air - 1) / k_air) - 1)
        denom2 = (1 + g_T) * eta_k * eta_T * eta_t * eta_m * (1 - d)
        if abs(denom2) < 1e-10:
            return nan_y, {}
        term2 = compress_work2 / denom2
        L_sv = lambda_heat * (term1 - term2)
        if not (np.isfinite(L_sv) and L_sv > 0):
            return nan_y, {}

        # (9) 最优自由能分配系数 x_pc
        V2_term  = m * V_flight**2
        num_xpc  = 1 + V2_term / (2 * L_sv * eta_tv * eta_v * eta_c2)
        den_xpc  = 1 + (m * eta_tv * eta_v * eta_c2) / (eta_c1 * lam)
        if abs(den_xpc) < 1e-10:
            return nan_y, {}
        x_pc = num_xpc / den_xpc
        if not np.isfinite(x_pc) or x_pc <= 0:
            return nan_y, {}

        # (10) 比推力 R_ud
        inner_sq1 = 2 * eta_c1 * lam * x_pc * L_sv
        if inner_sq1 < 0:
            return nan_y, {}
        V_j1 = (1 + g_T) * np.sqrt(inner_sq1) - V_flight
        inner_sq2 = 2 * (1 - x_pc) / m * L_sv * eta_tv * eta_v * eta_c2 + V_flight**2
        if inner_sq2 < 0:
            return nan_y, {}
        V_j2 = np.sqrt(inner_sq2) - V_flight
        R_ud = (1 / (1 + m)) * V_j1 + (m / (1 + m)) * V_j2
        if not (np.isfinite(R_ud) and R_ud > 0):
            return nan_y, {}

        # (11) 比油耗 C_ud
        denom_C = R_ud * (1 + m)
        if abs(denom_C) < 1e-10:
            return nan_y, {}
        C_ud = 3600 * g_T * (1 - d) / denom_C
        if not (np.isfinite(C_ud) and C_ud > 0):
            return nan_y, {}

        y = np.array([R_ud, C_ud])
        aux = dict(T_B=T_B, T_k=T_k, tau_v=tau_v, g_T=g_T,
                   lambda_heat=lambda_heat, sigma_bx=sigma_bx,
                   L_sv=L_sv, x_pc=x_pc, kT=kT, RT=RT, delta=d)
        return y, aux

    except Exception as e:
        return nan_y, {}


# ─────────────────────────────────────────────
# 3. 参数定义
# ─────────────────────────────────────────────

PARAM_NAMES = [
    'eta_k',      # 压气机绝热效率
    'eta_t',      # 涡轮绝热效率
    'eta_m',      # 机械效率
    'eta_v',      # 风扇效率
    'eta_tv',     # 风扇涡轮效率
    'eta_c1',     # 一次喷管效率
    'eta_c2',     # 二次喷管效率
    'sigma_cc',   # 进气道激波总压恢复系数
    'sigma_kan',  # 进气通道总压恢复系数
    'sigma_kask', # 压气机级间总压恢复系数
    'sigma_ks',   # 燃烧室总压恢复系数
    'eta_T',      # 燃烧放热系数
    'lambda',     # 风扇涡轮热恢复系数
]

PARAM_NAMES_CN = [
    '压气机绝热效率',
    '涡轮绝热效率',
    '机械效率',
    '风扇效率',
    '风扇涡轮效率',
    '一次喷管效率',
    '二次喷管效率',
    '进气道激波总压恢复系数',
    '进气通道总压恢复系数',
    '压气机级间总压恢复系数',
    '燃烧室总压恢复系数',
    '燃烧放热系数',
    '风扇涡轮热恢复系数',
]

LB = np.array([0.84, 0.86, 0.980, 0.85, 0.90, 0.94, 0.92,
               0.98, 0.98, 0.98, 0.94, 0.97, 1.02])
UB = np.array([0.86, 0.92, 0.995, 0.87, 0.92, 0.95, 0.94,
               1.00, 0.99, 0.99, 0.96, 0.99, 1.04])

THETA_TRUE = np.array([
    0.85,   # eta_k
    0.89,   # eta_t
    0.988,  # eta_m
    0.860,  # eta_v
    0.910,  # eta_tv
    0.945,  # eta_c1
    0.930,  # eta_c2
    0.990,  # sigma_cc
    0.985,  # sigma_kan
    0.985,  # sigma_kask
    0.950,  # sigma_ks
    0.980,  # eta_T
    1.030,  # lambda
])

COND = {
    'T_H':      288.0,
    'M_flight': 0.0,
    'm':        10.0,
    'pi_k':     33.0,
    'T_g':      1700.0,
}

# 观测噪声水平（与 MATLAB 一致）
NOISE_LEVEL = 0.01   # 1%
# 敏感性阈值：输出相对变化范围 < threshold 则认为不敏感
SENSITIVITY_THRESHOLD = NOISE_LEVEL  # 1%

N_POINTS = 100   # 每个参数的扫描点数


# ─────────────────────────────────────────────
# 4. 主敏感性分析函数
# ─────────────────────────────────────────────

def run_sensitivity_analysis(theta_true, lb, ub, cond,
                             n_points=N_POINTS,
                             threshold=SENSITIVITY_THRESHOLD):
    """
    逐一测试 (OAT) 敏感性分析。

    Returns
    -------
    results : list of dict  每个参数的详细结果
    y_base  : np.ndarray    基准输出 [R_ud, C_ud]
    """
    n_params = len(theta_true)

    # 基准输出（所有参数均取真值）
    y_base, aux_base = engine_forward(theta_true, cond)
    assert np.all(np.isfinite(y_base)), "真值参数前向模型输出无效！"
    R_base, C_base = y_base

    print(f"基准输出:")
    print(f"  R_ud_base = {R_base:.4f} [N·s/kg]")
    print(f"  C_ud_base = {C_base:.6f} [kg/(N·h)]")
    print()

    results = []

    for i in range(n_params):
        # 在先验区间内均匀取 n_points 个点
        param_values = np.linspace(lb[i], ub[i], n_points)
        R_vals = np.full(n_points, np.nan)
        C_vals = np.full(n_points, np.nan)

        for j, pv in enumerate(param_values):
            theta_test = theta_true.copy()
            theta_test[i] = pv
            y, _ = engine_forward(theta_test, cond)
            if np.all(np.isfinite(y)):
                R_vals[j] = y[0]
                C_vals[j] = y[1]

        # 有效点
        valid = np.isfinite(R_vals) & np.isfinite(C_vals)
        n_valid = valid.sum()

        if n_valid < 2:
            # 无法计算变化幅度
            R_range_rel = np.nan
            C_range_rel = np.nan
            max_rel     = np.nan
            sensitive   = None
        else:
            R_range = np.nanmax(R_vals) - np.nanmin(R_vals)
            C_range = np.nanmax(C_vals) - np.nanmin(C_vals)
            R_range_rel = R_range / abs(R_base)   # 相对于基准的变化范围
            C_range_rel = C_range / abs(C_base)
            max_rel = max(R_range_rel, C_range_rel)
            sensitive = bool(max_rel >= threshold)

        results.append({
            'index':        i,
            'name':         PARAM_NAMES[i],
            'name_cn':      PARAM_NAMES_CN[i],
            'lb':           lb[i],
            'ub':           ub[i],
            'true_value':   theta_true[i],
            'param_values': param_values.tolist(),
            'R_vals':       R_vals.tolist(),
            'C_vals':       C_vals.tolist(),
            'valid_count':  int(n_valid),
            'R_range_rel':  float(R_range_rel) if np.isfinite(R_range_rel) else None,
            'C_range_rel':  float(C_range_rel) if np.isfinite(C_range_rel) else None,
            'max_rel':      float(max_rel)      if np.isfinite(max_rel)    else None,
            'sensitive':    sensitive,
        })

        status = "敏感  ✓" if sensitive else ("不敏感 ✗" if sensitive is False else "未知")
        print(f"[{i+1:02d}/{n_params}] {PARAM_NAMES[i]:12s} ({PARAM_NAMES_CN[i]})")
        print(f"         区间: [{lb[i]:.4f}, {ub[i]:.4f}]   真值: {theta_true[i]:.4f}")
        if n_valid >= 2:
            print(f"         ΔR_ud/R_base = {R_range_rel*100:.4f}%   "
                  f"ΔC_ud/C_base = {C_range_rel*100:.4f}%   "
                  f"max = {max_rel*100:.4f}%  →  {status}")
        else:
            print(f"         有效点不足 ({n_valid})，无法评估")
        print()

    return results, y_base


# ─────────────────────────────────────────────
# 5. 可视化
# ─────────────────────────────────────────────

def plot_oat_curves(results, y_base, out_dir='.'):
    """为每个参数绘制 OAT 扫描曲线，并保存。"""
    R_base, C_base = y_base
    n_params = len(results)

    # ── 5a. 每个参数一张子图（汇总图）
    n_cols = 3
    n_rows = int(np.ceil(n_params / n_cols))

    fig_R, axes_R = plt.subplots(n_rows, n_cols,
                                 figsize=(6 * n_cols, 4 * n_rows))
    fig_C, axes_C = plt.subplots(n_rows, n_cols,
                                 figsize=(6 * n_cols, 4 * n_rows))

    for idx, res in enumerate(results):
        row, col = divmod(idx, n_cols)
        ax_R = axes_R[row, col]
        ax_C = axes_C[row, col]

        pvs   = np.array(res['param_values'])
        R_arr = np.array(res['R_vals'])
        C_arr = np.array(res['C_vals'])

        # R_ud 子图
        ax_R.plot(pvs, R_arr, 'b-', linewidth=1.2)
        ax_R.axvline(res['true_value'], color='g', linestyle='--',
                     linewidth=1.2, label='真值')
        ax_R.axhline(R_base, color='r', linestyle=':', linewidth=1.0,
                     label='基准')
        ax_R.set_xlabel(res['name'])
        ax_R.set_ylabel('R_ud [N·s/kg]')
        title_r = (f"{res['name']} ({res['name_cn']})\n"
                   f"ΔR/R={res['R_range_rel']*100:.3f}%"
                   if res['R_range_rel'] is not None else res['name'])
        color_r = 'red' if res['sensitive'] else 'gray'
        ax_R.set_title(title_r, fontsize=8, color=color_r)
        if idx == 0:
            ax_R.legend(fontsize=7)

        # C_ud 子图
        ax_C.plot(pvs, C_arr, 'm-', linewidth=1.2)
        ax_C.axvline(res['true_value'], color='g', linestyle='--',
                     linewidth=1.2, label='真值')
        ax_C.axhline(C_base, color='r', linestyle=':', linewidth=1.0,
                     label='基准')
        ax_C.set_xlabel(res['name'])
        ax_C.set_ylabel('C_ud [kg/(N·h)]')
        title_c = (f"{res['name']} ({res['name_cn']})\n"
                   f"ΔC/C={res['C_range_rel']*100:.3f}%"
                   if res['C_range_rel'] is not None else res['name'])
        color_c = 'red' if res['sensitive'] else 'gray'
        ax_C.set_title(title_c, fontsize=8, color=color_c)
        if idx == 0:
            ax_C.legend(fontsize=7)

    # 隐藏多余子图
    for idx in range(n_params, n_rows * n_cols):
        row, col = divmod(idx, n_cols)
        axes_R[row, col].set_visible(False)
        axes_C[row, col].set_visible(False)

    fig_R.suptitle('OAT 敏感性分析 — R_ud 扫描曲线\n(红色标题=敏感，灰色=不敏感)',
                   fontsize=12)
    fig_C.suptitle('OAT 敏感性分析 — C_ud 扫描曲线\n(红色标题=敏感，灰色=不敏感)',
                   fontsize=12)

    fig_R.tight_layout()
    fig_C.tight_layout()

    path_R = os.path.join(out_dir, 'sensitivity_R_ud_curves.png')
    path_C = os.path.join(out_dir, 'sensitivity_C_ud_curves.png')
    fig_R.savefig(path_R, dpi=150)
    fig_C.savefig(path_C, dpi=150)
    plt.close(fig_R)
    plt.close(fig_C)
    print(f"已保存: {path_R}")
    print(f"已保存: {path_C}")

    # ── 5b. 敏感性排名条形图
    names    = [r['name'] for r in results]
    max_rels = [r['max_rel'] * 100 if r['max_rel'] is not None else 0.0
                for r in results]
    colors   = ['tomato' if r['sensitive'] else 'steelblue' for r in results]

    sort_idx = np.argsort(max_rels)[::-1]
    names_s    = [names[i]    for i in sort_idx]
    max_rels_s = [max_rels[i] for i in sort_idx]
    colors_s   = [colors[i]   for i in sort_idx]

    fig_bar, ax_bar = plt.subplots(figsize=(10, 6))
    bars = ax_bar.barh(names_s, max_rels_s, color=colors_s, edgecolor='k',
                       linewidth=0.5)
    ax_bar.axvline(SENSITIVITY_THRESHOLD * 100, color='black',
                   linestyle='--', linewidth=1.5,
                   label=f'阈值 {SENSITIVITY_THRESHOLD*100:.0f}%')
    ax_bar.set_xlabel('max(ΔR/R, ΔC/C) [%]')
    ax_bar.set_title('各参数敏感性排名\n(红色=敏感，蓝色=不敏感)')
    ax_bar.legend()

    # 标注数值
    for bar, val in zip(bars, max_rels_s):
        ax_bar.text(val + 0.02, bar.get_y() + bar.get_height() / 2,
                    f'{val:.3f}%', va='center', fontsize=8)

    fig_bar.tight_layout()
    path_bar = os.path.join(out_dir, 'sensitivity_ranking.png')
    fig_bar.savefig(path_bar, dpi=150)
    plt.close(fig_bar)
    print(f"已保存: {path_bar}")

    return path_R, path_C, path_bar


def plot_individual_curves(results, y_base, out_dir='.'):
    """为每个参数单独保存一张详细图（R 和 C 并排）。"""
    R_base, C_base = y_base
    ind_dir = os.path.join(out_dir, 'individual_curves')
    os.makedirs(ind_dir, exist_ok=True)

    for res in results:
        pvs   = np.array(res['param_values'])
        R_arr = np.array(res['R_vals'])
        C_arr = np.array(res['C_vals'])

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))

        # R_ud
        ax1.plot(pvs, R_arr, 'b-o', markersize=2, linewidth=1.2)
        ax1.axvline(res['true_value'], color='g', linestyle='--',
                    linewidth=1.5, label=f"真值={res['true_value']:.4f}")
        ax1.axhline(R_base, color='r', linestyle=':', linewidth=1.2,
                    label=f'基准={R_base:.4f}')
        ax1.set_xlabel(f"{res['name']}  [{res['lb']:.4f}, {res['ub']:.4f}]")
        ax1.set_ylabel('R_ud [N·s/kg]')
        rel_r = f"{res['R_range_rel']*100:.4f}%" if res['R_range_rel'] is not None else 'N/A'
        ax1.set_title(f"R_ud 对 {res['name']} 的敏感性\nΔR/R={rel_r}")
        ax1.legend(fontsize=8)
        ax1.grid(True, alpha=0.3)

        # C_ud
        ax2.plot(pvs, C_arr, 'm-o', markersize=2, linewidth=1.2)
        ax2.axvline(res['true_value'], color='g', linestyle='--',
                    linewidth=1.5, label=f"真值={res['true_value']:.4f}")
        ax2.axhline(C_base, color='r', linestyle=':', linewidth=1.2,
                    label=f'基准={C_base:.6f}')
        ax2.set_xlabel(f"{res['name']}  [{res['lb']:.4f}, {res['ub']:.4f}]")
        ax2.set_ylabel('C_ud [kg/(N·h)]')
        rel_c = f"{res['C_range_rel']*100:.4f}%" if res['C_range_rel'] is not None else 'N/A'
        ax2.set_title(f"C_ud 对 {res['name']} 的敏感性\nΔC/C={rel_c}")
        ax2.legend(fontsize=8)
        ax2.grid(True, alpha=0.3)

        status = '敏感' if res['sensitive'] else '不敏感'
        fig.suptitle(f"参数 {res['index']+1:02d}: {res['name']} ({res['name_cn']})  →  {status}",
                     fontsize=11, color='red' if res['sensitive'] else 'gray')
        fig.tight_layout()

        fname = os.path.join(ind_dir,
                             f"param_{res['index']+1:02d}_{res['name']}.png")
        fig.savefig(fname, dpi=120)
        plt.close(fig)

    print(f"已保存各参数独立图到: {ind_dir}/")
    return ind_dir


# ─────────────────────────────────────────────
# 6. 结果汇总与保存
# ─────────────────────────────────────────────

def print_summary(results, threshold):
    """打印敏感性分析汇总表。"""
    sensitive_params   = [r for r in results if r['sensitive'] is True]
    insensitive_params = [r for r in results if r['sensitive'] is False]

    print("=" * 80)
    print("敏感性分析汇总")
    print("=" * 80)
    print(f"{'#':>3}  {'参数名':12s}  {'中文名':18s}  "
          f"{'ΔR/R [%]':>10}  {'ΔC/C [%]':>10}  {'最大 [%]':>10}  {'结论':6s}")
    print("-" * 80)

    # 按 max_rel 降序排列
    sorted_results = sorted(results,
                            key=lambda r: r['max_rel'] if r['max_rel'] is not None else 0,
                            reverse=True)

    for r in sorted_results:
        mr  = f"{r['max_rel']*100:.4f}" if r['max_rel'] is not None else 'N/A'
        rr  = f"{r['R_range_rel']*100:.4f}" if r['R_range_rel'] is not None else 'N/A'
        cr  = f"{r['C_range_rel']*100:.4f}" if r['C_range_rel'] is not None else 'N/A'
        tag = '敏感  ✓' if r['sensitive'] else ('不敏感 ✗' if r['sensitive'] is False else '未知')
        print(f"{r['index']+1:>3}  {r['name']:12s}  {r['name_cn']:18s}  "
              f"{rr:>10}  {cr:>10}  {mr:>10}  {tag}")

    print("=" * 80)
    print(f"\n[结果]  阈值 = {threshold*100:.1f}%（与观测噪声水平一致）")
    print(f"  敏感参数   ({len(sensitive_params)}个): "
          + ", ".join(r['name'] for r in sensitive_params))
    print(f"  不敏感参数 ({len(insensitive_params)}个): "
          + ", ".join(r['name'] for r in insensitive_params))
    print()
    print("建议：")
    print("  ✓ 保留敏感参数作为 MCMC 待估变量")
    print("  ✗ 将不敏感参数固定为先验均值（或物理合理值），减少参数维度")
    print()


def save_results_json(results, y_base, out_dir='.'):
    """将完整结果保存为 JSON 文件。"""
    output = {
        'timestamp':          datetime.now().isoformat(),
        'n_points_per_param': N_POINTS,
        'threshold_pct':      SENSITIVITY_THRESHOLD * 100,
        'noise_level_pct':    NOISE_LEVEL * 100,
        'y_base': {
            'R_ud': float(y_base[0]),
            'C_ud': float(y_base[1]),
        },
        'cond': COND,
        'parameters': results,
    }
    path = os.path.join(out_dir, 'sensitivity_results.json')
    with open(path, 'w', encoding='utf-8') as f:
        json.dump(output, f, ensure_ascii=False, indent=2)
    print(f"已保存完整结果: {path}")
    return path


def save_results_csv(results, out_dir='.'):
    """将汇总结果保存为 CSV 文件。"""
    path = os.path.join(out_dir, 'sensitivity_summary.csv')
    with open(path, 'w', encoding='utf-8-sig') as f:
        f.write("序号,参数名,中文名,下界,上界,真值,"
                "ΔR/R[%],ΔC/C[%],max_rel[%],敏感性\n")
        for r in sorted(results,
                        key=lambda x: x['max_rel'] if x['max_rel'] else 0,
                        reverse=True):
            rr  = f"{r['R_range_rel']*100:.6f}" if r['R_range_rel'] is not None else ''
            cr  = f"{r['C_range_rel']*100:.6f}" if r['C_range_rel'] is not None else ''
            mr  = f"{r['max_rel']*100:.6f}"     if r['max_rel']     is not None else ''
            tag = '敏感' if r['sensitive'] else ('不敏感' if r['sensitive'] is False else '未知')
            f.write(f"{r['index']+1},{r['name']},{r['name_cn']},"
                    f"{r['lb']},{r['ub']},{r['true_value']},"
                    f"{rr},{cr},{mr},{tag}\n")
    print(f"已保存汇总 CSV: {path}")
    return path


# ─────────────────────────────────────────────
# 7. 主入口
# ─────────────────────────────────────────────

def main():
    out_dir = 'sensitivity_output'
    os.makedirs(out_dir, exist_ok=True)

    print("=" * 80)
    print("涡扇发动机参数敏感性分析（逐一测试 OAT）")
    print("=" * 80)
    print(f"扫描点数/参数: {N_POINTS}")
    print(f"敏感性阈值:    {SENSITIVITY_THRESHOLD*100:.1f}%（观测噪声水平）")
    print(f"工况: T_H={COND['T_H']}K, M={COND['M_flight']}, "
          f"m={COND['m']}, pi_k={COND['pi_k']}, T_g={COND['T_g']}K")
    print()

    # 运行分析
    results, y_base = run_sensitivity_analysis(
        THETA_TRUE, LB, UB, COND,
        n_points=N_POINTS,
        threshold=SENSITIVITY_THRESHOLD
    )

    # 打印汇总
    print_summary(results, SENSITIVITY_THRESHOLD)

    # 保存数据
    save_results_json(results, y_base, out_dir)
    save_results_csv(results, out_dir)

    # 绘图
    print()
    print("正在绘制图形...")
    plot_oat_curves(results, y_base, out_dir)
    plot_individual_curves(results, y_base, out_dir)

    print()
    print("敏感性分析完成。所有结果保存在:", out_dir)


if __name__ == '__main__':
    main()
