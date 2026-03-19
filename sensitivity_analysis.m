%% sensitivity_analysis.m
% 参数敏感性分析 — 逐一测试法 (One-At-a-Time, OAT)
%
% 方法：
%   固定其余12个参数为真值，依次让每个参数在其先验区间内均匀扫描
%   N_POINTS 个点，运行前向模型，计算输出 (R_ud, C_ud) 的相对变化幅度。
%   以观测噪声水平（1%）为阈值：
%     max(ΔR/R, ΔC/C) >= 1%  → 敏感，保留为待估变量
%     max(ΔR/R, ΔC/C) <  1%  → 不敏感，可固定为先验均值
%
% 输出：
%   - 控制台逐参数详细结果 + 汇总表
%   - 每个参数独立扫描图（保存到 sensitivity_output/）
%   - R_ud / C_ud 总览图
%   - 敏感性排名条形图
%   - sensitivity_summary.csv

clear; clc; close all;

%% ──────────────────────────────────────────────────────────────
%  1. 参数定义（与 xiuzhengcanshu.m 完全一致）
%% ──────────────────────────────────────────────────────────────
param_names = {
    'eta\_k',      ...  % 压气机绝热效率
    'eta\_t',      ...  % 涡轮绝热效率
    'eta\_m',      ...  % 机械效率
    'eta\_v',      ...  % 风扇效率
    'eta\_tv',     ...  % 风扇涡轮效率
    'eta\_c1',     ...  % 一次喷管效率
    'eta\_c2',     ...  % 二次喷管效率
    'sigma\_cc',   ...  % 进气道激波总压恢复系数
    'sigma\_kan',  ...  % 进气通道总压恢复系数
    'sigma\_kask', ...  % 压气机级间总压恢复系数
    'sigma\_ks',   ...  % 燃烧室总压恢复系数
    'eta\_T',      ...  % 燃烧放热系数
    'lambda'        ...  % 风扇涡轮热恢复系数
};

param_names_cn = {
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
};

lb = [0.84, 0.86, 0.980, 0.85, 0.90, 0.94, 0.92, 0.98, 0.98, 0.98, 0.94, 0.97, 1.02];
ub = [0.86, 0.92, 0.995, 0.87, 0.92, 0.95, 0.94, 1.00, 0.99, 0.99, 0.96, 0.99, 1.04];

theta_true = [
    0.85,   ...  % eta_k
    0.89,   ...  % eta_t
    0.988,  ...  % eta_m
    0.860,  ...  % eta_v
    0.910,  ...  % eta_tv
    0.945,  ...  % eta_c1
    0.930,  ...  % eta_c2
    0.990,  ...  % sigma_cc
    0.985,  ...  % sigma_kan
    0.985,  ...  % sigma_kask
    0.950,  ...  % sigma_ks
    0.980,  ...  % eta_T
    1.030   ...  % lambda
];

n_params = length(lb);

%% ──────────────────────────────────────────────────────────────
%  2. 工况（与 xiuzhengcanshu.m 完全一致）
%% ──────────────────────────────────────────────────────────────
cond.T_H      = 288.0;
cond.M_flight = 0.0;
cond.m        = 10.0;
cond.pi_k     = 33.0;
cond.T_g      = 1700.0;

%% ──────────────────────────────────────────────────────────────
%  3. 分析配置
%% ──────────────────────────────────────────────────────────────
N_POINTS  = 100;    % 每个参数的扫描点数
THRESHOLD = 0.01;   % 敏感性阈值（= 观测噪声水平 1%）

out_dir = 'sensitivity_output';
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
ind_dir = fullfile(out_dir, 'individual_curves');
if ~exist(ind_dir, 'dir'), mkdir(ind_dir); end

%% ──────────────────────────────────────────────────────────────
%  4. 基准输出
%% ──────────────────────────────────────────────────────────────
[y_base, aux_base] = engine_forward(theta_true, cond);
if ~all(isfinite(y_base))
    error('真值参数前向模型输出无效，请检查 engine_forward。');
end
R_base = y_base(1);
C_base = y_base(2);

fprintf('========================================\n');
fprintf('涡扇发动机参数敏感性分析（逐一测试 OAT）\n');
fprintf('========================================\n');
fprintf('扫描点数/参数: %d\n', N_POINTS);
fprintf('敏感性阈值:    %.1f%%（观测噪声水平）\n', THRESHOLD * 100);
fprintf('工况: T_H=%.0fK, M=%.1f, m=%.0f, pi_k=%.0f, T_g=%.0fK\n', ...
    cond.T_H, cond.M_flight, cond.m, cond.pi_k, cond.T_g);
fprintf('\n基准输出（theta_true）:\n');
fprintf('  R_ud_base = %.4f [N·s/kg]\n', R_base);
fprintf('  C_ud_base = %.6f [kg/(N·h)]\n', C_base);
fprintf('\n');

%% ──────────────────────────────────────────────────────────────
%  5. 逐一测试 OAT
%% ──────────────────────────────────────────────────────────────

% 存储每个参数的结果
SA.param_values = cell(n_params, 1);   % 扫描点
SA.R_vals       = cell(n_params, 1);   % R_ud 值
SA.C_vals       = cell(n_params, 1);   % C_ud 值
SA.R_range_rel  = zeros(n_params, 1);  % ΔR/R
SA.C_range_rel  = zeros(n_params, 1);  % ΔC/C
SA.max_rel      = zeros(n_params, 1);  % max(ΔR/R, ΔC/C)
SA.sensitive    = false(n_params, 1);  % 是否敏感

for i = 1:n_params
    pv = linspace(lb(i), ub(i), N_POINTS);
    R_arr = NaN(1, N_POINTS);
    C_arr = NaN(1, N_POINTS);

    for j = 1:N_POINTS
        theta_test    = theta_true;
        theta_test(i) = pv(j);
        [y_j, ~] = engine_forward(theta_test, cond);
        if all(isfinite(y_j))
            R_arr(j) = y_j(1);
            C_arr(j) = y_j(2);
        end
    end

    SA.param_values{i} = pv;
    SA.R_vals{i}       = R_arr;
    SA.C_vals{i}       = C_arr;

    valid = isfinite(R_arr) & isfinite(C_arr);
    n_valid = sum(valid);

    if n_valid >= 2
        R_range = max(R_arr(valid)) - min(R_arr(valid));
        C_range = max(C_arr(valid)) - min(C_arr(valid));
        SA.R_range_rel(i) = R_range / abs(R_base);
        SA.C_range_rel(i) = C_range / abs(C_base);
        SA.max_rel(i)     = max(SA.R_range_rel(i), SA.C_range_rel(i));
        SA.sensitive(i)   = SA.max_rel(i) >= THRESHOLD;
    else
        SA.R_range_rel(i) = NaN;
        SA.C_range_rel(i) = NaN;
        SA.max_rel(i)     = NaN;
        SA.sensitive(i)   = false;
    end

    % 控制台打印
    if SA.sensitive(i)
        status_str = '敏感  [✓]';
    else
        status_str = '不敏感 [✗]';
    end
    % 去掉下划线转义以便 fprintf
    raw_name = strrep(param_names{i}, '\_', '_');
    fprintf('[%02d/%02d] %-12s (%s)\n', i, n_params, raw_name, param_names_cn{i});
    fprintf('         区间: [%.4f, %.4f]   真值: %.4f\n', lb(i), ub(i), theta_true(i));
    if n_valid >= 2
        fprintf('         ΔR_ud/R_base = %.4f%%   ΔC_ud/C_base = %.4f%%   max = %.4f%%  →  %s\n', ...
            SA.R_range_rel(i)*100, SA.C_range_rel(i)*100, SA.max_rel(i)*100, status_str);
    else
        fprintf('         有效点不足 (%d)，无法评估\n', n_valid);
    end
    fprintf('\n');
end

%% ──────────────────────────────────────────────────────────────
%  6. 汇总表（按 max_rel 降序）
%% ──────────────────────────────────────────────────────────────
[~, sort_idx] = sort(SA.max_rel, 'descend');

fprintf('================================================================================\n');
fprintf('敏感性分析汇总（按影响大小降序排列）\n');
fprintf('================================================================================\n');
fprintf('%-4s  %-12s  %-20s  %10s  %10s  %10s  %s\n', ...
    '#', '参数名', '中文名', 'ΔR/R [%]', 'ΔC/C [%]', '最大 [%]', '结论');
fprintf('%s\n', repmat('-', 1, 80));

for k = 1:n_params
    i = sort_idx(k);
    raw_name = strrep(param_names{i}, '\_', '_');
    if SA.sensitive(i)
        tag = '敏感  [✓]';
    else
        tag = '不敏感 [✗]';
    end
    fprintf('%-4d  %-12s  %-20s  %10.4f  %10.4f  %10.4f  %s\n', ...
        i, raw_name, param_names_cn{i}, ...
        SA.R_range_rel(i)*100, SA.C_range_rel(i)*100, SA.max_rel(i)*100, tag);
end

fprintf('================================================================================\n');
sensitive_list   = find(SA.sensitive);
insensitive_list = find(~SA.sensitive);

fprintf('\n[结果]  阈值 = %.1f%%（观测噪声水平）\n', THRESHOLD*100);
fprintf('  敏感参数   (%d 个):', numel(sensitive_list));
for k = sensitive_list
    fprintf(' %s', strrep(param_names{k}, '\_', '_'));
end
fprintf('\n');
fprintf('  不敏感参数 (%d 个):', numel(insensitive_list));
for k = insensitive_list
    fprintf(' %s', strrep(param_names{k}, '\_', '_'));
end
fprintf('\n\n');
fprintf('建议：\n');
fprintf('  [✓] 保留敏感参数作为 MCMC 待估变量\n');
fprintf('  [✗] 将不敏感参数固定为先验均值，减少参数维度\n\n');

%% ──────────────────────────────────────────────────────────────
%  7. CSV 输出（完整结果）
%% ──────────────────────────────────────────────────────────────
csv_path = fullfile(out_dir, 'sensitivity_summary.csv');
fid = fopen(csv_path, 'w', 'n', 'UTF-8');
fprintf(fid, '序号,参数名,中文名,下界,上界,真值,ΔR/R[%%],ΔC/C[%%],max_rel[%%],敏感性\n');
for k = 1:n_params
    i = sort_idx(k);
    raw_name = strrep(param_names{i}, '\_', '_');
    if SA.sensitive(i), tag2 = '敏感'; else, tag2 = '不敏感'; end
    fprintf(fid, '%d,%s,%s,%.4f,%.4f,%.4f,%.6f,%.6f,%.6f,%s\n', ...
        i, raw_name, param_names_cn{i}, ...
        lb(i), ub(i), theta_true(i), ...
        SA.R_range_rel(i)*100, SA.C_range_rel(i)*100, SA.max_rel(i)*100, tag2);
end
fclose(fid);
fprintf('已保存汇总 CSV: %s\n', csv_path);

% 同时保存各参数逐点扫描数据（每个参数一个 CSV）
for i = 1:n_params
    raw_name = strrep(param_names{i}, '\_', '_');
    pdata_path = fullfile(out_dir, sprintf('scan_param_%02d_%s.csv', i, raw_name));
    fid2 = fopen(pdata_path, 'w', 'n', 'UTF-8');
    fprintf(fid2, '参数值,R_ud,C_ud\n');
    pv    = SA.param_values{i};
    R_arr = SA.R_vals{i};
    C_arr = SA.C_vals{i};
    for j = 1:N_POINTS
        fprintf(fid2, '%.8f,%.8f,%.10f\n', pv(j), R_arr(j), C_arr(j));
    end
    fclose(fid2);
end
fprintf('已保存各参数逐点扫描 CSV（%d 个）到: %s/\n', n_params, out_dir);

%% ──────────────────────────────────────────────────────────────
%  8. 图形输出
%% ──────────────────────────────────────────────────────────────
fprintf('\n正在绘制图形...\n');

% ── 8a. 每个参数独立图（R 和 C 并排）
for i = 1:n_params
    raw_name = strrep(param_names{i}, '\_', '_');
    pv    = SA.param_values{i};
    R_arr = SA.R_vals{i};
    C_arr = SA.C_vals{i};

    fig = figure('Visible', 'off', 'Position', [100 100 1100 420]);

    % R_ud
    subplot(1, 2, 1);
    plot(pv, R_arr, 'b-o', 'MarkerSize', 3, 'LineWidth', 1.2); hold on;
    xline(theta_true(i), 'g--', 'LineWidth', 1.5, ...
        'Label', sprintf('真值=%.4f', theta_true(i)), 'LabelVerticalAlignment', 'bottom');
    yline(R_base, 'r:', 'LineWidth', 1.2, ...
        'Label', sprintf('基准=%.4f', R_base), 'LabelHorizontalAlignment', 'right');
    xlabel(sprintf('%s  [%.4f, %.4f]', raw_name, lb(i), ub(i)));
    ylabel('R_{ud} [N·s/kg]');
    grid on;
    rel_r = SA.R_range_rel(i) * 100;
    title(sprintf('R_{ud} 对 %s 的敏感性\n\\DeltaR/R = %.4f%%', raw_name, rel_r));

    % C_ud
    subplot(1, 2, 2);
    plot(pv, C_arr, 'm-o', 'MarkerSize', 3, 'LineWidth', 1.2); hold on;
    xline(theta_true(i), 'g--', 'LineWidth', 1.5, ...
        'Label', sprintf('真值=%.4f', theta_true(i)), 'LabelVerticalAlignment', 'bottom');
    yline(C_base, 'r:', 'LineWidth', 1.2, ...
        'Label', sprintf('基准=%.6f', C_base), 'LabelHorizontalAlignment', 'right');
    xlabel(sprintf('%s  [%.4f, %.4f]', raw_name, lb(i), ub(i)));
    ylabel('C_{ud} [kg/(N·h)]');
    grid on;
    rel_c = SA.C_range_rel(i) * 100;
    title(sprintf('C_{ud} 对 %s 的敏感性\n\\DeltaC/C = %.4f%%', raw_name, rel_c));

    if SA.sensitive(i)
        status_title = '敏感 [✓]';
        clr = [0.8 0 0];
    else
        status_title = '不敏感 [✗]';
        clr = [0.4 0.4 0.4];
    end
    sgtitle(sprintf('参数 %02d: %s (%s)   →   %s', ...
        i, raw_name, param_names_cn{i}, status_title), ...
        'Color', clr, 'FontSize', 11, 'FontWeight', 'bold');

    fname = fullfile(ind_dir, sprintf('param_%02d_%s.png', i, raw_name));
    saveas(fig, fname);
    close(fig);
end
fprintf('已保存各参数独立图（%d 张）到: %s/\n', n_params, ind_dir);

% ── 8b. R_ud 总览图（4行×4列，按参数序号）
n_cols = 4;
n_rows = ceil(n_params / n_cols);

fig_R = figure('Visible', 'off', 'Name', 'OAT-R_ud', ...
               'Position', [50 50 1600 900]);
for i = 1:n_params
    subplot(n_rows, n_cols, i);
    pv    = SA.param_values{i};
    R_arr = SA.R_vals{i};
    raw_name = strrep(param_names{i}, '\_', '_');

    plot(pv, R_arr, 'b-', 'LineWidth', 1.2); hold on;
    xline(theta_true(i), 'g--', 'LineWidth', 1.2);
    yline(R_base, 'r:',  'LineWidth', 1.0);
    xlabel(raw_name, 'FontSize', 8);
    ylabel('R_{ud}', 'FontSize', 8);
    grid on; box on;

    rel_r = SA.R_range_rel(i) * 100;
    if SA.sensitive(i)
        ttl_clr = [0.8 0 0];
    else
        ttl_clr = [0.4 0.4 0.4];
    end
    title(sprintf('%s\n\\DeltaR/R=%.3f%%', raw_name, rel_r), ...
        'FontSize', 7.5, 'Color', ttl_clr);
end
sgtitle('OAT 敏感性分析 — R_{ud} 扫描曲线（红色=敏感，灰色=不敏感）', 'FontSize', 11);
saveas(fig_R, fullfile(out_dir, 'sensitivity_R_ud_overview.png'));
close(fig_R);
fprintf('已保存: %s\n', fullfile(out_dir, 'sensitivity_R_ud_overview.png'));

% ── 8c. C_ud 总览图
fig_C = figure('Visible', 'off', 'Name', 'OAT-C_ud', ...
               'Position', [50 50 1600 900]);
for i = 1:n_params
    subplot(n_rows, n_cols, i);
    pv    = SA.param_values{i};
    C_arr = SA.C_vals{i};
    raw_name = strrep(param_names{i}, '\_', '_');

    plot(pv, C_arr, 'm-', 'LineWidth', 1.2); hold on;
    xline(theta_true(i), 'g--', 'LineWidth', 1.2);
    yline(C_base, 'r:',  'LineWidth', 1.0);
    xlabel(raw_name, 'FontSize', 8);
    ylabel('C_{ud}', 'FontSize', 8);
    grid on; box on;

    rel_c = SA.C_range_rel(i) * 100;
    if SA.sensitive(i)
        ttl_clr = [0.8 0 0];
    else
        ttl_clr = [0.4 0.4 0.4];
    end
    title(sprintf('%s\n\\DeltaC/C=%.3f%%', raw_name, rel_c), ...
        'FontSize', 7.5, 'Color', ttl_clr);
end
sgtitle('OAT 敏感性分析 — C_{ud} 扫描曲线（红色=敏感，灰色=不敏感）', 'FontSize', 11);
saveas(fig_C, fullfile(out_dir, 'sensitivity_C_ud_overview.png'));
close(fig_C);
fprintf('已保存: %s\n', fullfile(out_dir, 'sensitivity_C_ud_overview.png'));

% ── 8d. 敏感性排名条形图
max_rels_pct = SA.max_rel * 100;
[sorted_mr, sidx] = sort(max_rels_pct, 'descend');
sorted_names = cellfun(@(s) strrep(s,'\_','_'), param_names(sidx), 'UniformOutput', false);
bar_colors = zeros(n_params, 3);
for k = 1:n_params
    if SA.sensitive(sidx(k))
        bar_colors(k,:) = [0.85 0.33 0.10];   % 红橙（敏感）
    else
        bar_colors(k,:) = [0.30 0.55 0.80];   % 蓝（不敏感）
    end
end

fig_bar = figure('Visible', 'off', 'Position', [100 100 900 560]);
hb = barh(1:n_params, sorted_mr, 'FaceColor', 'flat');
hb.CData = bar_colors;
hold on;
xline(THRESHOLD * 100, 'k--', 'LineWidth', 2.0, ...
    'Label', sprintf('阈值 %.0f%%', THRESHOLD*100), ...
    'LabelVerticalAlignment', 'bottom', 'FontSize', 9);
set(gca, 'YTick', 1:n_params, 'YTickLabel', sorted_names, 'FontSize', 9);
xlabel('max(\DeltaR/R, \DeltaC/C) [%]', 'FontSize', 10);
title('各参数敏感性排名（橙红=敏感，蓝=不敏感）', 'FontSize', 11);
% 标注数值
for k = 1:n_params
    text(sorted_mr(k) + 0.01, k, sprintf('%.3f%%', sorted_mr(k)), ...
        'VerticalAlignment', 'middle', 'FontSize', 8);
end
grid on; box on;
saveas(fig_bar, fullfile(out_dir, 'sensitivity_ranking.png'));
close(fig_bar);
fprintf('已保存: %s\n', fullfile(out_dir, 'sensitivity_ranking.png'));

fprintf('\n========================================\n');
fprintf('敏感性分析完成。全部结果保存在: %s/\n', out_dir);
fprintf('========================================\n');

%% ══════════════════════════════════════════════════════════════
%  本地函数（与 xiuzhengcanshu.m 完全一致）
%% ══════════════════════════════════════════════════════════════

function kT = piecewise_kT(T_g)
if T_g > 1600
    kT = 1.25;
elseif T_g > 1400
    kT = 1.30;
elseif T_g > 800
    kT = 1.33;
else
    kT = 1.33;
end
end

function RT = piecewise_RT(T_g)
if T_g > 1600
    RT = 288.6;
elseif T_g > 1400
    RT = 288.0;
elseif T_g > 800
    RT = 287.6;
else
    RT = 287.6;
end
end

function d = delta_cooling(T_g)
d = 0.02 + (T_g - 1200) / 100 * 0.02;
d = max(0.0, min(d, 0.15));
end

function [y, aux] = engine_forward(theta, cond)
eta_k     = theta(1);
eta_t     = theta(2);
eta_m     = theta(3);
eta_v     = theta(4);
eta_tv    = theta(5);
eta_c1    = theta(6);
eta_c2    = theta(7);
sigma_cc  = theta(8);
sigma_kan = theta(9);
sigma_kask= theta(10);
sigma_ks  = theta(11);
eta_T     = theta(12);
lambda    = theta(13);

T_H      = cond.T_H;
M_flight = cond.M_flight;
m        = cond.m;
pi_k     = cond.pi_k;
T_g      = cond.T_g;

y   = [NaN, NaN];
aux = struct();

try
    % (1) 基本常数
    k_air = 1.4;
    R_air = 287.3;
    a = sqrt(k_air * R_air * T_H);
    if ~isfinite(a) || a <= 0, return; end
    V_flight = a * M_flight;

    % (2) 分段函数值
    kT = piecewise_kT(T_g);
    RT = piecewise_RT(T_g);
    d  = delta_cooling(T_g);

    % (3) 进口总压比与压气机入口温度
    inner = 1 + V_flight^2 / (2 * (k_air/(k_air-1)) * R_air * T_H);
    if inner <= 0, return; end
    tau_v = inner^(k_air / (k_air - 1));
    T_B   = T_H * (inner^k_air);
    if ~isfinite(T_B) || T_B <= 0, return; end

    % (4) 压气机出口温度
    pi_k_ratio = pi_k^((k_air-1)/k_air);
    if ~isfinite(pi_k_ratio) || pi_k_ratio < 1, return; end
    T_k = T_B * (1 + (pi_k_ratio - 1) / eta_k);
    if ~isfinite(T_k) || T_k <= 0, return; end

    % (5) 相对耗油量
    g_T = 3e-5 * T_g - 2.69e-5 * T_k - 0.003;
    if ~isfinite(g_T) || g_T <= 0, return; end

    % (6) 热恢复系数
    compress_work = (k_air/(k_air-1)) * R_air * T_B * (pi_k_ratio - 1);
    gas_enthalpy  = (kT/(kT-1)) * RT * T_g;
    if abs(gas_enthalpy) < 1e-6, return; end
    num_lambda  = 1 - compress_work / (gas_enthalpy * eta_k);
    den_lambda  = 1 - compress_work / (gas_enthalpy * eta_k * eta_t);
    if abs(den_lambda) < 1e-10, return; end
    lambda_heat = num_lambda / den_lambda;
    if ~isfinite(lambda_heat), return; end

    % (7) 进口总压恢复系数
    sigma_bx = sigma_cc * sigma_kan;

    % (8) 单位自由能
    exp_T = (kT - 1) / kT;
    expansion_pr_denom = tau_v * sigma_bx * pi_k * sigma_kask * sigma_ks;
    if expansion_pr_denom <= 0, return; end
    expansion_term = (1.0 / expansion_pr_denom)^exp_T;
    if ~isfinite(expansion_term), return; end
    term1 = (kT / (kT - 1)) * RT * T_g * (1 - expansion_term);
    compress_work2 = (k_air / (k_air - 1)) * R_air * T_B * (pi_k^((k_air-1)/k_air) - 1);
    denom2 = (1 + g_T) * eta_k * eta_T * eta_t * eta_m * (1 - d);
    if abs(denom2) < 1e-10, return; end
    term2 = compress_work2 / denom2;
    L_sv  = lambda_heat * (term1 - term2);
    if ~isfinite(L_sv) || L_sv <= 0, return; end

    % (9) 最优自由能分配系数
    V2_term = m * V_flight^2;
    num_xpc = 1 + V2_term / (2 * L_sv * eta_tv * eta_v * eta_c2);
    den_xpc = 1 + (m * eta_tv * eta_v * eta_c2) / (eta_c1 * lambda);
    if abs(den_xpc) < 1e-10, return; end
    x_pc = num_xpc / den_xpc;
    if ~isfinite(x_pc) || x_pc <= 0, return; end

    % (10) 比推力
    inner_sq1 = 2 * eta_c1 * lambda * x_pc * L_sv;
    if inner_sq1 < 0, return; end
    V_j1 = (1 + g_T) * sqrt(inner_sq1) - V_flight;
    inner_sq2 = 2 * (1 - x_pc) / m * L_sv * eta_tv * eta_v * eta_c2 + V_flight^2;
    if inner_sq2 < 0, return; end
    V_j2  = sqrt(inner_sq2) - V_flight;
    R_ud  = (1/(1+m)) * V_j1 + (m/(1+m)) * V_j2;
    if ~isfinite(R_ud) || R_ud <= 0, return; end

    % (11) 比油耗
    denom_C = R_ud * (1 + m);
    if abs(denom_C) < 1e-10, return; end
    C_ud = 3600 * g_T * (1 - d) / denom_C;
    if ~isfinite(C_ud) || C_ud <= 0, return; end

    y = [R_ud, C_ud];
    aux.T_B        = T_B;
    aux.T_k        = T_k;
    aux.tau_v      = tau_v;
    aux.g_T        = g_T;
    aux.lambda_heat= lambda_heat;
    aux.sigma_bx   = sigma_bx;
    aux.L_sv       = L_sv;
    aux.x_pc       = x_pc;
    aux.kT         = kT;
    aux.RT         = RT;
    aux.delta      = d;

catch ME
    warning('engine_forward caught: %s', ME.message);
end
end
