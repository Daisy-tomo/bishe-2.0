%% parameter_sensitivity_analysis.m
% 涡扇发动机参数敏感性分析脚本
%
% 功能：
%   1. 基于 xiuzhengcanshu.m 中的前向热力学模型
%   2. 采用 OAT（One-At-a-Time）方法逐个扰动参数
%   3. 计算全范围扫描信噪比 (SNR) 和局部归一化敏感性指数 (SI)
%   4. 与 1% 噪声水平比较，标记强/中等/弱敏感参数
%   5. 生成可视化图形并保存结果
%
% 敏感性分级标准（基于信噪比 SNR = 全范围输出变化量 / 噪声标准差）：
%   强敏感：  SNR >= 5   （输出变化超过 5 倍噪声，参数可识别性强）
%   中等敏感：1 <= SNR < 5
%   弱敏感：  SNR <  1   （输出变化淹没在噪声中，参数难以辨识）
%
% 输出文件：
%   sensitivity_SNR.png      — SNR 条形图（全范围）
%   sensitivity_SI.png       — 归一化弹性系数图（局部）
%   sensitivity_sweep.png    — 各参数全范围扫描曲线
%   sensitivity_heatmap.png  — 敏感性热图
%   sensitivity_results.mat  — 数值结果（结构体）

clear; clc; close all;
rng(42);

fprintf('================================================\n');
fprintf('  涡扇发动机参数敏感性分析\n');
fprintf('================================================\n\n');

%% ============================================================
%  1. 参数定义（与 xiuzhengcanshu.m 保持一致）
%% ============================================================
param_names = {
    'eta_k',      'eta_t',     'eta_m',      'eta_v',    'eta_tv', ...
    'eta_c1',     'eta_c2',    'sigma_cc',   'sigma_kan','sigma_kask', ...
    'sigma_ks',   'eta_T',     'lambda'
};

param_desc = {
    '压气机绝热效率',       '涡轮绝热效率',         '机械效率', ...
    '风扇效率',             '风扇涡轮效率',         '一次喷管效率', ...
    '二次喷管效率',         '进气道激波总压恢复系数','进气通道总压恢复系数', ...
    '压气机级间总压恢复系数','燃烧室总压恢复系数',  '燃烧放热系数', ...
    '风扇涡轮热恢复系数'
};

% 参数下界与上界
lb = [0.84, 0.86, 0.980, 0.85, 0.90, 0.94, 0.92, 0.98, 0.98, 0.98, 0.94, 0.97, 1.02];
ub = [0.86, 0.92, 0.995, 0.87, 0.92, 0.95, 0.94, 1.00, 0.99, 0.99, 0.96, 0.99, 1.04];
n_params = length(lb);

% 真值参数（分析基准点）
theta_true = [
    0.85, 0.89, 0.988, 0.860, 0.910, 0.945, 0.930, ...
    0.990, 0.985, 0.985, 0.950, 0.980, 1.030
];

% 工况条件
cond.T_H      = 288;    % 环境温度 [K]
cond.M_flight = 0.0;    % 飞行 Mach 数（地面静止）
cond.m        = 10.0;   % 涵道比
cond.pi_k     = 33.0;   % 总压比
cond.T_g      = 1700.0; % 涡轮前燃气温度 [K]

% 噪声水平（与虚拟数据生成保持一致）
noise_level = 0.01;  % 1%

%% ============================================================
%  2. 基准输出计算
%% ============================================================
[y_base, aux_base] = engine_forward(theta_true, cond);
if ~all(isfinite(y_base))
    error('基准前向模型计算失败，请检查 theta_true 和 cond 设置');
end

R_base  = y_base(1);   % 比推力基准值 [N·s/kg]
C_base  = y_base(2);   % 比油耗基准值 [kg/(N·h)]
sigma_R = noise_level * abs(R_base);  % R_ud 噪声标准差
sigma_C = noise_level * abs(C_base);  % C_ud 噪声标准差

fprintf('工况条件：\n');
fprintf('  T_H = %.0f K,  M_flight = %.1f,  m = %.0f,  pi_k = %.0f,  T_g = %.0f K\n\n', ...
    cond.T_H, cond.M_flight, cond.m, cond.pi_k, cond.T_g);

fprintf('基准输出（theta_true）：\n');
fprintf('  R_ud  = %10.4f  N·s/kg\n',    R_base);
fprintf('  C_ud  = %12.6f  kg/(N·h)\n',  C_base);
fprintf('  sigma_R = %.4f  (噪声标准差)\n', sigma_R);
fprintf('  sigma_C = %.6f  (噪声标准差)\n', sigma_C);
fprintf('\n');

%% ============================================================
%  3. OAT 敏感性分析
%     3a. 全范围扫描：将参数从 lb 扫到 ub，记录输出变化量
%     3b. 局部差分：  双侧有限差分，估计局部导数及弹性系数
%% ============================================================
n_sweep   = 60;      % 每个参数的扫描点数（越大越精细）
delta_frac = 0.10;   % 局部差分扰动量 = delta_frac × 参数范围

% 结果存储
SNR_R_range = zeros(n_params, 1);   % 全范围扫描：R_ud 最大变化 / sigma_R
SNR_C_range = zeros(n_params, 1);   % 全范围扫描：C_ud 最大变化 / sigma_C
SI_R        = zeros(n_params, 1);   % 局部弹性系数（%/% 对 R_ud）
SI_C        = zeros(n_params, 1);   % 局部弹性系数（%/% 对 C_ud）
dR_abs      = zeros(n_params, 1);   % 全范围 R_ud 峰峰变化量
dC_abs      = zeros(n_params, 1);   % 全范围 C_ud 峰峰变化量
sweep_data  = cell(n_params, 1);    % 扫描原始数据（用于绘图）

fprintf('正在执行全范围参数扫描（%d 点/参数）...\n', n_sweep);

for i = 1:n_params
    range_i    = ub(i) - lb(i);
    theta_vals = linspace(lb(i), ub(i), n_sweep);
    R_vals     = nan(n_sweep, 1);
    C_vals     = nan(n_sweep, 1);

    % --- 全范围扫描 ---
    for j = 1:n_sweep
        th = theta_true;
        th(i) = theta_vals(j);
        [y_j, ~] = engine_forward(th, cond);
        if all(isfinite(y_j))
            R_vals(j) = y_j(1);
            C_vals(j) = y_j(2);
        end
    end

    sweep_data{i}.theta_vals = theta_vals;
    sweep_data{i}.R_vals     = R_vals;
    sweep_data{i}.C_vals     = C_vals;

    % 全范围峰峰变化
    dR_abs(i) = max(R_vals, [], 'omitnan') - min(R_vals, [], 'omitnan');
    dC_abs(i) = max(C_vals, [], 'omitnan') - min(C_vals, [], 'omitnan');
    SNR_R_range(i) = dR_abs(i) / sigma_R;
    SNR_C_range(i) = dC_abs(i) / sigma_C;

    % --- 局部双侧差分（弹性系数） ---
    d = delta_frac * range_i;
    th_p = theta_true; th_p(i) = min(theta_true(i) + d/2, ub(i));
    th_m = theta_true; th_m(i) = max(theta_true(i) - d/2, lb(i));
    [yp, ~] = engine_forward(th_p, cond);
    [ym, ~] = engine_forward(th_m, cond);
    act_d = th_p(i) - th_m(i);

    if all(isfinite(yp)) && all(isfinite(ym)) && act_d > 1e-12
        SI_R(i) = ((yp(1) - ym(1)) / R_base) / (act_d / theta_true(i)) * 100;
        SI_C(i) = ((yp(2) - ym(2)) / C_base) / (act_d / theta_true(i)) * 100;
    end

    fprintf('  [%2d/13] %-15s  SNR_R=%7.3f  SNR_C=%7.3f  SI_R=%+7.3f%%/%%  SI_C=%+7.3f%%/%%\n', ...
        i, param_names{i}, SNR_R_range(i), SNR_C_range(i), SI_R(i), SI_C(i));
end

fprintf('\n');

%% ============================================================
%  4. 敏感性分级
%% ============================================================
threshold_strong = 5.0;   % SNR 强敏感阈值
threshold_weak   = 1.0;   % SNR 弱敏感阈值

% 综合两个输出的最大 SNR 作为分级依据
SNR_max = max(SNR_R_range, SNR_C_range);

label       = cell(n_params, 1);
label_color = zeros(n_params, 3);
label_id    = zeros(n_params, 1);  % 1=强, 2=中, 3=弱（用于热图颜色映射）

for i = 1:n_params
    if SNR_max(i) >= threshold_strong
        label{i}       = '强敏感';
        label_color(i,:) = [0.82, 0.18, 0.18];
        label_id(i)    = 3;
    elseif SNR_max(i) < threshold_weak
        label{i}       = '弱敏感';
        label_color(i,:) = [0.20, 0.60, 0.25];
        label_id(i)    = 1;
    else
        label{i}       = '中等敏感';
        label_color(i,:) = [0.90, 0.65, 0.10];
        label_id(i)    = 2;
    end
end

strong_idx = find(SNR_max >= threshold_strong);
medium_idx = find(SNR_max >= threshold_weak & SNR_max < threshold_strong);
weak_idx   = find(SNR_max <  threshold_weak);

%% ============================================================
%  5. 打印汇总结果
%% ============================================================
fprintf('============================================================\n');
fprintf('  参数敏感性分析汇总结果\n');
fprintf('  （SNR = 全范围输出峰峰变化 / 噪声标准差，强>%.0f，弱<%.0f）\n', ...
    threshold_strong, threshold_weak);
fprintf('============================================================\n');
fprintf('%-15s  %-20s  %8s  %8s  %10s  %9s  %9s  %10s\n', ...
    '参数', '物理含义', 'SNR_R', 'SNR_C', 'SNR_max', 'SI_R(%%/%%)', 'SI_C(%%/%%)', '敏感性等级');
fprintf('%s\n', repmat('-', 1, 100));
for i = 1:n_params
    fprintf('%-15s  %-20s  %8.3f  %8.3f  %10.3f  %9.3f  %9.3f  %10s\n', ...
        param_names{i}, param_desc{i}, ...
        SNR_R_range(i), SNR_C_range(i), SNR_max(i), ...
        SI_R(i), SI_C(i), label{i});
end
fprintf('\n');

fprintf('强敏感参数  (%d 个): %s\n', length(strong_idx), strjoin(param_names(strong_idx), ', '));
fprintf('中等敏感参数 (%d 个): %s\n', length(medium_idx), strjoin(param_names(medium_idx), ', '));
fprintf('弱敏感参数  (%d 个): %s\n\n', length(weak_idx),   strjoin(param_names(weak_idx),   ', '));

%% ============================================================
%  6. 可视化
%% ============================================================

% ---- 图 1：SNR 条形图（全范围扫描）----
fig1 = figure('Name','参数敏感性 SNR 对比', 'Position',[50 50 1300 520]);

subplot(1,2,1);
b1 = bar(SNR_R_range, 'FaceColor','flat');
for i = 1:n_params; b1.CData(i,:) = label_color(i,:); end
hold on;
yl1 = yline(threshold_strong,'r--','LineWidth',2,'Label','强敏感阈值 (5)','LabelHorizontalAlignment','left');
yl2 = yline(threshold_weak,  'b--','LineWidth',2,'Label','弱敏感阈值 (1)','LabelHorizontalAlignment','left');
set(gca,'XTick',1:n_params,'XTickLabel',param_names,'XTickLabelRotation',45,'FontSize',9);
ylabel('SNR  =  \DeltaR_{ud} / \sigma_R');
title('比推力 R_{ud} 的参数敏感性（全范围 SNR）');
grid on; box on;
% 标注标签
for i = 1:n_params
    text(i, SNR_R_range(i)+0.08*max(SNR_R_range), label{i}(1), ...
        'HorizontalAlignment','center','FontSize',7,'Color',label_color(i,:));
end

subplot(1,2,2);
b2 = bar(SNR_C_range, 'FaceColor','flat');
for i = 1:n_params; b2.CData(i,:) = label_color(i,:); end
hold on;
yline(threshold_strong,'r--','LineWidth',2,'Label','强敏感阈值 (5)','LabelHorizontalAlignment','left');
yline(threshold_weak,  'b--','LineWidth',2,'Label','弱敏感阈值 (1)','LabelHorizontalAlignment','left');
set(gca,'XTick',1:n_params,'XTickLabel',param_names,'XTickLabelRotation',45,'FontSize',9);
ylabel('SNR  =  \DeltaC_{ud} / \sigma_C');
title('比油耗 C_{ud} 的参数敏感性（全范围 SNR）');
grid on; box on;
for i = 1:n_params
    text(i, SNR_C_range(i)+0.08*max(SNR_C_range), label{i}(1), ...
        'HorizontalAlignment','center','FontSize',7,'Color',label_color(i,:));
end

% 图例（统一在右下角）
legend_labels = {'强敏感 (红)', '中等敏感 (黄)', '弱敏感 (绿)'};
legend_colors = [0.82,0.18,0.18; 0.90,0.65,0.10; 0.20,0.60,0.25];
for k = 1:3
    patch('XData',[],'YData',[],'FaceColor',legend_colors(k,:),'DisplayName',legend_labels{k});
end
legend('show','Location','northeast');

sgtitle('OAT 全范围扫描参数敏感性分析（信噪比 SNR）','FontSize',13,'FontWeight','bold');

% ---- 图 2：局部归一化弹性系数 ----
fig2 = figure('Name','局部归一化弹性系数', 'Position',[100 120 1100 480]);

x = 1:n_params;
bw = 0.35;
b3 = bar(x - bw/2, SI_R, bw, 'FaceColor',[0.25 0.50 0.87], 'DisplayName','R_{ud}');
hold on;
b4 = bar(x + bw/2, SI_C, bw, 'FaceColor',[0.87 0.35 0.25], 'DisplayName','C_{ud}');
yline(0,'k-','LineWidth',1.0);
set(gca,'XTick',1:n_params,'XTickLabel',param_names,'XTickLabelRotation',45,'FontSize',9);
ylabel('弹性系数 SI (%/%)  — (Δy/y₀) / (Δθ/θ₀)');
legend('Location','northeast','FontSize',9);
title('局部双侧差分归一化敏感性指数（弹性系数）','FontSize',12);
grid on; box on;

% 颜色背景带（正/负区域）
ax = gca;
yl = ylim;
fill([0.5 n_params+0.5 n_params+0.5 0.5], [0 0 yl(2) yl(2)], ...
    [0.95 0.98 0.95],'EdgeColor','none','FaceAlpha',0.3);
fill([0.5 n_params+0.5 n_params+0.5 0.5], [yl(1) yl(1) 0 0], ...
    [0.98 0.95 0.95],'EdgeColor','none','FaceAlpha',0.3);
uistack(b3,'top'); uistack(b4,'top');

% ---- 图 3：全范围扫描曲线（每个参数一个子图）----
n_cols_plt = 4;
n_rows_plt = ceil(n_params / n_cols_plt);
fig3 = figure('Name','全范围扫描输出曲线', 'Position',[80 30 1500 900]);

for i = 1:n_params
    subplot(n_rows_plt, n_cols_plt, i);
    tv  = sweep_data{i}.theta_vals;
    Rv  = sweep_data{i}.R_vals;
    Cv  = sweep_data{i}.C_vals;

    yyaxis left;
    plot(tv, Rv, '-','Color',[0.15 0.45 0.85],'LineWidth',1.8);
    hold on;
    % ±1σ 噪声带（蓝色虚线）
    yline(R_base + sigma_R,'--','Color',[0.50 0.70 0.95],'LineWidth',0.9);
    yline(R_base - sigma_R,'--','Color',[0.50 0.70 0.95],'LineWidth',0.9);
    yline(R_base,':','Color',[0.15 0.45 0.85],'LineWidth',0.7);
    ylabel('R_{ud}','Color',[0.15 0.45 0.85],'FontSize',8);
    ax_l = gca; ax_l.YColor = [0.15 0.45 0.85];

    yyaxis right;
    plot(tv, Cv,'-','Color',[0.85 0.25 0.15],'LineWidth',1.8);
    yline(C_base + sigma_C,'--','Color',[0.95 0.65 0.60],'LineWidth',0.9);
    yline(C_base - sigma_C,'--','Color',[0.95 0.65 0.60],'LineWidth',0.9);
    yline(C_base,':','Color',[0.85 0.25 0.15],'LineWidth',0.7);
    ylabel('C_{ud}','Color',[0.85 0.25 0.15],'FontSize',8);
    ax_r = gca; ax_r.YColor = [0.85 0.25 0.15];

    % 真值竖线
    xline(theta_true(i),'k--','LineWidth',1.2);
    xlabel(param_names{i},'FontSize',8);
    title(sprintf('%s  [%s]', param_names{i}, label{i}), ...
        'Color', label_color(i,:), 'FontSize', 9, 'FontWeight','bold');
    grid on;
end

sgtitle(sprintf('各参数全范围扫描对输出的影响（虚线：基准值 ± 1\\sigma 噪声带，竖线：真值）'), ...
    'FontSize', 11, 'FontWeight','bold');

% ---- 图 4：敏感性热图 ----
fig4 = figure('Name','敏感性热图', 'Position',[150 150 850 320]);

hmap = [SNR_R_range'; SNR_C_range'];
imagesc(hmap);
colormap(parula);
cb = colorbar;
cb.Label.String = 'SNR 值';
set(gca,'XTick',1:n_params,'XTickLabel',param_names,'XTickLabelRotation',45, ...
    'YTick',1:2,'YTickLabel',{'R_{ud}','C_{ud}'},'FontSize',9);
title('参数敏感性热图（SNR：全范围输出变化 / 噪声标准差）','FontSize',11);
xlabel('参数');

% 数值标注
for i = 1:n_params
    for j = 1:2
        val = hmap(j,i);
        clr = 'white';
        if val < max(hmap(:))*0.4; clr = 'black'; end
        text(i, j, sprintf('%.2f\n%s', val, label{i}(1:2)), ...
            'HorizontalAlignment','center','VerticalAlignment','middle', ...
            'Color',clr,'FontSize',7,'FontWeight','bold');
    end
end

%% ============================================================
%  7. 保存图形与数值结果
%% ============================================================
figs   = {fig1, fig2, fig3, fig4};
fnames = {'sensitivity_SNR.png','sensitivity_SI.png', ...
          'sensitivity_sweep.png','sensitivity_heatmap.png'};

fprintf('保存图形...\n');
for k = 1:4
    try
        saveas(figs{k}, fnames{k});
        fprintf('  ✓ %s\n', fnames{k});
    catch ME
        fprintf('  ✗ %s 保存失败: %s\n', fnames{k}, ME.message);
    end
end

% 保存数值结果
res.param_names      = param_names;
res.param_desc       = param_desc;
res.theta_true       = theta_true;
res.lb               = lb;
res.ub               = ub;
res.R_base           = R_base;
res.C_base           = C_base;
res.sigma_R          = sigma_R;
res.sigma_C          = sigma_C;
res.noise_level      = noise_level;
res.SNR_R_range      = SNR_R_range;
res.SNR_C_range      = SNR_C_range;
res.SNR_max          = SNR_max;
res.SI_R             = SI_R;
res.SI_C             = SI_C;
res.label            = label;
res.threshold_strong = threshold_strong;
res.threshold_weak   = threshold_weak;
res.strong_params    = param_names(strong_idx);
res.medium_params    = param_names(medium_idx);
res.weak_params      = param_names(weak_idx);

save('sensitivity_results.mat', '-struct', 'res');
fprintf('  ✓ sensitivity_results.mat\n\n');
fprintf('分析完成。\n');

%% ============================================================
%  辅助函数（与 xiuzhengcanshu.m 相同，独立复制以保证脚本可单独运行）
%% ============================================================

% 燃气绝热指数（分段函数）
function kT = piecewise_kT(T_g)
    if     T_g > 800  && T_g <= 1400
        kT = 1.33;
    elseif T_g > 1400 && T_g <= 1600
        kT = 1.30;
    elseif T_g > 1600
        kT = 1.25;
    else
        kT = 1.33;
    end
end

% 燃气气体常数（分段函数）
function RT = piecewise_RT(T_g)
    if     T_g > 800  && T_g <= 1400
        RT = 287.6;
    elseif T_g > 1400 && T_g <= 1600
        RT = 288.0;
    elseif T_g > 1600
        RT = 288.6;
    else
        RT = 287.6;
    end
end

% 涡轮冷却引气系数
function d = delta_cooling(T_g)
    d = 0.02 + (T_g - 1200) / 100 * 0.02;
    d = max(0.0, min(d, 0.15));
end

% 前向热力学模型（涡扇发动机）
function [y, aux] = engine_forward(theta, cond)
    eta_k    = theta(1);   eta_t    = theta(2);   eta_m  = theta(3);
    eta_v    = theta(4);   eta_tv   = theta(5);   eta_c1 = theta(6);
    eta_c2   = theta(7);   sigma_cc = theta(8);   sigma_kan  = theta(9);
    sigma_kask = theta(10); sigma_ks = theta(11);  eta_T  = theta(12);
    lambda   = theta(13);

    T_H      = cond.T_H;
    M_flight = cond.M_flight;
    m        = cond.m;
    pi_k     = cond.pi_k;
    T_g      = cond.T_g;

    y   = [NaN, NaN];
    aux = struct();

    try
        k_air = 1.4; R_air = 287.3;
        a = sqrt(k_air * R_air * T_H);
        if ~isfinite(a) || a <= 0; return; end
        V_flight = a * M_flight;

        kT = piecewise_kT(T_g);
        RT = piecewise_RT(T_g);
        d  = delta_cooling(T_g);

        inner = 1 + V_flight^2 / (2 * (k_air/(k_air-1)) * R_air * T_H);
        if inner <= 0; return; end
        tau_v = inner^(k_air / (k_air - 1));
        T_B   = T_H * (inner^k_air);
        if ~isfinite(T_B) || T_B <= 0; return; end

        pi_k_ratio = pi_k^((k_air-1)/k_air);
        if ~isfinite(pi_k_ratio) || pi_k_ratio < 1; return; end
        T_k = T_B * (1 + (pi_k_ratio - 1) / eta_k);
        if ~isfinite(T_k) || T_k <= 0; return; end

        g_T = 3e-5 * T_g - 2.69e-5 * T_k - 0.003;
        if ~isfinite(g_T) || g_T <= 0; return; end

        compress_work = (k_air/(k_air-1)) * R_air * T_B * (pi_k_ratio - 1);
        gas_enthalpy  = (kT/(kT-1)) * RT * T_g;
        if abs(gas_enthalpy) < 1e-6; return; end

        num_lambda = 1 - compress_work / (gas_enthalpy * eta_k);
        den_lambda = 1 - compress_work / (gas_enthalpy * eta_k * eta_t);
        if abs(den_lambda) < 1e-10; return; end
        lambda_heat = num_lambda / den_lambda;
        if ~isfinite(lambda_heat); return; end

        sigma_bx = sigma_cc * sigma_kan;

        exp_T = (kT - 1) / kT;
        expansion_pr_denom = tau_v * sigma_bx * pi_k * sigma_kask * sigma_ks;
        if expansion_pr_denom <= 0; return; end
        expansion_term = (1.0 / expansion_pr_denom)^exp_T;
        if ~isfinite(expansion_term); return; end

        term1 = (kT / (kT - 1)) * RT * T_g * (1 - expansion_term);
        compress_work2 = (k_air / (k_air - 1)) * R_air * T_B * (pi_k^((k_air-1)/k_air) - 1);
        denom2 = (1 + g_T) * eta_k * eta_T * eta_t * eta_m * (1 - d);
        if abs(denom2) < 1e-10; return; end
        term2 = compress_work2 / denom2;
        L_sv  = lambda_heat * (term1 - term2);
        if ~isfinite(L_sv) || L_sv <= 0; return; end

        V2_term = m * V_flight^2;
        num_xpc = 1 + V2_term / (2 * L_sv * eta_tv * eta_v * eta_c2);
        den_xpc = 1 + (m * eta_tv * eta_v * eta_c2) / (eta_c1 * lambda);
        if abs(den_xpc) < 1e-10; return; end
        x_pc = num_xpc / den_xpc;
        if ~isfinite(x_pc) || x_pc <= 0; return; end

        inner_sq1 = 2 * eta_c1 * lambda * x_pc * L_sv;
        if inner_sq1 < 0; return; end
        V_j1 = (1 + g_T) * sqrt(inner_sq1) - V_flight;

        inner_sq2 = 2 * (1 - x_pc) / m * L_sv * eta_tv * eta_v * eta_c2 + V_flight^2;
        if inner_sq2 < 0; return; end
        V_j2 = sqrt(inner_sq2) - V_flight;

        R_ud = (1/(1+m)) * V_j1 + (m/(1+m)) * V_j2;
        if ~isfinite(R_ud) || R_ud <= 0; return; end

        denom_C = R_ud * (1 + m);
        if abs(denom_C) < 1e-10; return; end
        C_ud = 3600 * g_T * (1 - d) / denom_C;
        if ~isfinite(C_ud) || C_ud <= 0; return; end

        y = [R_ud, C_ud];

        aux.T_B       = T_B;   aux.T_k     = T_k;
        aux.tau_v     = tau_v; aux.g_T     = g_T;
        aux.lambda_heat = lambda_heat;
        aux.sigma_bx  = sigma_bx;
        aux.L_sv      = L_sv;  aux.x_pc    = x_pc;
        aux.kT        = kT;    aux.RT      = RT;
        aux.delta     = d;

    catch ME
        warning('engine_forward caught: %s', ME.message);
    end
end
