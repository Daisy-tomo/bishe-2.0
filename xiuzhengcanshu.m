
clear; clc; close all;
rng(42);  
 

param_names = {
    'eta_k',       ...  % 压气机绝热效率
    'eta_t',       ...  % 涡轮绝热效率
    'eta_m',       ...  % 机械效率
    'eta_v',       ...  % 风扇效率
    'eta_tv',      ...  % 风扇涡轮效率
    'eta_c1',      ...  % 一次喷管效率
    'eta_c2',      ...  % 二次喷管效率
    'sigma_cc',    ...  % 进气道激波总压恢复系数
    'sigma_kan',   ...  % 进气通道总压恢复系数
    'sigma_kask',  ...  % 压气机级间总压恢复系数
    'sigma_ks',    ...  % 燃烧室总压恢复系数
    'eta_T',       ...  % 燃烧放热系数
    'lambda'        ...  % 风扇涡轮热恢复系数
};

lb = [0.84, 0.86, 0.980, 0.85, 0.90, 0.94, 0.92, 0.98, 0.98, 0.98, 0.94, 0.97, 1.02];
ub = [0.86, 0.92, 0.995, 0.87, 0.92, 0.95, 0.94, 1.00, 0.99, 0.99, 0.96, 0.99, 1.04];

n_params = length(lb);


cond.T_H = 288;      
cond.M_flight = 0.0;     
cond.m = 10.0;    
cond.pi_k = 33.0;    
cond.T_g = 1700.0;  

fprintf('T_H = %.2f K\n', cond.T_H);
fprintf('M_flight = %.2f\n', cond.M_flight);
fprintf('m = %.1f\n', cond.m);
fprintf('pi_k = %.1f\n', cond.pi_k);
fprintf('T_g = %.1f K\n', cond.T_g);
fprintf('\n');


theta_true = [
    0.85,    ...  % eta_k
    0.89,    ...  % eta_t
    0.988,   ...  % eta_m
    0.860,   ...  % eta_v
    0.910,   ...  % eta_tv
    0.945,   ...  % eta_c1
    0.930,   ...  % eta_c2
    0.990,   ...  % sigma_cc
    0.985,   ...  % sigma_kan
    0.985,   ...  % sigma_kask
    0.950,   ...  % sigma_ks
    0.980,   ...  % eta_T
    1.030    ...  % lambda
];

assert(all(theta_true >= lb) && all(theta_true <= ub), ...
    'theta_true 超出参数范围');

theta0 = 0.5 * (lb + ub);

fprintf('真值参数 theta_true：\n');
for i = 1:n_params
    fprintf('  %-15s = %.4f  [取值区间: %.4f ~ %.4f]\n', ...
        param_names{i}, theta_true(i), lb(i), ub(i));
end
fprintf('\n初始值 theta0：\n');
for i = 1:n_params
    fprintf('  %-15s = %.4f\n', param_names{i}, theta0(i));
end
fprintf('\n');


fprintf('前向模型验证：\n');
[y_true, aux_true] = engine_forward(theta_true, cond);
if ~all(isfinite(y_true))
    error('前向模型在 theta_true 处输出非有限值');
end
fprintf('theta_true 前向计算结果：\n');
fprintf('R_ud_true = %.4f [N·s/kg]\n', y_true(1));
fprintf('C_ud_true = %.6f [kg/(N·h)]\n', y_true(2));
fprintf('中间变量：\n');
fprintf('T_B = %.2f K, T_k = %.2f K\n', aux_true.T_B, aux_true.T_k);
fprintf('g_T = %.6f\n', aux_true.g_T);
fprintf('lambda_heat = %.4f\n', aux_true.lambda_heat);
fprintf('L_sv = %.2f J/kg\n', aux_true.L_sv);
fprintf('x_pc = %.4f\n', aux_true.x_pc);
fprintf('tau_v = %.4f\n', aux_true.tau_v);
fprintf('\n');


noise_level_R = 0.01;  
noise_level_C = 0.01;   
rng_seed = 123;   

data = generate_virtual_data(theta_true, cond, noise_level_R, noise_level_C, rng_seed);

fprintf('试验数据：\n');
fprintf('R_true = %.4f, R_obs = %.4f, sigma_R = %.4f\n', ...
    data.R_true, data.R_obs, data.sigma_R);
fprintf('C_true = %.6f, C_obs = %.6f, sigma_C = %.6f\n', ...
    data.C_true, data.C_obs, data.sigma_C);
fprintf('\n');


opts.n_samples = 10000;  % 总采样步数
opts.burn_in = 2000;   % 烧入期步数
opts.proposal_sd = 0.012;  % 标准差
opts.adapt_start = 0;    % 开始
opts.adapt_end = 2000;   % 结束
opts.adapt_interval = 100;    % 自适应调整间隔
opts.target_accept_low = 0.18;   % 目标接受率下限
opts.target_accept_high = 0.35;   % 目标接受率上限

fprintf('MCMC 配置：\n');
fprintf('  总采样步数:%d\n', opts.n_samples);
fprintf('  烧入期:%d\n', opts.burn_in);
fprintf('  后验链长度:%d\n', opts.n_samples - opts.burn_in);
fprintf('  初始标准差:%.4f\n', opts.proposal_sd);
fprintf('\n');

%% 运行 MCMC
results = run_mcmc(data, cond, lb, ub, theta0, opts);
fprintf('后验结果：\n');
fprintf('总体接受率: %.3f\n\n', results.accept_rate);

fprintf('%-10s  %8s  %8s  %8s  %8s  %8s  %8s\n', ...
    '参数名', '真值', 'Post均值', 'MAP', 'CI95低', 'CI95高', '先验中点');
fprintf('%s\n', repmat('-', 1, 75));
for i = 1:n_params
    fprintf('%-15s  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f\n', ...
        param_names{i}, ...
        theta_true(i), ...
        results.theta_mean(i), ...
        results.theta_map(i), ...
        results.theta_ci95(i,1), ...
        results.theta_ci95(i,2), ...
        theta0(i));
end
fprintf('\n');

%% 后验预测
[y_mean, ~] = engine_forward(results.theta_mean, cond);
[y_map,  ~] = engine_forward(results.theta_map,  cond);

R_pred_mean = y_mean(1); 
C_pred_mean = y_mean(2);
R_pred_map= y_map(1);   
C_pred_map = y_map(2);

fprintf('后验预测对比：\n');
fprintf('  %-20s  %10s  %10s\n', '性能指标', 'R_ud [N·s/kg]', 'C_ud [kg/(N·h)]');
fprintf('  %-20s  %10.4f  %10.6f\n', '真值', data.R_true,data.C_true);
fprintf('  %-20s  %10.4f  %10.6f\n', '观测值', data.R_obs,data.C_obs);
fprintf('  %-20s  %10.4f  %10.6f\n', '后验均值预测',R_pred_mean,C_pred_mean);
fprintf('  %-20s  %10.4f  %10.6f\n', 'MAP预测',R_pred_map,C_pred_map);
fprintf('\n');


err_mean_R = abs(R_pred_mean - data.R_true) / abs(data.R_true) * 100;
err_mean_C = abs(C_pred_mean - data.C_true) / abs(data.C_true) * 100;
err_map_R  = abs(R_pred_map  - data.R_true) / abs(data.R_true) * 100;
err_map_C  = abs(C_pred_map  - data.C_true) / abs(data.C_true) * 100;
fprintf('后验均值预测相对误差: R_ud=%.3f%%, C_ud=%.3f%%\n', err_mean_R, err_mean_C);
fprintf('MAP预测相对误差: R_ud=%.3f%%, C_ud=%.3f%%\n', err_map_R,  err_map_C);
fprintf('\n');

plot_results(results, theta_true, lb, ub, param_names);

%  燃气绝热指数
function kT = piecewise_kT(T_g)
if T_g > 800 && T_g <= 1400
    kT = 1.33;
elseif T_g > 1400 && T_g <= 1600
    kT = 1.30;
elseif T_g > 1600
    kT = 1.25;
else
    kT = 1.33;
end
end

% 燃气气体常数 
function RT = piecewise_RT(T_g)
if T_g > 800 && T_g <= 1400
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


function [y, aux] = engine_forward(theta, cond)
eta_k= theta(1);   
eta_t= theta(2);  
eta_m= theta(3);   
eta_v= theta(4);   
eta_tv= theta(5);   
eta_c1= theta(6);   
eta_c2= theta(7);   
sigma_cc= theta(8);   
sigma_kan= theta(9);   
sigma_kask= theta(10);  
sigma_ks= theta(11);  
eta_T= theta(12);  
lambda= theta(13);  


T_H = cond.T_H;
M_flight = cond.M_flight;
m = cond.m;
pi_k = cond.pi_k;
T_g  = cond.T_g;


y = [NaN, NaN];
aux = struct();

try
 

    % (1) 基本常数
    k_air = 1.4;       % 空气绝热指数
    R_air = 287.3;     % 空气气体常数
    a = sqrt(k_air * R_air * T_H);
    if ~isfinite(a) || a <= 0, return; end
    V_flight = a * M_flight;  

   
    % (2) 分段函数值
    kT = piecewise_kT(T_g);    
    RT = piecewise_RT(T_g);    
    d = delta_cooling(T_g);  


    % (3) 进口总压比与压气机入口温度
    inner = 1 + V_flight^2 / (2 * (k_air/(k_air-1)) * R_air * T_H);
    if inner <= 0, return; end
    tau_v = inner^(k_air / (k_air - 1));
    T_B = T_H * (inner^k_air);  
    if ~isfinite(T_B) || T_B <= 0, return; end

    % (4) 压气机出口温度
    pi_k_ratio = pi_k^((k_air-1)/k_air);  
    if ~isfinite(pi_k_ratio) || pi_k_ratio < 1, return; end
    T_k = T_B * (1 + (pi_k_ratio - 1) / eta_k);
    if ~isfinite(T_k) || T_k <= 0, return; end

   
    % (5) 相对耗油量
    g_T = 3e-5 * T_g - 2.69e-5 * T_k - 0.003;
    if ~isfinite(g_T) || g_T <= 0
        return;
    end

    % (6) 热恢复系数 
    compress_work = (k_air/(k_air-1)) * R_air * T_B * (pi_k_ratio - 1);
    gas_enthalpy = (kT/(kT-1)) * RT * T_g;
    if abs(gas_enthalpy) < 1e-6, return; end

    num_lambda = 1 - compress_work / (gas_enthalpy * eta_k );
    den_lambda = 1 - compress_work / (gas_enthalpy * eta_k * eta_t);

    if abs(den_lambda) < 1e-10, return; end
    lambda_heat = num_lambda / den_lambda;
    if ~isfinite(lambda_heat), return; end

    % (7) 进口总压恢复系数 
    sigma_bx = sigma_cc * sigma_kan;


   %(8) 单位自由能 
exp_T = (kT - 1) / kT;
expansion_pr_denom = tau_v * sigma_bx * pi_k * sigma_kask * sigma_ks;
if expansion_pr_denom <= 0
    return;
end
expansion_term = (1.0 / expansion_pr_denom)^exp_T;
if ~isfinite(expansion_term)
    return;
end
term1 = (kT / (kT - 1)) * RT * T_g * (1 - expansion_term);
compress_work = (k_air / (k_air - 1)) * R_air * T_B * (pi_k^((k_air - 1) / k_air) - 1);
denom2 = (1 + g_T) * eta_k * eta_T * eta_t * eta_m * (1 - d);
if abs(denom2) < 1e-10
    return;
end
term2 = compress_work / denom2;
L_sv = lambda_heat * (term1 - term2);
if ~isfinite(L_sv) || L_sv <= 0
    return;
end

  % (9) 最优自由能分配系数 
    if L_sv <= 0, return; end
    V2_term = m * V_flight^2;   % = 0 when V_flight=0
    num_xpc = 1 + V2_term / (2 * L_sv * eta_tv * eta_v * eta_c2);
    den_xpc = 1 + (m * eta_tv * eta_v * eta_c2) / (eta_c1 * lambda);
    if abs(den_xpc) < 1e-10, return; end
    x_pc = num_xpc / den_xpc;
    if ~isfinite(x_pc) || x_pc <= 0 || x_pc >= 1
        if x_pc <= 0, return; end
    end

    % (10) 比推力 
    inner_sq1 = 2 * eta_c1 * lambda * x_pc * L_sv;
    if inner_sq1 < 0, return; end
    V_j1 = (1 + g_T) * sqrt(inner_sq1) - V_flight;
    inner_sq2 = 2 * (1 - x_pc) / m * L_sv * eta_tv * eta_v * eta_c2 + V_flight^2;
    if inner_sq2 < 0, return; end
    V_j2 = sqrt(inner_sq2) - V_flight;
    R_ud = (1/(1+m)) * V_j1 + (m/(1+m)) * V_j2;
    if ~isfinite(R_ud) || R_ud <= 0, return; end

    % (11) 比油耗 
    denom_C = R_ud * (1 + m);
    if abs(denom_C) < 1e-10, return; end
    C_ud = 3600 * g_T * (1 - d) / denom_C;
    if ~isfinite(C_ud) || C_ud <= 0, return; end


 
    y = [R_ud, C_ud];

    aux.T_B= T_B;
    aux.T_k= T_k;
    aux.tau_v= tau_v;
    aux.g_T= g_T;
    aux.lambda_heat= lambda_heat;
    aux.sigma_bx= sigma_bx;
    aux.L_sv= L_sv;
    aux.x_pc= x_pc;
    aux.kT= kT;
    aux.RT= RT;
    aux.delta= d;

catch ME
    warning('engine_forward caught: %s', ME.message);
end
end



%%
function data = generate_virtual_data(theta_true, cond, noise_level_R, noise_level_C, rng_seed)
rng(rng_seed);  

[y_true, ~] = engine_forward(theta_true, cond);
if ~all(isfinite(y_true))
    error('generate_virtual_data: 真值前向模型返回非有限值');
end
R_true = y_true(1);
C_true = y_true(2);

sigma_R = noise_level_R * abs(R_true);
sigma_C = noise_level_C * abs(C_true);

R_obs = R_true + sigma_R * randn();
C_obs = C_true + sigma_C * randn();

data.R_true = R_true;
data.C_true = C_true;
data.R_obs = R_obs;
data.C_obs = C_obs;
data.sigma_R = sigma_R;
data.sigma_C = sigma_C;
end


function lp = log_prior(theta, lb, ub) %对数先验
if all(theta >= lb) && all(theta <= ub)
    lp = 0.0;
else
    lp = -Inf;
end
end


function ll = log_likelihood(theta, data, cond) %对数似然
ll = -Inf;  

[y_pred, ~] = engine_forward(theta, cond);
if ~all(isfinite(y_pred))
    return;
end

R_pred = y_pred(1);
C_pred = y_pred(2);
R_obs = data.R_obs;
C_obs = data.C_obs;
sigma_R = data.sigma_R;
sigma_C = data.sigma_C;

if sigma_R <= 0 || sigma_C <= 0
    return;
end

res_R = (R_obs - R_pred) / sigma_R;
res_C = (C_obs - C_pred) / sigma_C;
ll = -0.5 * (res_R^2 + log(2*pi*sigma_R^2) + res_C^2 + log(2*pi*sigma_C^2));

if ~isfinite(ll)
    ll = -Inf;
end
end


function lpost = log_posterior_fn(theta, data, cond, lb, ub) %对数后验
lp = log_prior(theta, lb, ub);
if isinf(lp)
    lpost = -Inf;
    return;
end
ll = log_likelihood(theta, data, cond);
lpost = lp + ll;
end

%%  
function z_ref = reflect(z)
z_ref = z;
for i = 1:length(z)
    zi = z(i);

    max_iter = 10;
    iter = 0;
    while (zi < 0 || zi > 1) && iter < max_iter
        if zi < 0
            zi = -zi;         
        elseif zi > 1
            zi = 2 - zi;        
        end
        iter = iter + 1;
    end
    z_ref(i) = max(0, min(1, zi));
end
end



%% 
function results = run_mcmc(data, cond, lb, ub, theta0, opts)
n_params = length(lb);
n_samples = opts.n_samples;
burn_in = opts.burn_in;
prop_sd = opts.proposal_sd;

chain_full = zeros(n_samples, n_params);
logpost_full = zeros(n_samples, 1);

theta_curr = theta0;
theta_curr = max(lb, min(ub, theta_curr));
z_curr = (theta_curr - lb) ./ (ub - lb);
lpost_curr = log_posterior_fn(theta_curr, data, cond, lb, ub);

if ~isfinite(lpost_curr)
    warning('初始点对数后验无效，改用先验中点');
    theta_curr = 0.5*(lb+ub);
    z_curr = 0.5 * ones(1, n_params);
    lpost_curr = log_posterior_fn(theta_curr, data, cond, lb, ub);
    if ~isfinite(lpost_curr)
        error('先验中点也无效\n');
    end
end

n_accept_total = 0;
n_accept_window = 0;

for s = 1:n_samples
    z_prop = z_curr + prop_sd * randn(1, n_params);
    z_prop = reflect(z_prop);
    theta_prop = lb + z_prop .* (ub - lb);
    lpost_prop = log_posterior_fn(theta_prop, data, cond, lb, ub);

    
    log_alpha = lpost_prop - lpost_curr;
    if log(rand()) < log_alpha
        z_curr = z_prop;
        theta_curr = theta_prop;
        lpost_curr = lpost_prop;
        n_accept_total = n_accept_total  + 1;
        n_accept_window = n_accept_window + 1;
    end

    
    chain_full(s, :) = theta_curr;
    logpost_full(s)  = lpost_curr;

    
    if s >= opts.adapt_start && s <= opts.adapt_end
        if mod(s - opts.adapt_start, opts.adapt_interval) == 0 && s > opts.adapt_start
            local_rate = n_accept_window / opts.adapt_interval;
            if local_rate < opts.target_accept_low
                prop_sd = prop_sd * 0.9;
            elseif local_rate > opts.target_accept_high
                prop_sd = prop_sd * 1.1;
            end
            prop_sd = max(1e-4, min(prop_sd, 0.5));
            n_accept_window = 0;
        end
    end
end


chain_post = chain_full(burn_in+1:end, :);
logpost_post = logpost_full(burn_in+1:end);

theta_mean = mean(chain_post, 1);
theta_std = std(chain_post, 0, 1);

[~, best_idx_post] = max(logpost_post);
theta_map = chain_post(best_idx_post, :);

theta_ci95 = zeros(n_params, 2);
for i = 1:n_params
    theta_ci95(i, :) = quantile(chain_post(:,i), [0.025, 0.975]);
end

results.chain_full = chain_full;
results.chain_post = chain_post;
results.logpost_full = logpost_full;
results.logpost_post = logpost_post;
results.accept_rate = n_accept_total / n_samples;
results.theta_mean = theta_mean;
results.theta_std = theta_std;
results.theta_map = theta_map;
results.theta_ci95 = theta_ci95;
results.best_idx = best_idx_post;
results.prop_sd_final = prop_sd;
end

%% 
function plot_results(results, theta_true, lb, ub, param_names)

chain_post = results.chain_post;
chain_full = results.chain_full;
n_params   = size(chain_post, 2);
n_post     = size(chain_post, 1);
n_full     = size(chain_full, 1);



fig1 = figure('Name','链轨迹图','Position',[50 50 1400 900]);
n_cols = 3;
n_rows = ceil(n_params / n_cols);
for i = 1:n_params
    subplot(n_rows, n_cols, i);
    plot(1:n_full, chain_full(:,i), 'b-', 'LineWidth', 0.3);
    hold on;
    xline(size(chain_full,1) - n_post, 'r--', 'LineWidth', 1.5);  % 去掉 Label
    yline(theta_true(i), 'g-', 'LineWidth', 1.5);
    xlabel('迭代步数');
    ylabel(param_names{i});
    title(sprintf('%s 链轨迹', param_names{i}));
    ylim([lb(i)-0.005*(ub(i)-lb(i)), ub(i)+0.005*(ub(i)-lb(i))]);
end

ax_leg = axes('Position', [0.72 0.88 0.001 0.001], 'Visible', 'off');
hold(ax_leg, 'on');
plot(ax_leg, NaN, NaN, 'b-',  'LineWidth', 1.0);
plot(ax_leg, NaN, NaN, 'r--', 'LineWidth', 1.5);
plot(ax_leg, NaN, NaN, 'g-',  'LineWidth', 1.5);
lgd = legend(ax_leg, '链', 'burn-in边界', '真值', 'FontSize', 8);
lgd.Units    = 'normalized';
lgd.Position = [0.72 0.87 0.08 0.08];

sgtitle('MCMC 链轨迹图');


fig2 = figure('Name','边缘后验直方图','Position',[100 50 1400 900]);
for i = 1:n_params
    subplot(n_rows, n_cols, i);
    histogram(chain_post(:,i), 40, 'Normalization', 'pdf', ...
        'FaceColor', [0.4 0.6 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.7);
    hold on;
    xline(theta_true(i), 'g-',  'LineWidth', 2);
    xline(results.theta_mean(i), 'r--', 'LineWidth', 2);
    xline(results.theta_map(i),  'm:',  'LineWidth', 2);
    xlabel(param_names{i});
    ylabel('概率密度');
    title(sprintf('%s\n真值=%.4f, 均值=%.4f, MAP=%.4f', ...
        param_names{i}, theta_true(i), results.theta_mean(i), results.theta_map(i)));
    xlim([lb(i), ub(i)]);

    if i == 1
        legend('后验', '真值', '后验均值', 'MAP', 'Location', 'best', 'FontSize', 8);
    end
end
sgtitle('边缘后验分布直方图');


try
    saveas(fig1, '链轨迹图.png');
    saveas(fig2, '边缘后验图.png');
catch
    fprintf('图形保存失败\n');
end
end



