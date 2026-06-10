%% 主动声呐信号分析：CW与LFM脉冲复现
% 作者：不知周（由AI辅助完成）
% 日期：2026-06-03
% 中心频率：80kHz
clear; clc; close all;

%% ====================== 全局参数设置 ======================
f0 = 80e3;       % 中心频率 80kHz
T = 1e-3;        % 脉冲宽度 1ms（可修改）
Fs = 1e6;        % 采样频率 1MHz（满足奈奎斯特采样定理）
A = 1;           % 信号幅度
dt = 1/Fs;       % 采样间隔

%% ====================== 第一部分：CW单频矩形脉冲 ======================
fprintf('正在生成CW脉冲信号...\n');

% 1. 时间函数与波形（对应图2.4）
t_cw = 0:dt:T-dt;  % 时间轴 [0, T)
s_cw = A * exp(1j * 2 * pi * f0 * t_cw);  % 复信号形式

figure('Name','图2.4 CW脉冲信号波形','Position',[100,100,800,400]);
plot(t_cw*1e3, real(s_cw), 'LineWidth',1.2);
xlabel('时间 t (ms)','FontSize',12);
ylabel('Re[s(t)]','FontSize',12);
title('CW脉冲信号波形（实部）','FontSize',14);
grid on;
axis tight;
ylim([-1.1*A, 1.1*A]);

% 2. 频谱函数（对应图2.5）
N_cw = length(s_cw);
f_cw = (-Fs/2:Fs/N_cw:Fs/2-Fs/N_cw);  % 频率轴(Hz)
S_cw = fftshift(fft(s_cw)) * dt;      % 幅度修正

figure('Name','图2.5 CW脉冲信号频谱','Position',[100,100,800,400]);
amp = abs(S_cw);
amp_norm = amp / max(amp);   % 正确归一化
plot(f_cw, amp_norm, 'LineWidth',1.2);

xlabel('频率 f (Hz)','FontSize',12);
ylabel('归一化|S(f)|','FontSize',12);
title('CW脉冲信号频谱','FontSize',14);
grid on;

axis tight;                  % 第一步：自适应全图
xlim([75000,85000]);         % 第二步：截取83k~87kHz，不会被覆盖

% 3. 模糊函数计算（对应图2.6）
tau_range = linspace(-T, T, 201);  % 时延范围 [-T, T]
xi_range = linspace(-2/T, 2/T, 201);  % 频移范围 [-2/T, 2/T]
[TAU, XI] = meshgrid(tau_range, xi_range);
chi_cw = zeros(size(TAU));

for i = 1:length(xi_range)
    for j = 1:length(tau_range)
        tau = tau_range(j);
        xi = xi_range(i);
        if abs(tau) > T
            chi_cw(i,j) = 0;
        else
            arg = pi * xi * (T - abs(tau));
            if arg == 0  % 处理分母为0的极限情况
                chi_cw(i,j) = T - abs(tau);
            else
                chi_cw(i,j) = abs(sin(arg)/arg * (T - abs(tau)));
            end
        end
    end
end

% (a) 模糊函数三维图
figure('Name','图2.6(a) CW模糊函数三维图','Position',[100,100,800,600]);
surf(TAU*1e3, XI*T, chi_cw/T, 'EdgeColor','none');  % 幅度归一化到T
xlabel('时延 τ ','FontSize',12);
ylabel('频移 ξ ','FontSize',12);
zlabel('|χ(τ,ξ)|','FontSize',12);
title('CW脉冲信号模糊函数','FontSize',14);
shading interp;
colormap jet;
view(30, 30);  % 视角与教材一致
colorbar;

% (b) ξ=0截面（时延截面）
chi_cw_xi0 = T - abs(tau_range);
figure('Name','图2.6(b) CW模糊函数ξ=0截面','Position',[100,100,800,400]);
plot(tau_range/T, chi_cw_xi0/T, 'LineWidth',1.2);
xlabel('时延 τ / T','FontSize',12);
ylabel('|χ(τ,0)| / T','FontSize',12);
title('CW模糊函数 ξ=0 截面','FontSize',14);
grid on;
hold on;
% 标记0.707T对应的τ≈0.3T（教材标注）
tau_07 = 0.3*T;
plot([-tau_07*1e3, tau_07*1e3], [0.7, 0.7], 'r--', 'LineWidth',1.2);
plot([-tau_07*1e3, -tau_07*1e3], [0, 0.7], 'r--', 'LineWidth',1.2);
plot([tau_07*1e3, tau_07*1e3], [0, 0.7], 'r--', 'LineWidth',1.2);
text(-tau_07*1e3-0.08, 0.72, '0.7','Color','red','FontSize',10);
text(-tau_07*1e3, -0.08, '-0.3','Color','red','FontSize',10);
text(tau_07*1e3, -0.08, '0.3','Color','red','FontSize',10);
hold off;
axis tight;
ylim([0, 1.05]);

% (c) τ=0截面（频移截面）
chi_cw_tau0 = abs(sin(pi*xi_range*T)./(pi*xi_range*T)) * T;
chi_cw_tau0(xi_range==0) = T;  % 处理xi=0的情况
figure('Name','图2.6(c) CW模糊函数τ=0截面','Position',[100,100,800,400]);
plot(xi_range*T, chi_cw_tau0/T, 'LineWidth',1.2);
xlabel('频移 ξ * T','FontSize',12);
ylabel('|χ(0,ξ)| / T','FontSize',12);
title('CW模糊函数 τ=0 截面','FontSize',14);
grid on;
hold on;
% 标记0.707T对应的ξ≈0.44/T（教材标注）
xi_07 = 0.44/T;
plot([-xi_07*T, xi_07*T], [0.707, 0.707], 'r--', 'LineWidth',1.2);
plot([-xi_07*T, -xi_07*T], [0, 0.707], 'r--', 'LineWidth',1.2);
plot([xi_07*T, xi_07*T], [0, 0.707], 'r--', 'LineWidth',1.2);
text(-xi_07*T-0.1, 0.75, '0.7','Color','red','FontSize',10);
text(-xi_07*T+0.05, +0.08, '-0.44','Color','red','FontSize',10);
text(xi_07*T+0.05, +0.08, '0.44','Color','red','FontSize',10);
hold off;
axis tight;
ylim([0, 1.05]);

% (d) 模糊度图（等高线图，自动计算-3dB轮廓数值）
figure('Name','图2.6(d) 80kHz CW脉冲模糊度图','Position',[100,100,600,600]);
% 绘制-3dB（0.707归一化幅度）等高线，同时获取等高线数据矩阵C
[C, h] = contour(TAU/T, XI*T, chi_cw/T, [0.707, 0.707], 'LineWidth',1.2);
xlabel('时延 τ/T（无量纲）','FontSize',12);
ylabel('频移 ξ·T（无量纲）','FontSize',12);
title('80kHz CW脉冲信号模糊度图（-3dB轮廓）','FontSize',14);
grid on;
axis equal;  % 保持坐标轴比例一致，椭圆不变形
axis tight;

% ---------------- 核心：自动计算-3dB轮廓的数值范围 ----------------
% 解析MATLAB等高线矩阵C的格式：
% C(1,1) = 等高线数值；C(2,1) = 该等高线的点数；后续列依次为(x,y)坐标
contour_level = C(1,1);  % 提取等高线数值（应为0.707）
num_points = C(2,1);     % 提取该等高线的总点数
x_coords = C(1, 2:num_points+1);  % 所有点的x坐标（τ/T）
y_coords = C(2, 2:num_points+1);  % 所有点的y坐标（ξ·T）

% 计算-3dB轮廓的极值（即椭圆的长半轴和短半轴）
tau_min = min(x_coords);
tau_max = max(x_coords);
xi_min = min(y_coords);
xi_max = max(y_coords);

% 打印计算结果到命令行（直接显示数值，无需手动量图）
fprintf('\n=== 80kHz CW信号模糊度图-3dB轮廓自动计算结果 ===\n');
fprintf('等高线归一化幅度：%.3f（理论值0.707）\n', contour_level);
fprintf('时延范围 τ/T：%.3f ~ %.3f（教材理论值±0.3）\n', tau_min, tau_max);
fprintf('频移范围 ξ·T：%.3f ~ %.3f（教材理论值±0.44）\n', xi_min, xi_max);
fprintf('===============================================\n');

% 在图上标注计算出的数值（红色加粗，和教材标注位置一致）
text(tau_min-0.06, 0, sprintf('%.2f', tau_min), 'Color','red','FontSize',10,'FontWeight','bold');
text(tau_max+0.02, 0, sprintf('%.2f', tau_max), 'Color','red','FontSize',10,'FontWeight','bold');
text(0, xi_max+0.04, sprintf('%.2f', xi_max), 'Color','red','FontSize',10,'FontWeight','bold');
text(0, xi_min-0.06, sprintf('%.2f', xi_min), 'Color','red','FontSize',10,'FontWeight','bold');





% (d) 模糊度图（等高线图，全量纲自动计算-3dB轮廓数值）
figure('Name','图2.6(d) 80kHz CW脉冲模糊度图','Position',[100,100,600,600]);
% 绘制-3dB（0.707归一化幅度）等高线，同时获取等高线数据矩阵C
[C, h] = contour(TAU*1e3, XI, chi_cw/T, [0.707, 0.707], 'LineWidth',1.2);
xlabel('时延 τ (ms)','FontSize',12);  % 量纲：毫秒
ylabel('频移 ξ (Hz)','FontSize',12);  % 量纲：赫兹
title('80kHz CW脉冲信号模糊度图（-3dB轮廓）','FontSize',14);
grid on;
axis tight;

% ---------------- 核心：全量纲自动计算-3dB轮廓数值 ----------------
% 解析MATLAB等高线矩阵C
contour_level = C(1,1);  % 提取等高线数值（应为0.707）
num_points = C(2,1);     % 提取该等高线的总点数
x_coords = C(1, 2:num_points+1);  % 所有点的x坐标（τ，单位：ms）
y_coords = C(2, 2:num_points+1);  % 所有点的y坐标（ξ，单位：Hz）

% 计算-3dB轮廓的极值（椭圆的长半轴和短半轴）
tau_min_ms = min(x_coords);
tau_max_ms = max(x_coords);
xi_min_hz = min(y_coords);
xi_max_hz = max(y_coords);

% 打印全量纲计算结果到命令行（直接显示工程常用单位）
fprintf('\n=== 80kHz CW信号模糊度图-3dB轮廓全量纲计算结果 ===\n');
fprintf('等高线归一化幅度：%.3f（理论值0.707）\n', contour_level);
fprintf('时延范围 τ：%.2f ~ %.2f ms（教材理论值±%.2f ms）\n', ...
    tau_min_ms, tau_max_ms, 0.3*T*1e3);
fprintf('频移范围 ξ：%.1f ~ %.1f Hz（教材理论值±%.1f Hz）\n', ...
    xi_min_hz, xi_max_hz, 0.44/T);
fprintf('对应距离误差：%.2f ~ %.2f m（声速c=1500m/s）\n', ...
    tau_min_ms*1e-3*1500/2, tau_max_ms*1e-3*1500/2);
fprintf('对应速度误差：%.3f ~ %.3f m/s（波长λ=18.75mm）\n', ...
    xi_min_hz*18.75e-3/2, xi_max_hz*18.75e-3/2);
fprintf('===============================================\n');

% 在图上标注全量纲数值（红色加粗，工程常用格式）
text(tau_min_ms-0.03, 0, sprintf('%.2f', tau_min_ms), 'Color','red','FontSize',10,'FontWeight','bold');
text(tau_max_ms+0.01, 0, sprintf('%.2f', tau_max_ms), 'Color','red','FontSize',10,'FontWeight','bold');
text(0, xi_max_hz+20, sprintf('%.1f', xi_max_hz), 'Color','red','FontSize',10,'FontWeight','bold');
text(0, xi_min_hz-30, sprintf('%.1f', xi_min_hz), 'Color','red','FontSize',10,'FontWeight','bold');

%% ====================== 第二部分：LFM线性调频脉冲 ======================
fprintf('正在生成LFM脉冲信号...\n');

F = 10e3;         % 调频带宽 10kHz（可修改，建议BT≥10）
k = F/T;          % 调频斜率 Hz/s
B = F;            % 信号带宽（BT>>1时近似等于调频带宽）

% 1. 时间函数与瞬时频率（对应图2.7）
t_lfm = -T/2:dt:T/2-dt;  % 时间轴 [-T/2, T/2)
s_lfm = A * exp(1j * (2*pi*f0*t_lfm + pi*k*t_lfm.^2));  % LFM复信号
f_inst = f0 + k*t_lfm;  % 瞬时频率

figure('Name','图2.7 LFM脉冲信号波形与瞬时频率','Position',[100,100,800,600]);
subplot(2,1,1);
plot(t_lfm*1e3, real(s_lfm), 'LineWidth',1.2);
xlabel('时间 t (ms)','FontSize',12);
ylabel('Re[s(t)]','FontSize',12);
title('LFM脉冲信号波形（实部）','FontSize',14);
grid on;
axis tight;
ylim([-1.1*A, 1.1*A]);

subplot(2,1,2);
plot(t_lfm*1e3, f_inst/1e3, 'LineWidth',1.2, 'Color','r');
xlabel('时间 t (ms)','FontSize',12);
ylabel('瞬时频率 f(t) (kHz)','FontSize',12);
title('LFM脉冲信号瞬时频率','FontSize',14);
grid on;
axis tight;
ylim([(f0-F/2)/1e3-0.5, (f0+F/2)/1e3+0.5]);

% 2. 频谱函数（对应图2.8）
N_lfm = length(s_lfm);
f_lfm = (-Fs/2:Fs/N_lfm:Fs/2-Fs/N_lfm);
S_lfm = fftshift(fft(s_lfm)) * dt;

% 幅度归一化：除以最大值，峰值=1
amp_norm = abs(S_lfm)/max(abs(S_lfm));

figure('Name','图2.8 LFM脉冲信号频谱','Position',[100,100,800,400]);
plot(f_lfm/1e3, amp_norm, 'LineWidth',1.2);
xlabel('频率 f (kHz)','FontSize',12);
ylabel('|S(f)| 归一化','FontSize',12);
title('LFM脉冲信号频谱','FontSize',14);
grid on;

axis tight;
xlim([60, 100]);  % 显示完整带宽范围

% 3. 模糊函数计算（对应图2.9）
tau_range_lfm = linspace(-T, T, 201);
xi_range_lfm = linspace(-2*B, 2*B, 201);  % 频移范围 [-2B, 2B]
[TAU_LFM, XI_LFM] = meshgrid(tau_range_lfm, xi_range_lfm);
chi_lfm = zeros(size(TAU_LFM));

for i = 1:length(xi_range_lfm)
    for j = 1:length(tau_range_lfm)
        tau = tau_range_lfm(j);
        xi = xi_range_lfm(i);
        if abs(tau) > T
            chi_lfm(i,j) = 0;
        else
            arg = pi * (k*tau + xi) * (T - abs(tau));
            if arg == 0  % 处理分母为0的极限情况
                chi_lfm(i,j) = T - abs(tau);
            else
                chi_lfm(i,j) = abs(sin(arg)/arg * (T - abs(tau)));
            end
        end
    end
end

% (a) 模糊函数三维图
figure('Name','图2.9(a) LFM模糊函数三维图','Position',[100,100,800,600]);
surf(TAU_LFM*1e3, XI_LFM/B, chi_lfm/T, 'EdgeColor','none');
xlabel('时延 τ ','FontSize',12);
ylabel('频移 ξ ','FontSize',12);
zlabel('|χ(τ,ξ)|','FontSize',12);
title('LFM脉冲信号模糊函数','FontSize',14);
shading interp;
colormap jet;
view(30, 30);
colorbar;

% (b) ξ=0截面（时延截面）
chi_lfm_xi0 = zeros(size(tau_range_lfm));
for j = 1:length(tau_range_lfm)
    tau = tau_range_lfm(j);
    if abs(tau) > T
        chi_lfm_xi0(j) = 0;
    else
        arg = pi * k * tau * (T - abs(tau));
        if arg == 0
            chi_lfm_xi0(j) = T - abs(tau);
        else
            chi_lfm_xi0(j) = abs(sin(arg)/arg * (T - abs(tau)));
        end
    end
end  

figure('Name','图2.9(b) LFM模糊函数ξ=0截面','Position',[100,100,800,400]);
plot(tau_range_lfm*B, chi_lfm_xi0/T, 'LineWidth',1.2);  % ✅ 横轴：τ*B（无量纲）
xlabel('时延 \tau \cdot B','FontSize',12);
ylabel('|\chi(\tau,0)| / T','FontSize',12);
title('LFM模糊函数 \xi=0 截面','FontSize',14);
grid on;
hold on;

% ---------------- 关键修正 ----------------
tau_07_lfm = 0.44/B;               % 3dB 对应时延
x_07 = tau_07_lfm * B;             % ✅ 转换成无量纲 τ*B 坐标

% 画水平虚线 y=0.707
plot([-x_07, x_07], [0.707, 0.707], 'r--', 'LineWidth',1.2);

% 画左右两条竖虚线
plot([-x_07, -x_07], [0, 0.707], 'r--', 'LineWidth',1.2);
plot([x_07, x_07], [0, 0.707], 'r--', 'LineWidth',1.2);

% 文字标注（全部用无量纲 τ*B 坐标）
text(-x_07+0.1, 0.75, '0.707','Color','red','FontSize',10);
text(-x_07-1, -0.05, sprintf('-%.2f',x_07),'Color','red','FontSize',10);
text( x_07, -0.05, sprintf('%.2f',x_07), 'Color','red','FontSize',10);

hold off;
ylim([0, 1.05]);

% (c) kτ+ξ=0截面（斜截面）
chi_lfm_kxi0 = T - abs(tau_range_lfm);
figure('Name','图2.9(c) LFM模糊函数kτ+ξ=0截面','Position',[100,100,800,400]);
plot(tau_range_lfm*1e3, chi_lfm_kxi0/T, 'LineWidth',1.2);
xlabel('时延τ / T','FontSize',12);
ylabel('|χ(τ,ξ)| / T','FontSize',12);
title('LFM模糊函数 kτ+ξ=0 截面','FontSize',14);
grid on;
hold on;
% 标记0.707T对应的τ≈0.3T（教材标注）
tau_07_kxi0 = 0.3*T;
plot([-tau_07_kxi0*1e3, tau_07_kxi0*1e3], [0.707, 0.707], 'r--', 'LineWidth',1.2);
plot([-tau_07_kxi0*1e3, -tau_07_kxi0*1e3], [0, 0.707], 'r--', 'LineWidth',1.2);
plot([tau_07_kxi0*1e3, tau_07_kxi0*1e3], [0, 0.707], 'r--', 'LineWidth',1.2);
text(-tau_07_kxi0*1e3-0.1, 0.72, '0.707','Color','red','FontSize',10);
text(-tau_07_kxi0*1e3, -0.08, '-0.3','Color','red','FontSize',10);
text(tau_07_kxi0*1e3, -0.08, '0.3','Color','red','FontSize',10);
hold off;
axis tight;
ylim([0, 1.05]);

% (d) τ=0截面（频移截面）
chi_lfm_tau0 = abs(sin(pi*xi_range_lfm*T)./(pi*xi_range_lfm*T)) * T;
chi_lfm_tau0(xi_range_lfm==0) = T;
figure('Name','图2.9(d) LFM模糊函数τ=0截面','Position',[100,100,800,400]);
plot(xi_range_lfm * T, chi_lfm_tau0/T, 'LineWidth',1.2);  % ✅ 横轴：ξ·T（无量纲）
xlabel('频移 \xi \cdot T','FontSize',12);
ylabel('|\chi(0,\xi)| / T','FontSize',12);
title('LFM模糊函数 \tau=0 截面','FontSize',14);
grid on;
hold on;

% ---------------- 核心修正 ----------------
xi_07_lfm = 0.44/T;          % 3dB 频移点
x_07 = xi_07_lfm * T;        % ✅ 转换成无量纲坐标：ξ·T = 0.44

% 画水平虚线
plot([-x_07, x_07], [0.707, 0.707], 'r--', 'LineWidth',1.2);

% 画左右竖线
plot([-x_07, -x_07], [0, 0.707], 'r--', 'LineWidth',1.2);
plot([x_07, x_07], [0, 0.707], 'r--', 'LineWidth',1.2);

% 文字标注（全部用无量纲坐标）
text(-x_07-3, 0.72, '0.707','Color','red','FontSize',10);
text(-x_07-2, -0.05, sprintf('-%.2f',x_07),'Color','red','FontSize',10);
text( x_07, -0.05, sprintf('%.2f',x_07), 'Color','red','FontSize',10);

hold off;
axis tight;
ylim([0, 1.05]);

% (e) 模糊度图（等高线图）
figure('Name','图2.9(e) LFM脉冲模糊度图','Position',[100,100,600,600]);
% x=TAU_LFM(τ,s)，y=XI_LFM(ξ,Hz)
contour(TAU_LFM*1000, XI_LFM, chi_lfm/T, [0.707, 0.707], 'LineWidth',1.2);
xlabel('时延(ms) \tau','FontSize',12);
ylabel('频移 \xi','FontSize',12);
title('LFM脉冲信号模糊度图','FontSize',14);
grid on;

axis tight;
hold on;

% kτ+ξ=0斜虚线，原生坐标不变
tau_line = linspace(-T, T, 100);
xi_line = -k*tau_line;
plot(tau_line*1000, xi_line, 'k--', 'LineWidth',1.5);
text(0.1*T,0.1*B,'k\tau+\xi=0','Color','black','FontSize',10,'FontWeight','bold');
hold off;

fprintf('所有图形生成完成！\n');
fprintf('CW脉冲参数：中心频率%.1f kHz，脉宽%.1f ms，带宽%.1f kHz\n', f0/1e3, T*1e3, 1/T/1e3);
fprintf('LFM脉冲参数：中心频率%.1f kHz，脉宽%.1f ms，带宽%.1f kHz，BT=%.1f\n', f0/1e3, T*1e3, F/1e3, B*T);




% (e) 模糊度图（等高线图，全量纲自动计算-3dB轮廓数值）
figure('Name','图2.9(e) 80kHz LFM脉冲模糊度图','Position',[100,100,600,600]);
[C_LFM, h_LFM] = contour(TAU_LFM*1000, XI_LFM, chi_lfm/T, [0.707, 0.707], 'LineWidth',1.2);
xlabel('时延 τ (ms)','FontSize',12);  % 量纲：毫秒
ylabel('频移 ξ (Hz)','FontSize',12);  % 量纲：赫兹
title('80kHz LFM脉冲信号模糊度图（-3dB轮廓）','FontSize',14);
grid on;
axis tight;
hold on;

% ---------------- 核心：全量纲自动计算-3dB轮廓数值 ----------------
contour_level_lfm = C_LFM(1,1);
num_points_lfm = C_LFM(2,1);
x_coords_lfm = C_LFM(1, 2:num_points_lfm+1);  % τ，单位：ms
y_coords_lfm = C_LFM(2, 2:num_points_lfm+1);  % ξ，单位：Hz

tau_min_lfm_ms = min(x_coords_lfm);
tau_max_lfm_ms = max(x_coords_lfm);
xi_min_lfm_hz = min(y_coords_lfm);
xi_max_lfm_hz = max(y_coords_lfm);

% 打印全量纲计算结果到命令行
fprintf('\n=== 80kHz LFM信号模糊度图-3dB轮廓全量纲计算结果 ===\n');
fprintf('等高线归一化幅度：%.3f（理论值0.707）\n', contour_level_lfm);
fprintf('时延范围 τ：%.4f ~ %.4f ms（教材理论值±%.4f ms）\n', ...
    tau_min_lfm_ms, tau_max_lfm_ms, 0.44/B*1e3);
fprintf('频移范围 ξ：%.1f ~ %.1f Hz（教材理论值±%.1f Hz）\n', ...
    xi_min_lfm_hz, xi_max_lfm_hz, 0.44/T);
fprintf('对应距离误差：%.4f ~ %.4f m（声速c=1500m/s）\n', ...
    tau_min_lfm_ms*1e-3*1500/2, tau_max_lfm_ms*1e-3*1500/2);
fprintf('对应速度误差：%.3f ~ %.3f m/s（波长λ=18.75mm）\n', ...
    xi_min_lfm_hz*18.75e-3/2, xi_max_lfm_hz*18.75e-3/2);
fprintf('===============================================\n');

% 绘制kτ+ξ=0斜虚线（教材标注）
tau_line = linspace(-T, T, 100);
xi_line = -k*tau_line;
plot(tau_line*1000, xi_line, 'k--', 'LineWidth',1.5);
text(0.1*T*1000,0.1*B,'kτ+ξ=0','Color','black','FontSize',10,'FontWeight','bold');

% 在图上标注全量纲数值
text(tau_min_lfm_ms-0.002, 0, sprintf('%.4f', tau_min_lfm_ms), 'Color','red','FontSize',10,'FontWeight','bold');
text(tau_max_lfm_ms+0.001, 0, sprintf('%.4f', tau_max_lfm_ms), 'Color','red','FontSize',10,'FontWeight','bold');
text(0, xi_max_lfm_hz+20, sprintf('%.1f', xi_max_lfm_hz), 'Color','red','FontSize',10,'FontWeight','bold');
text(0, xi_min_lfm_hz-30, sprintf('%.1f', xi_min_lfm_hz), 'Color','red','FontSize',10,'FontWeight','bold');

hold off;




% (e) 模糊度图（等高线图，全量纲自动计算-3dB轮廓与坐标轴交点）
figure('Name','图2.9(e) 80kHz LFM脉冲模糊度图','Position',[100,100,600,600]);
[C_LFM, h_LFM] = contour(TAU_LFM*1000, XI_LFM, chi_lfm/T, [0.707, 0.707], 'LineWidth',1.2);
xlabel('时延 τ (ms)','FontSize',12);  % 量纲：毫秒
ylabel('频移 ξ (Hz)','FontSize',12);  % 量纲：赫兹
title('80kHz LFM脉冲信号模糊度图（-3dB轮廓）','FontSize',14);
grid on;
axis tight;
hold on;

% ---------------- 核心：精确计算-3dB轮廓与坐标轴的交点 ----------------
% 解析等高线矩阵C_LFM
contour_level_lfm = C_LFM(1,1);
num_points_lfm = C_LFM(2,1);
x_coords = C_LFM(1, 2:num_points_lfm+1);  % 所有点的x坐标（τ，ms）
y_coords = C_LFM(2, 2:num_points_lfm+1);  % 所有点的y坐标（ξ，Hz）

% ---------------- 函数：计算等高线与指定直线的交点 ----------------
function intersections = find_intersections(x, y, line_type, line_value)
    intersections = [];
    for i = 1:length(x)-1
        x1 = x(i); y1 = y(i);
        x2 = x(i+1); y2 = y(i+1);
        
        if strcmp(line_type, 'horizontal')  % 水平线 y=line_value
            if (y1 - line_value) * (y2 - line_value) <= 0
                t = (line_value - y1) / (y2 - y1);
                x_intersect = x1 + t * (x2 - x1);
                intersections = [intersections, x_intersect];
            end
        elseif strcmp(line_type, 'vertical')  % 垂直线 x=line_value
            if (x1 - line_value) * (x2 - line_value) <= 0
                t = (line_value - x1) / (x2 - x1);
                y_intersect = y1 + t * (y2 - y1);
                intersections = [intersections, y_intersect];
            end
        end
    end
end

% 1. 计算与x轴（ξ=0，即y=0）的交点 → 时延测量精度
tau_intersections = find_intersections(x_coords, y_coords, 'horizontal', 0);
tau_min_ms = min(tau_intersections);
tau_max_ms = max(tau_intersections);

% 2. 计算与y轴（τ=0，即x=0）的交点 → 频移测量精度
xi_intersections = find_intersections(x_coords, y_coords, 'vertical', 0);
xi_min_hz = min(xi_intersections);
xi_max_hz = max(xi_intersections);

% 打印全量纲计算结果到命令行
fprintf('\n=== 80kHz LFM信号模糊度图-3dB轮廓与坐标轴交点计算结果 ===\n');
fprintf('等高线归一化幅度：%.3f（理论值0.707）\n', contour_level_lfm);
fprintf('与x轴交点（时延τ）：%.4f ~ %.4f ms（教材理论值±%.4f ms）\n', ...
    tau_min_ms, tau_max_ms, 0.44/B*1e3);
fprintf('与y轴交点（频移ξ）：%.1f ~ %.1f Hz（教材理论值±%.1f Hz）\n', ...
    xi_min_hz, xi_max_hz, 0.44/T);
fprintf('对应距离误差：%.4f ~ %.4f m（声速c=1500m/s）\n', ...
    tau_min_ms*1e-3*1500/2, tau_max_ms*1e-3*1500/2);
fprintf('对应速度误差：%.3f ~ %.3f m/s（波长λ=18.75mm）\n', ...
    xi_min_hz*18.75e-3/2, xi_max_hz*18.75e-3/2);
fprintf('=======================================================\n');

