% ==================== 均匀线列阵波束形成定位实验 ====================
% 内容：
% 1. 计算不同窗函数下的阵列指向性（方向图）
% 2. 模拟远场窄带信号，利用常规波束形成估计目标方向
% 3. 分析阵元数目 N 和窗函数对主瓣宽度、旁瓣抑制的影响
% 作者：不知周大王（AI辅助完成）
% ===================================================================

clear; clc; close all;

%% 1. 系统参数设置
c0 = 1500;          % 声速 (m/s)
f0 = 300;           % 信号频率 (Hz)
lambda = c0 / f0;   % 波长 (m)
d = lambda / 2;     % 阵元间距 = 半波长 (避免栅瓣)
% 本实验采用半波长间距，使波束形成效果明显

N_list = [11, 21, 31, 41];          % 待分析的阵元数目
angle_scan = -90:0.5:90;        % 扫描角度 (度)
theta_scan = angle_scan * pi/180;

%% 2. 窗函数定义 (用于后续方向图计算和波束形成)
win_type = {'rectwin', 'hamming', 'hanning', 'chebwin'};
win_color = {'b', 'r', 'g', 'm'};
win_name = {'矩形窗', '汉明窗', '汉宁窗', '切比雪夫窗(-40dB)'};

%% 3. 分析阵元数目 N 对指向性的影响 (采用矩形窗)
figure('Name', '阵元数目对指向性的影响', 'Position', [100 100 800 600]);
for i = 1:length(N_list)
    N = N_list(i);
    % 计算矩形窗的阵列因子 (归一化)
    w = ones(1, N);                         % 矩形窗
    AF = array_factor(w, d, lambda, theta_scan);
    
    subplot(2,2,i);
    plot(angle_scan, 20*log10(abs(AF)+eps), 'LineWidth', 1.5);
    grid on; xlim([-90 90]); ylim([-50 0]);
    xlabel('角度 (度)'); ylabel('归一化幅度 (dB)');
    title(sprintf('N = %d', N));
end
sgtitle('不同阵元数下的指向性 (矩形窗)');

%% 4. 分析不同窗函数对指向性的影响 (固定 N = 21)
N_fixed = 21;
figure('Name', '窗函数对指向性的影响', 'Position', [100 100 800 600]);
for i = 1:length(win_type)
    w = generate_window(win_type{i}, N_fixed);
    AF = array_factor(w, d, lambda, theta_scan);
    
    subplot(2,2,i);
    plot(angle_scan, 20*log10(abs(AF)+eps), win_color{i}, 'LineWidth', 1.5);
    grid on; xlim([-90 90]); ylim([-50 0]);
    xlabel('角度 (度)'); ylabel('归一化幅度 (dB)');
    title(win_name{i});
    hold on;
    % 叠加矩形窗的曲线作为对比 (灰色虚线)
    w_rect = ones(1, N_fixed);
    AF_rect = array_factor(w_rect, d, lambda, theta_scan);
    plot(angle_scan, 20*log10(abs(AF_rect)+eps), 'k--', 'LineWidth', 1);
    legend('当前窗', '矩形窗', 'Location', 'south');
end
sgtitle(sprintf('不同窗函数指向性对比 (N = %d)', N_fixed));

%% 5. 模拟接收信号，利用波束形成进行目标定位
% 参数设置
fs = 1000;          % 采样率 (Hz)
T = 1;              % 信号时长 (s)
t = 0:1/fs:T-1/fs;  % 时间向量
L = length(t);      % 快拍数

% 目标参数 (可设置多个目标)
targets = [30, 0.8;   % [方向(度), 幅度]
            -20, 0.5]; % 第二个目标
M_target = size(targets, 1);

% 生成阵列接收信号 (每个阵元接收多个平面波的叠加 + 噪声)
SNR_dB = 15;        % 信噪比 (dB)
N_sel = 21;         % 选用阵元数
d_sel = lambda/2;
array_pos = (0:N_sel-1) * d_sel;   % 阵元位置

% 生成信号源 (窄带复指数，简化为单频)
s = exp(1j*2*pi*f0*t);   % 参考信号 (幅度归一化)

% 构造接收数据矩阵 X (N阵元 × L快拍)
X = zeros(N_sel, L);
for k = 1:M_target
    angle_deg = targets(k,1);
    angle_rad = angle_deg * pi/180;
    amp = targets(k,2);
    % 各阵元相对于参考点的时延
    tau = array_pos * sin(angle_rad) / c0;
    % 相位延迟 (窄带近似)
    steering_vec = exp(-1j*2*pi*f0 * tau);
    X = X + amp * (steering_vec.' * s);
end

% 加入高斯白噪声
signal_power = mean(abs(X(:)).^2);
noise_power = signal_power / (10^(SNR_dB/10));
X = X + sqrt(noise_power/2) * (randn(N_sel, L) + 1j*randn(N_sel, L));

% 波束形成扫描，估计空间谱
theta_scan_deg = -90:0.5:90;
theta_scan_rad = theta_scan_deg * pi/180;
P_cbf = zeros(length(theta_scan_deg), 1);   % 常规波束形成输出功率
P_cbf_hamming = zeros(length(theta_scan_deg), 1);

% 不加窗 (矩形)
w_rect = ones(N_sel, 1);
% 加汉明窗
w_hamming = hamming(N_sel);

for i = 1:length(theta_scan_rad)
    % 扫描方向的导向矢量
    steering = exp(-1j*2*pi*f0 * array_pos' * sin(theta_scan_rad(i)) / c0);
    % 波束形成输出 (y = w^H * X)
    y_rect = w_rect' * X;
    y_hamming = w_hamming' * X;
    % 功率 (时间平均)
    P_cbf(i) = mean(abs(y_rect).^2);
    P_cbf_hamming(i) = mean(abs(y_hamming).^2);
end



%% 6. 辅助函数
% ------------------------------------------------------------------
% 计算均匀线列阵的阵列因子 (归一化)
% 输入: w - 权系数向量 (1×N 或 N×1)
%       d - 阵元间距
%       lambda - 波长
%       theta_scan - 扫描角度向量 (弧度)
% 输出: AF - 阵列因子 (与theta_scan同维度)
function AF = array_factor(w, d, lambda, theta_scan)
    N = length(w);
    w = w(:).';                     % 转为行向量
    k = 2*pi/lambda;               % 波数
    % 构造阵列流形矩阵 (N × M), M = length(theta_scan)
    n = (0:N-1)';
    steering = exp(1j * k * d * n * sin(theta_scan));
    AF = (w * steering) / sum(w);   % 归一化
end

% 生成指定类型的窗函数 (返回行向量)
function w = generate_window(type, N)
    switch lower(type)
        case 'rectwin'
            w = ones(1, N);
        case 'hamming'
            w = hamming(N)';
        case 'hanning'
            w = hanning(N)';
        case 'chebwin'
            % 切比雪夫窗，旁瓣电平设为 -40 dB
            w = chebwin(N, 40)';
        otherwise
            error('未知窗类型');
    end
end