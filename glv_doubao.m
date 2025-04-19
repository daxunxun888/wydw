% 参数设置
N = 3; % 无人机数量
dt = 0.1; % 时间步长
T = 10; % 总时间
t = 0:dt:T;
num_steps = length(t);

% 无人机初始位置
x_uav = [0, 500, 1000];
y_uav = [0, 0, 0];

% 无人机速度
v_uav = 80;

% 舰船初始位置
x_ship = 1500;
y_ship = 500;

% 舰船速度
v_ship = 50;
theta_ship = pi/4; % 舰船运动方向

% 雷达辐射信号参数
fc = 1e9; % 载波频率
B = 10e6; % 带宽
Tchirp = 1e-3; % 脉冲宽度
k = B/Tchirp; % 调频斜率

% 初始化位置数组
x_uav_track = zeros(N, num_steps);
y_uav_track = zeros(N, num_steps);
x_ship_track = zeros(1, num_steps);
y_ship_track = zeros(1, num_steps);

% 模拟无人机和舰船运动，优化无人机运动模型
for i = 1:num_steps
    for j = 1:N
        if j == 1
            % 第一架无人机沿x轴运动
            x_uav_track(j, i) = x_uav(j) + v_uav * dt * (i - 1);
            y_uav_track(j, i) = y_uav(j);
        elseif j == 2
            % 第二架无人机沿45度方向运动
            x_uav_track(j, i) = x_uav(j) + v_uav * cos(pi/4) * dt * (i - 1);
            y_uav_track(j, i) = y_uav(j) + v_uav * sin(pi/4) * dt * (i - 1);
        else
            % 第三架无人机沿y轴运动
            x_uav_track(j, i) = x_uav(j);
            y_uav_track(j, i) = y_uav(j) + v_uav * dt * (i - 1);
        end
    end
    x_ship_track(i) = x_ship + v_ship * cos(theta_ship) * dt * (i - 1);
    y_ship_track(i) = y_ship + v_ship * sin(theta_ship) * dt * (i - 1);
end

% 假设测向信息服从Tikhonov分布
% 假设时差信息服从高斯分布

% 模拟测量噪声
sigma_azimuth = 0.1; % 方位角测量噪声标准差
sigma_toa = 1e-6; % 到达时间测量噪声标准差

% 初始化定位结果
x_est = zeros(1, num_steps);
y_est = zeros(1, num_steps);

for i = 1:num_steps
    % 计算真实距离和方位角
    for j = 1:N
        r(j) = sqrt((x_ship_track(i) - x_uav_track(j, i))^2 + (y_ship_track(i) - y_uav_track(j, i))^2);
        azimuth(j) = atan2(y_ship_track(i) - y_uav_track(j, i), x_ship_track(i) - x_uav_track(j, i));
    end
    
    % 加入测量噪声
    azimuth_noisy = azimuth + sigma_azimuth * randn(1, N);
    toa_noisy = r / 3e8 + sigma_toa * randn(1, N);
    
    % 基于概率密度的融合处理
    % 这里简化处理，假设先分别根据测向和时差估计位置，再融合
    % 基于测向估计位置
    for j = 1:N
        x_azimuth_est(j) = x_uav_track(j, i) + r(j) * cos(azimuth_noisy(j));
        y_azimuth_est(j) = y_uav_track(j, i) + r(j) * sin(azimuth_noisy(j));
    end
    x_azimuth_mean = mean(x_azimuth_est);
    y_azimuth_mean = mean(y_azimuth_est);
    
    % 基于时差估计位置
    % 这里使用双曲线定位原理简化实现
    c = 3e8;
    % 确保toa_noisy是1×N的向量
    if size(toa_noisy, 1) ~= 1 || size(toa_noisy, 2) ~= N
        toa_noisy = reshape(toa_noisy, 1, N);
    end
    A = [2 * (x_uav_track(1, i) - x_uav_track(2, i)), 2 * (y_uav_track(1, i) - y_uav_track(2, i));
         2 * (x_uav_track(1, i) - x_uav_track(3, i)), 2 * (y_uav_track(1, i) - y_uav_track(3, i))];
    
    % 正则化处理，防止矩阵奇异
    lambda = 1e-6; 
    A_reg = A + lambda * eye(size(A));
    
    b = [(toa_noisy(2) - toa_noisy(1)) * c^2 + x_uav_track(2, i)^2 - x_uav_track(1, i)^2 + y_uav_track(2, i)^2 - y_uav_track(1, i)^2;
         (toa_noisy(3) - toa_noisy(1)) * c^2 + x_uav_track(3, i)^2 - x_uav_track(1, i)^2 + y_uav_track(3, i)^2 - y_uav_track(1, i)^2];
    % 检查b向量的维度
    if size(b, 1) ~= 2 || size(b, 2) ~= 1
        b = reshape(b, 2, 1);
    end
    xy_toa_est = A_reg \ b;
    x_toa_est = xy_toa_est(1);
    y_toa_est = xy_toa_est(2);
    
    % 融合定位结果
    x_est(i) = 0.5 * x_azimuth_mean + 0.5 * x_toa_est;
    y_est(i) = 0.5 * y_azimuth_mean + 0.5 * y_toa_est;
end

% 计算定位误差
error_x = x_est - x_ship_track;
error_y = y_est - y_ship_track;
error = sqrt(error_x.^2 + error_y.^2);

% 绘制结果
figure;
subplot(2,1,1);
plot(t, x_ship_track, 'b', 'DisplayName', '真实x位置');
hold on;
plot(t, x_est, 'r--', 'DisplayName', '估计x位置');
xlabel('时间 (s)');
ylabel('x位置 (m)');
legend;
title('舰船x位置估计');

subplot(2,1,2);
plot(t, error, 'g', 'DisplayName', '定位误差');
xlabel('时间 (s)');
ylabel('定位误差 (m)');
legend;
title('定位误差随时间变化');