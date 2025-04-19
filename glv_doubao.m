% ��������
N = 3; % ���˻�����
dt = 0.1; % ʱ�䲽��
T = 10; % ��ʱ��
t = 0:dt:T;
num_steps = length(t);

% ���˻���ʼλ��
x_uav = [0, 500, 1000];
y_uav = [0, 0, 0];

% ���˻��ٶ�
v_uav = 80;

% ������ʼλ��
x_ship = 1500;
y_ship = 500;

% �����ٶ�
v_ship = 50;
theta_ship = pi/4; % �����˶�����

% �״�����źŲ���
fc = 1e9; % �ز�Ƶ��
B = 10e6; % ����
Tchirp = 1e-3; % ������
k = B/Tchirp; % ��Ƶб��

% ��ʼ��λ������
x_uav_track = zeros(N, num_steps);
y_uav_track = zeros(N, num_steps);
x_ship_track = zeros(1, num_steps);
y_ship_track = zeros(1, num_steps);

% ģ�����˻��ͽ����˶����Ż����˻��˶�ģ��
for i = 1:num_steps
    for j = 1:N
        if j == 1
            % ��һ�����˻���x���˶�
            x_uav_track(j, i) = x_uav(j) + v_uav * dt * (i - 1);
            y_uav_track(j, i) = y_uav(j);
        elseif j == 2
            % �ڶ������˻���45�ȷ����˶�
            x_uav_track(j, i) = x_uav(j) + v_uav * cos(pi/4) * dt * (i - 1);
            y_uav_track(j, i) = y_uav(j) + v_uav * sin(pi/4) * dt * (i - 1);
        else
            % ���������˻���y���˶�
            x_uav_track(j, i) = x_uav(j);
            y_uav_track(j, i) = y_uav(j) + v_uav * dt * (i - 1);
        end
    end
    x_ship_track(i) = x_ship + v_ship * cos(theta_ship) * dt * (i - 1);
    y_ship_track(i) = y_ship + v_ship * sin(theta_ship) * dt * (i - 1);
end

% ���������Ϣ����Tikhonov�ֲ�
% ����ʱ����Ϣ���Ӹ�˹�ֲ�

% ģ���������
sigma_azimuth = 0.1; % ��λ�ǲ���������׼��
sigma_toa = 1e-6; % ����ʱ�����������׼��

% ��ʼ����λ���
x_est = zeros(1, num_steps);
y_est = zeros(1, num_steps);

for i = 1:num_steps
    % ������ʵ����ͷ�λ��
    for j = 1:N
        r(j) = sqrt((x_ship_track(i) - x_uav_track(j, i))^2 + (y_ship_track(i) - y_uav_track(j, i))^2);
        azimuth(j) = atan2(y_ship_track(i) - y_uav_track(j, i), x_ship_track(i) - x_uav_track(j, i));
    end
    
    % �����������
    azimuth_noisy = azimuth + sigma_azimuth * randn(1, N);
    toa_noisy = r / 3e8 + sigma_toa * randn(1, N);
    
    % ���ڸ����ܶȵ��ںϴ���
    % ����򻯴��������ȷֱ���ݲ����ʱ�����λ�ã����ں�
    % ���ڲ������λ��
    for j = 1:N
        x_azimuth_est(j) = x_uav_track(j, i) + r(j) * cos(azimuth_noisy(j));
        y_azimuth_est(j) = y_uav_track(j, i) + r(j) * sin(azimuth_noisy(j));
    end
    x_azimuth_mean = mean(x_azimuth_est);
    y_azimuth_mean = mean(y_azimuth_est);
    
    % ����ʱ�����λ��
    % ����ʹ��˫���߶�λԭ���ʵ��
    c = 3e8;
    % ȷ��toa_noisy��1��N������
    if size(toa_noisy, 1) ~= 1 || size(toa_noisy, 2) ~= N
        toa_noisy = reshape(toa_noisy, 1, N);
    end
    A = [2 * (x_uav_track(1, i) - x_uav_track(2, i)), 2 * (y_uav_track(1, i) - y_uav_track(2, i));
         2 * (x_uav_track(1, i) - x_uav_track(3, i)), 2 * (y_uav_track(1, i) - y_uav_track(3, i))];
    
    % ���򻯴�����ֹ��������
    lambda = 1e-6; 
    A_reg = A + lambda * eye(size(A));
    
    b = [(toa_noisy(2) - toa_noisy(1)) * c^2 + x_uav_track(2, i)^2 - x_uav_track(1, i)^2 + y_uav_track(2, i)^2 - y_uav_track(1, i)^2;
         (toa_noisy(3) - toa_noisy(1)) * c^2 + x_uav_track(3, i)^2 - x_uav_track(1, i)^2 + y_uav_track(3, i)^2 - y_uav_track(1, i)^2];
    % ���b������ά��
    if size(b, 1) ~= 2 || size(b, 2) ~= 1
        b = reshape(b, 2, 1);
    end
    xy_toa_est = A_reg \ b;
    x_toa_est = xy_toa_est(1);
    y_toa_est = xy_toa_est(2);
    
    % �ں϶�λ���
    x_est(i) = 0.5 * x_azimuth_mean + 0.5 * x_toa_est;
    y_est(i) = 0.5 * y_azimuth_mean + 0.5 * y_toa_est;
end

% ���㶨λ���
error_x = x_est - x_ship_track;
error_y = y_est - y_ship_track;
error = sqrt(error_x.^2 + error_y.^2);

% ���ƽ��
figure;
subplot(2,1,1);
plot(t, x_ship_track, 'b', 'DisplayName', '��ʵxλ��');
hold on;
plot(t, x_est, 'r--', 'DisplayName', '����xλ��');
xlabel('ʱ�� (s)');
ylabel('xλ�� (m)');
legend;
title('����xλ�ù���');

subplot(2,1,2);
plot(t, error, 'g', 'DisplayName', '��λ���');
xlabel('ʱ�� (s)');
ylabel('��λ��� (m)');
legend;
title('��λ�����ʱ��仯');