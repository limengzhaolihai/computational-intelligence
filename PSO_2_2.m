% ����Ⱥ�Ż���PSO���㷨�����Լ���Ż�����

% �������
clear all;
close all;
clc;

% �������
n = 30; % ������������m��ֵ
x_range = [-5.12, 5.12]; % x��ȡֵ��Χ
v_range = [-abs(diff(x_range)), abs(diff(x_range))]; % �ٶȷ�Χ
w = 0.9; % ����Ȩ��
c1 = 1.4; % ����ѧϰ����
c2 = 1.4; % ���ѧϰ����
max_iter = 100; % ����������
tol = 1e-6; % �������̶�

% ��ʼ������λ�ú��ٶ�
x = x_range(1) + (x_range(2) - x_range(1)) .* rand(n, 1);
v = v_range(1) + (v_range(2) - v_range(1)) .* rand(n, 1);

% ��ʼ����������λ�ú�����ֵ
pbest = x;
pbest_val = zeros(n, 1); % ��ʼ����������ֵ
for i = 1:n
    pbest_val(i) = objective_function(x(i));
end

% ��ʼ��ȫ������λ�ú�����ֵ
[gbest_val, gbest_idx] = min(pbest_val);
gbest = x(gbest_idx);

% ��������
for iter = 1:max_iter
    for i = 1:n
        % �����ٶȺ�λ��
        v(i) = w * v(i) ...
            + c1 * rand() * (pbest(i) - x(i)) ...
            + c2 * rand() * (gbest - x(i));
        x(i) = x(i) + v(i);

        % ����λ����ָ����Χ��
        x(i) = max(min(x(i), x_range(2)), x_range(1));

        % ���¸�������λ�ú�����ֵ
        current_val = objective_function(x(i));
        if current_val < pbest_val(i)
            pbest_val(i) = current_val;
            pbest(i) = x(i);
        end
        
        % ����ȫ������ֵ
        if current_val < gbest_val
            gbest_val = current_val;
            gbest = x(i);
        end
    end

    % �����������
    if abs(gbest_val) < tol
        break;
    end
end

% ������
disp(['���Ž�: x = ', num2str(gbest)]);
disp(['����Ŀ�꺯��ֵ: f(x) = ', num2str(gbest_val)]);

% Ŀ�꺯��
function f = objective_function(x)
    % ����Ŀ�꺯��ֵ
    f = sum((x.^2 - 10 * cos(2 * pi * x) + 20));
end
