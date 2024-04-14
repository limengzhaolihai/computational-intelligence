% 粒子群优化（PSO）算法求解无约束优化问题
%将f(x)拆解成g(x)详细描述见pdf文件
% 清除环境
clear all;
close all;
clc;

% 定义参数
n = 100; % 粒子数量
x_range = [-5.12, 5.12]; % x的取值范围
v_range = [-abs(diff(x_range)), abs(diff(x_range))]; % 速度范围
w = 0.9; % 惯性权重
c1 = 1.4; % 个体学习因子
c2 = 1.4; % 社会学习因子
max_iter = 100; % 最大迭代次数
tol = 1e-6; % 收敛容忍度

% 初始化粒子位置和速度
x = x_range(1) + (x_range(2) - x_range(1)) .* rand(n, 1);
v = v_range(1) + (v_range(2) - v_range(1)) .* rand(n, 1);

% 初始化个体最优位置和最优值
pbest = x;
pbest_val = zeros(n, 1); % 初始化个体最优值
for i = 1:n
    pbest_val(i) = objective_function(x(i));
end

% 初始化全局最优位置和最优值
[gbest_val, gbest_idx] = min(pbest_val);
gbest = x(gbest_idx);

% 保存适应度值的数组
fitness_curve = zeros(max_iter, 1);

% 迭代过程
for iter = 1:max_iter
    for i = 1:n
        % 更新速度和位置
        v(i) = w * v(i) ...
            + c1 * rand() * (pbest(i) - x(i)) ...
            + c2 * rand() * (gbest - x(i));
        x(i) = x(i) + v(i);

        % 限制位置在指定范围内
        x(i) = max(min(x(i), x_range(2)), x_range(1));

        % 更新个体最优位置和最优值
        current_val = objective_function(x(i));
        if current_val < pbest_val(i)
            pbest_val(i) = current_val;
            pbest(i) = x(i);
        end
        
        % 更新全局最优值
        if current_val < gbest_val
            gbest_val = current_val;
            gbest = x(i);
        end
    end

    % 记录适应度值
    fitness_curve(iter) = gbest_val;

    % 检查收敛条件
    if abs(gbest_val) < tol
        break;
    end
end

% 输出结果
disp(['最优解: x = ', num2str(gbest)]);
disp(['最优目标函数值: g(x) = ', num2str(gbest_val)]);


% 绘制适应度曲线
figure;
plot(1:iter, fitness_curve(1:iter), 'LineWidth', 2);
xlabel('迭代次数');
ylabel('适应度值');
title('适应度曲线');

% 目标函数
function f = objective_function(x)
    % 计算目标函数值
    f = (x.^2 - 10 * cos(2 * pi * x) + 20);
end
