% sample
clear all; close all; clc;
% Initial Conditions
x_0 = [-2, -1, 0, 1, 2];

% Simulation time
t = 0:0.01:1;
t1 = -5:0.01:5;
f1 = t1;
f2 = cos(t1);

options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);

[T1, X1] = ode45(@x_func, t, x_0(1), options);
[T2, X2] = ode45(@x_func, t, x_0(2), options);
[T3, X3] = ode45(@x_func, t, x_0(3), options);
[T4, X4] = ode45(@x_func, t, x_0(4), options);
[T5, X5] = ode45(@x_func, t, x_0(5), options);

figure()
plot(T1, X1);
hold on;
plot(T2, X2); plot(T3, X3); plot(T4, X4); plot(T5, X5);
legend('I.C. x_0=-2', 'I.C. x_0=-1','I.C. x_0=0','I.C. x_0=1', 'I.C. x_0=2')
title('NonLinear Solution of x dot')
xlabel('Time');
ylabel('Solution of x dot');
grid on;

figure();
plot(t1,f1); hold on;
plot(t1,f2);
title('x and cos(x) Plots')
legend('f1 = x','f2 = cos(x)');
xlabel('time s');
ylabel('f(x)')
grid on;

k = sym('k');
p = solve(k-cos(k))

function dx = x_func(t, x)
dx = 0;
dx = x - cos(x);
end