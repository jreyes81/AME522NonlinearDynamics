% Homework 1
% AME 552 - Jeovanny Reyes

%% Problem 2.2.4
% xdot = sin(x) * exp(-x)
close all; clear all; clc

% Simulation time
t = 0:0.01:1;
xin = -10:0.01:10;
f1 = sin(xin);
f2 = exp(-xin);
f12 = f1.*f2;
% Initial Conditions
x_0 = [-10,-9,-8, -7, -6, -5, -4,-3,-2, -1, -0.5];

options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);

[T1, X1] = ode45(@x_func, t, x_0(1), options);
[T2, X2] = ode45(@x_func, t, x_0(2), options);
[T3, X3] = ode45(@x_func, t, x_0(3), options);
[T4, X4] = ode45(@x_func, t, x_0(4), options);
[T5, X5] = ode45(@x_func, t, x_0(5), options);
[T6, X6] = ode45(@x_func, t, x_0(6), options);
[T7, X7] = ode45(@x_func, t, x_0(7), options);
[T8, X8] = ode45(@x_func, t, x_0(8), options);
[T9, X9] = ode45(@x_func, t, x_0(9), options);
[T10, X10] = ode45(@x_func, t, x_0(10), options);
[T11, X11] = ode45(@x_func, t, x_0(11), options);

figure()
plot(xin,f12);
title('Phase Potrait of exp(-x) and sin(x)')
ylim([-20 180])
xlabel('x')
ylabel('xdot')
legend('f = exp(-x)*sin(x)')
grid on;

figure()
plot(T1, X1);
hold on;
plot(T2, X2); plot(T3, X3); plot(T4, X4); plot(T5, X5);plot(T6, X6);
plot(T7, X7); plot(T8, X8); plot(T9, X9); plot(T10, X10); plot(T11, X11);
legend('I.C. x_0=-10','I.C. x_0=-9','I.C. x_0=-8', 'I.C. x_0=-7','I.C. x_0=-6','I.C. x_0=-5', 'I.C. x_0=-4',.....
    'I.C. x_0=-3','I.C. x_0=-2','I.C. x_0=-1','I.C. x_0=-0.5');
title('NonLinear Solution of x dot = exp(-x)*sin(x)')
xlabel('Time');
ylabel('Solution of x dot');
grid on;

%% problem 2.2.7
% xdot = exp(x) - cos(x)
t2 = -17:0.01:17;
fa = exp(t2);
fb = cos(t2);
fab = fa-fb;
%Initial Conditions
x_1 = [-10,-9,-8, -7, -6, -5, -4,-3,-2, -1, -0.5, -0.1];

figure()
plot(t2,fa); hold on;
plot(t2,fb); plot(t2,fab);
title('Phase Potrait of exp(x) and cos(x)')
ylim([-2 2])
xlabel('x')
ylabel('xdot')
legend('fa=exp(x)','fb=cos(x)', 'fab = exp(x)-cos(x)')
grid on;

options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);

[Ta, Xa] = ode45(@x2_func, t2, x_1(1), options);
[Tb, Xb] = ode45(@x2_func, t2, x_1(2), options);
[Tc, Xc] = ode45(@x2_func, t2, x_1(3), options);
[Td, Xd] = ode45(@x2_func, t2, x_1(4), options);
[Te, Xe] = ode45(@x2_func, t2, x_1(5), options);
[Tf, Xf] = ode45(@x2_func, t2, x_1(6), options);
[Tg, Xg] = ode45(@x2_func, t2, x_1(7), options);
[Th, Xh] = ode45(@x2_func, t2, x_1(8), options);
[Ti, Xi] = ode45(@x2_func, t2, x_1(9), options);
[Tj, Xj] = ode45(@x2_func, t2, x_1(10), options);
[Tk, Xk] = ode45(@x2_func, t2, x_1(11), options);
[TL, XL] = ode45(@x2_func, t2, x_1(12), options);

figure()
plot(Ta, Xa);
hold on;
plot(Tb, Xb); plot(Tc, Xc); plot(Td, Xd); plot(Te, Xe); plot(Tf, Xf); plot(Tg, Xg);
plot(Th, Xh); plot(Ti, Xi); plot(Tj, Xj); plot(Tk, Xk); plot(TL, XL);
legend('I.C. x_0=-10','I.C. x_0=-9','I.C. x_0=-8', 'I.C. x_0=-7','I.C. x_0=-6','I.C. x_0=-5', 'I.C. x_0=-4',.....
    'I.C. x_0=-3','I.C. x_0=-2','I.C. x_0=-1','I.C. x_0=-0.5', 'I.C. x_0=-0.1');
title('NonLinear Solution of x dot = exp(x) - cos(x)')
xlabel('Time');
ylabel('Solution of x dot');
grid on;

%% Problem 3.1.1
% For the following exercise, sketch all the qualitatively different vector
% fields that occir as r is varied. Show that a saddle-node bifurcation
% occurs at a critical value of r, to be determined. Finally, sketch the
% bifuraction diagram of fixed points x* verus r

% xdot = 1+rx+x^2
r= [-3; -2; -1; 0; 1; 2; 3];
x = -5:0.01:5;
f_x1 = 1+r(1)*x+x.^2;
f_x2 = 1+r(2)*x+x.^2;
f_x3 = 1+r(3)*x+x.^2;
f_x4 = 1+r(4)*x+x.^2;
f_x5 = 1+r(5)*x+x.^2;
f_x6 = 1+r(6)*x+x.^2;
f_x7 = 1+r(7)*x+x.^2;

figure();
plot(x, f_x1); hold on
plot(x, f_x2);
plot(x, f_x3);
plot(x, f_x4);
plot(x, f_x5); plot(x, f_x6); plot(x, f_x7);
ylabel('xdot')
xlabel('x')
title('Plots of xdot = 1+r*x+x^2')
legend('r=-3','r=-2','r=-1', 'r=0','r=1', 'r=2', 'r=3');
grid on;

figure()
rvals = -6:0.01:6;
xstarPos = -(rvals + sqrt(rvals.^2 -4))/2;
xstarNeg = -(rvals - sqrt(rvals.^2 -4))/2;
plot(rvals, xstarPos);
hold on; plot(rvals, xstarNeg);
xlabel('r');
ylabel('xstar');
legend('xstarPos', 'xstarNeg');
grid on;

%% Problem 3.2.4
% For the following exercise, sketch all the qualitatively different vector
% fields that occir as r is varied. Show that a transcritical bifurcation
% occurs at a critical value of r, to be determined. Finally, sketch the
% bifuraction diagram of fixed points x* verus r

% xdot = x*(r-exp(x))
r= [-3; -2; -1; 0; 1; 2; 3];
x = -5:0.01:5;
f_x1 = x.*(r(1) - exp(x));
f_x2 = x.*(r(2) - exp(x));
f_x3 = x.*(r(3) - exp(x));
f_x4 = x.*(r(4) - exp(x));
f_x5 = x.*(r(5) - exp(x));
f_x6 = x.*(r(6) - exp(x));
f_x7 = x.*(r(7) - exp(x));

figure();
plot(x, f_x1); hold on
plot(x, f_x2);
plot(x, f_x3);
plot(x, f_x4);
plot(x, f_x5); plot(x, f_x6); plot(x, f_x7);
ylabel('xdot')
xlabel('x')
title('Plots of xdot = x*(r-exp(x))')
legend('r=-3','r=-2','r=-1', 'r=0','r=1', 'r=2', 'r=3');
grid on;

function dx = x_func(t, x)
dx = 0;
dx = sin(x)*exp(-x);
end

function dx = x2_func(t, x)
dx = 0;
dx = exp(x) - cos(x);
end