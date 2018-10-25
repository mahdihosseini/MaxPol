clear all
close all
clc

%%
current_directory = pwd;
cd ..
cd ..
cd ..
addpath([cd, filesep, 'utilities'])
cd(current_directory)

%%
frequency = 2;    % Harmonic frequency
x_start = -1;
x_end = 1;
N = 128;    % number of discrete sampling
h_T = (x_end - x_start)/(N - 1);
x = [x_start: h_T: x_end]';
omega = 2*pi*frequency;

%%  creat analytical vector-valued functions
y = sin(omega*x);
dy = omega * cos(omega*x);
dy_half_sample_delay = omega * cos(omega*(x - h_T/2));

%%  approximate numerical derivative
l = 7;
P = 2*l - 1;
n_ord = 1;
nod = 'staggered';
sparse_flag = true;
sym_flag = false;
[D, D_forward] = derivmtx(l, P, n_ord, N, nod, sparse_flag, sym_flag);
dy_approximate = D*y/h_T^n_ord;

%%  error approximation
switch nod
    case 'staggered'
        dy_analytical = dy_half_sample_delay;
    case 'centralized'
        dy_analytical = dy;
end
err = abs(dy_analytical - dy_approximate);


%%  excute plots
figure('rend','painters','pos',[325, 700, 500, 250]);
plot(x, y)
title(['sinusoid signal y = sin(\omega x), here \omega = ', num2str(omega)])
xlabel('x')
ylabel('y')

figure('rend','painters','pos',[900, 700, 500, 250]);
plot(x, dy_analytical)
hold on
plot(x, dy_approximate, '--r')
title(['sinusoid derivative dy = \omega cos(\omega x)'])
legend('analytical', ['approximated, ', nod, ', l=', num2str(l)])
xlabel('x')
ylabel('y')

figure('rend','painters','pos',[900, 200, 500, 250]);
plot(x, err)
set(gca, 'YScale', 'log')
xlabel('x')
ylabel('|dy - D_{1}y|')