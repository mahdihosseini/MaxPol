clear all
close all
clc

%%
current_directory = pwd;
cd ..
cd ..
addpath([cd, filesep, 'utilities'])
cd(current_directory)

%%  lowpass filter with cutoff level P_x and P_y
n_ord = 1;  % order of differentiations to be evaluated
l = 13; % polynomial degree of the filter
p_x = 12;   % maxpol degree (cutoff) in x-axis
p_y = 12;   % maxpol degree (cutoff) in y-axis
theta = pi/4;   % rotation degree
N_fft = 513;
sym_flag = true;

%%  parameter for frequency response (fft) calculation
stp = 2*pi/(N_fft);
omega = [-pi: stp: pi-stp];

%%

[G] = derivdirec(l, p_x, p_y, n_ord, theta, sym_flag);
a = max(G(:));
b = min(G(:));
if a>abs(b)
    max_range = [-a, a];
else
    max_range = [b, -b];
end

%%  excute 2D fitler plot
figure('rend','painters','pos',[50 225 350 350]);
imagesc(G, max_range)
axis image
axis off
colormap gray
set(gca, 'Xtick', [], 'Ytick', [])
title(['l=',num2str(l), ', P_x=',num2str(p_x) , ', P_y=',num2str(p_y), ...
    ', derivative order n=', num2str(n_ord), ', \theta=', num2str(theta)])

%%  excute 2D fitler plot (Fourier Response)
F_G = fftshift(abs(fft2(G, N_fft, N_fft)));
a = max(F_G(:));
b = min(F_G(:));
if a>abs(b)
    max_range = [0, a];
else
    max_range = [0, -b];
end

%%
figure('rend','painters','pos',[450 225 350 350]);
imagesc(F_G, max_range)
axis image
axis off
colormap gray
set(gca, 'Xtick', [], 'Ytick', [])
title('2D Frequency response (magnitude)')

%%  contour level
figure('rend','painters','pos',[850 225 350 350]);
contour_levels = [1, .5, 10.^[-1 -2 -4 -7 -12]];
[c,h] = imcontour(omega, omega, F_G, contour_levels);
clabel(c,h)
axis image
set(gca, 'Xtick', [], 'Ytick', [])
title('Frequency magnitudes contour level')