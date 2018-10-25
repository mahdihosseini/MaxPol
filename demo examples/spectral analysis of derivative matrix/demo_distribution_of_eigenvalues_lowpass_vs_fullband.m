clear all
close all
clc

%% setup library directories
current_directory = pwd;
cd ..
cd ..
addpath([cd, filesep, 'utilities'])
addpath([cd, filesep, 'third party codes', filesep, 'Savitzky Golay'])
cd(current_directory)


%%  setup parameters
N = 256;
l_polynomial = 5;
P = 1:2: (2*l_polynomial-1);
Sig = [2.7:-.5:.7];
n_ord = 1;
nod = 'staggered';

%%  marker type setup
marker_type = {'>', '<', '^', 'v', 'o'};

%%  Color specification
col = colormap(hot);
col = col(1:end-20, :);
col = flipud(col);
stp = length(col)/(l_polynomial+1);
col = col(round(1: stp: length(col)), :);
close all

%%  recall maxderiv (proposed) differential matrix
figure('rend','painters','pos',[50, 50, 500, 400]);
iteration = 0;
for p = P    
    iteration = iteration + 1;
    maxpol_degree = p
    %
    [D_MaxPol(:,:,iteration)] = derivmtx(l_polynomial, p, n_ord, N, nod, false, true);
    maxpol_eigenvalues(:, iteration) = eig(-D_MaxPol(:, :, iteration));
    rank_MaxPol(iteration) = rank(D_MaxPol(:, :, iteration));
    %%  excute plots
    legend_string_MaxPol{iteration} = ['P=', num2str(p), ...
        ', r(D)=', num2str(rank_MaxPol(iteration))];
    plot(real(maxpol_eigenvalues(:, iteration)), ...
        imag(maxpol_eigenvalues(:, iteration)), marker_type{iteration}, ...
        'MarkerEdgeColor', col(iteration, :),...
        'MarkerSize', 5)
    if iteration == 1
        hold on
    end
end
%
real_lambda_maxpol = real(maxpol_eigenvalues);
imag_lambda_maxpol = imag(maxpol_eigenvalues);
MaxPol_axis_ratio = range(real_lambda_maxpol(:))/range(imag_lambda_maxpol(:));
y_size = 5;
x_size = y_size*MaxPol_axis_ratio;
%
axis([-2.5 0.1 -1 1])
xlabel('Re(\lambda)')
ylabel('Im(\lambda)')
legend(legend_string_MaxPol, 'Location', [.25 .445 .35 .15])
set(gca, 'FontSize', 10)
title(['MaxPol derivativem matrix D_', num2str(n_ord), ' eigenvalues, polynomial degree l = ', num2str(l_polynomial)])

%%  Savitzky-Golay
figure('rend','painters','pos',[750, 50, 600, 400]);
iteration = 0;
for p = P
    iteration = iteration + 1;
    savitzky_golay_degree = p
    %
    [D_SavitzkyGolay(:,:,iteration)] = derivmtx_SavitzkyGolay(l_polynomial, p, n_ord, N, nod, false, true);
    savitzkygolay_eigenvalues(:, iteration) = eig(-D_SavitzkyGolay(:, :, iteration));
    rank_SavitzyGolay(iteration) = rank(D_SavitzkyGolay(:, :, iteration));
    %
    %%  excute plots
    legend_string_SavtizkyGolay{iteration} = ['P=', num2str(p), ...
        ', r(D)=', num2str(rank_SavitzyGolay(iteration))];
    plot(real(savitzkygolay_eigenvalues(:, iteration)), ...
        imag(savitzkygolay_eigenvalues(:, iteration)), marker_type{iteration}, ...
        'MarkerEdgeColor', col(iteration, :),...
        'MarkerSize', 5)
    if iteration == 1
        hold on
    end
end
real_lambda_SavtizkyGolay = real(savitzkygolay_eigenvalues);
imag_lambda_SavtizkyGolay = imag(savitzkygolay_eigenvalues);
SavtizkyGolay_axis_ratio = range(real_lambda_SavtizkyGolay(:))/range(imag_lambda_SavtizkyGolay(:));
%
y_size = 5;
x_size = y_size*SavtizkyGolay_axis_ratio;
%
axis([-2.5 .6 -1 1])
xlabel('Re(\lambda)')
ylabel('Im(\lambda)')
set(gca, 'FontSize', 10)
legend(legend_string_MaxPol, 'Location', [.2 .445 .35 .15])
title(['Savtizky-Golay derivativem matrix D_', num2str(n_ord), ' eigenvalues, polynomial degree l = ', num2str(l_polynomial)])
