clear all
close all
clc

%% setup library directories
current_directory = pwd;
cd ..
cd ..
addpath([cd, filesep, 'utilities'])
addpath([cd, filesep, 'third party codes'])
cd(current_directory)

%%  setup parameters
possible_sizes_l = [2, 3:2:11]; % polynomial degree
N = 256;    % derivative matrix size
n_ord = 1;  % order of derivative

%%  marker type setup
marker_type = {'>', '<', '^', 'v', 'o', 'd'};

%%  Color specification
col = colormap(hot);
col = col(1:end-20, :);
col = flipud(col);
stp = length(col)/(numel(possible_sizes_l)+1);
col = col(round(1: stp: length(col)), :);
close all

%%  recall maxderiv (proposed) differential matrix for staggered coefficients with different matrix sizes
nod = 'staggered';
figure('rend','painters','pos',[50, 50, 424, 600]);
iteration = 0;
for l_polynomial = possible_sizes_l
    iteration = iteration + 1;
    staggered_polynomial_degree = l_polynomial
    %
    P = 2*l_polynomial-1;
    %
    [D_MaxPol_staggered{iteration}] = derivmtx(l_polynomial, P, n_ord, N, nod, false, false);
    maxpol_eigenvalues_staggered{iteration} = eig(-D_MaxPol_staggered{iteration});
    rank_MaxPol_staggered(iteration) = rank(D_MaxPol_staggered{iteration});
    
    %%
    range_of_eigenvlaues_staggered(iteration, :) = ...
        [min(real(maxpol_eigenvalues_staggered{iteration})),...
        max(real(maxpol_eigenvalues_staggered{iteration})), ...
        min(imag(maxpol_eigenvalues_staggered{iteration})), ...
        max(imag(maxpol_eigenvalues_staggered{iteration}))];
    %
    legend_string_MaxPol{iteration} = ['l=', num2str(l_polynomial)];
    %
    plot(real(maxpol_eigenvalues_staggered{iteration}), ...
        imag(maxpol_eigenvalues_staggered{iteration}), marker_type{iteration}, ...
        'MarkerEdgeColor', col(iteration, :),...
        'MarkerSize', 4)
    if iteration == 1
        hold on
    end
end

%%  recall maxderiv (proposed) differential matrix for staggered coefficients with different matrix sizes
nod = 'centralized';
figure('rend','painters','pos',[550, 50, 424, 600]);
iteration = 0;
for l_polynomial = possible_sizes_l
    iteration = iteration + 1;
    centralized_polynomial_degree = l_polynomial
    %
    P = 2*l_polynomial;
    %
    [D_MaxPol_centralized{iteration}] = derivmtx(l_polynomial, P, n_ord, N, nod, false, false);
    maxpol_eigenvalues_centralized{iteration} = eig(-D_MaxPol_centralized{iteration});
    rank_MaxPol_staggered(iteration) = rank(D_MaxPol_centralized{iteration});
    
    %%
    range_of_eigenvlaues_centralized(iteration, :) = ...
        [min(real(maxpol_eigenvalues_centralized{iteration})),...
        max(real(maxpol_eigenvalues_centralized{iteration})), ...
        min(imag(maxpol_eigenvalues_centralized{iteration})), ...
        max(imag(maxpol_eigenvalues_centralized{iteration}))];
    %
    legend_string_MaxPol{iteration} = ['l=', num2str(l_polynomial)];
    %
    plot(real(maxpol_eigenvalues_centralized{iteration}), ...
        imag(maxpol_eigenvalues_centralized{iteration}), marker_type{iteration}, ...
        'MarkerEdgeColor', col(iteration, :),...
        'MarkerSize', 4)
    if iteration == 1
        hold on
    end
end
%
range_of_eigenvlaues = [range_of_eigenvlaues_centralized;
    range_of_eigenvlaues_staggered];

ratio = (max(range_of_eigenvlaues(:,2))-min(range_of_eigenvlaues(:,1)))/...
    (max(range_of_eigenvlaues(:,4))-min(range_of_eigenvlaues(:,3)));
x_size = 5;
y_size = 5*ratio;

figure(1)
axis([min(range_of_eigenvlaues(:, 1)), ...
    max(range_of_eigenvlaues(:, 2)), ...
    min(range_of_eigenvlaues(:, 3)), ...
    max(range_of_eigenvlaues(:, 4))])
xlabel('Re(\lambda)')
ylabel('Im(\lambda)')
legend(legend_string_MaxPol, 'Location', 'NorthEast')
set(gca, 'FontSize', 12)
title('Staggered first order derivative matrix')

figure(2)
axis([min(range_of_eigenvlaues(:, 1)), ...
    max(range_of_eigenvlaues(:, 2)), ...
    min(range_of_eigenvlaues(:, 3)), ...
    max(range_of_eigenvlaues(:, 4))])
xlabel('Re(\lambda)')
ylabel('Im(\lambda)')
legend(legend_string_MaxPol, 'Location', 'NorthWest')
set(gca, 'FontSize', 12)
title('Centralized first order derivative matrix')
