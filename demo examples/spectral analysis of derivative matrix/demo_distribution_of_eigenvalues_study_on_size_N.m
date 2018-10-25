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
possible_sizes_N = [11:5:200];
l_polynomial = 5;

%%  marker type setup
marker_type = {'>', '<', '^', 'v', 'o'};

%%  Color specification
col = colormap(hot);
col = col(1:end-20, :);
col = flipud(col);
stp = length(col)/(l_polynomial+1);
col = col(round(1: stp: length(col)), :);
close all

%%  recall maxderiv (proposed) differential matrix for staggered coefficients with different matrix sizes
P = 2*l_polynomial-1;
n_ord = 1;
nod = 'staggered';
iteration = 0;
for N = possible_sizes_N
    iteration = iteration + 1;
    matrix_size = N
    %
    [D_MaxPol_staggered{iteration}] = derivmtx(l_polynomial, P, n_ord, N, nod, false, false);
    maxpol_eigenvalues_staggered{iteration} = eig(D_MaxPol_staggered{iteration});
    rank_MaxPol_staggered(iteration) = rank(D_MaxPol_staggered{iteration});
    
    %%
    range_of_eigenvlaues_staggered(iteration, :) = ...
        [min(real(maxpol_eigenvalues_staggered{iteration})),...
        max(real(maxpol_eigenvalues_staggered{iteration})), ...
        min(imag(maxpol_eigenvalues_staggered{iteration})), ...
        max(imag(maxpol_eigenvalues_staggered{iteration}))];
end

%%  recall maxderiv (proposed) differential matrix for centralized coefficients with different matrix sizes
P = 2*l_polynomial;
n_ord = 1;
nod = 'centralized';
iteration = 0;
for N = possible_sizes_N
    iteration = iteration + 1;
    matrix_size = N
    %
    [D_MaxPol_centralized{iteration}] = derivmtx(l_polynomial, P, n_ord, N, nod, false, false);
    maxpol_eigenvalues_centralized{iteration} = eig(D_MaxPol_centralized{iteration});
    rank_MaxPol_centralized(iteration) = rank(D_MaxPol_centralized{iteration});
    
    %%
    range_of_eigenvlaues_centralized(iteration, :) = ...
        [min(real(maxpol_eigenvalues_centralized{iteration})),...
        max(real(maxpol_eigenvalues_centralized{iteration})), ...
        min(imag(maxpol_eigenvalues_centralized{iteration})), ...
        max(imag(maxpol_eigenvalues_centralized{iteration}))];
end

%%
marker_size = 5;
figure('rend','painters','pos',[50, 100, 350, 350]);
plot(possible_sizes_N, range_of_eigenvlaues_staggered(:, 1), 'v', ...
    'MarkerEdgeColor', [.7 0 0 ],...
    'MarkerSize', marker_size)
hold on
plot(possible_sizes_N, range_of_eigenvlaues_staggered(:, 2), '^', ...
    'MarkerEdgeColor', [.7 0 0],...
    'MarkerSize', marker_size)
plot(possible_sizes_N, range_of_eigenvlaues_staggered(:, 3), 'v', ...
    'MarkerEdgeColor', [.7 .7 0],...
    'MarkerSize', marker_size)
plot(possible_sizes_N, range_of_eigenvlaues_staggered(:, 4), '^', ...
    'MarkerEdgeColor', [.7 .7 0],...
    'MarkerSize', marker_size)
grid
axis([min(possible_sizes_N), max(possible_sizes_N), ...
    1.1*min([range_of_eigenvlaues_centralized(:);range_of_eigenvlaues_staggered(:)]), ...
    1.1*max([range_of_eigenvlaues_centralized(:);range_of_eigenvlaues_staggered(:)])])
xlabel('Matrix size - (N)')
ylabel('Magnitude - (log scale)')
legend_string = {'min(Re(\lambda))', ...
    'max(Re(\lambda))', ...
    'min(Im(\lambda))', ...
    'max(Im(\lambda))'};
legend(legend_string, 'Location', [.6 .55 .15 .15])
set(gca, 'FontSize',12)
symlog(gca,'y',-1.5)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
marker_size = 5;
figure('rend','painters','pos',[450, 100, 350, 350]);
plot(possible_sizes_N, range_of_eigenvlaues_centralized(:, 1), 'v', ...
    'MarkerEdgeColor', [.7 0 0 ],...
    'MarkerSize', marker_size)
hold on
plot(possible_sizes_N, range_of_eigenvlaues_centralized(:, 2), '^', ...
    'MarkerEdgeColor', [.7 0 0],...
    'MarkerSize', marker_size)
plot(possible_sizes_N, range_of_eigenvlaues_centralized(:, 3), 'v', ...
    'MarkerEdgeColor', [.7 .7 0],...
    'MarkerSize', marker_size)
plot(possible_sizes_N, range_of_eigenvlaues_centralized(:, 4), '^', ...
    'MarkerEdgeColor', [.7 .7 0],...
    'MarkerSize', marker_size)
grid
axis([min(possible_sizes_N), max(possible_sizes_N), ...
    1.1*min([range_of_eigenvlaues_centralized(:);range_of_eigenvlaues_staggered(:)]), ...
    1.1*max([range_of_eigenvlaues_centralized(:);range_of_eigenvlaues_staggered(:)])])
xlabel('Matrix size - (N)')
ylabel('Magnitude - (log scale)')
legend_string = {'min(Re(\lambda))', ...
    'max(Re(\lambda))', ...
    'min(Im(\lambda))', ...
    'max(Im(\lambda))'};
legend(legend_string, 'Location', [.6 .45 .15 .15])
set(gca, 'FontSize',12)
symlog(gca,'y',-1.5)