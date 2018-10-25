clear all
close all
clc

%% setup library directories
current_directory = pwd;
cd ..
cd ..
addpath([cd, filesep, 'utilities'])
cd(current_directory)

%%  setup parameters
N = 256;
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

%%  eigenvlaue analysis of staggered differential matrix
P = 2*l_polynomial-1;
nod = 'staggered';
n_ord = 1;
[D_MaxPol_staggered_1] = derivmtx(l_polynomial, P, n_ord, N, nod, false, false);
maxpol_eigenvalues_staggered_1 = eig(-D_MaxPol_staggered_1,'nobalance');
rank_MaxPol_staggered_1 = rank(D_MaxPol_staggered_1);
n_ord = 2;
D_MaxPol_staggered_2 = derivmtx(l_polynomial, P, n_ord, N, nod, false, false);
maxpol_eigenvalues_staggered_2 = eig(D_MaxPol_staggered_2,'nobalance');
rank_MaxPol_staggered_2 = rank(D_MaxPol_staggered_2);
%
P = 2*l_polynomial;
nod = 'centralized';
n_ord = 1;
D_MaxPol_centralized_1 = derivmtx(l_polynomial, P, n_ord, N, nod, false, false);
maxpol_eigenvalues_centralized_1 = eig(D_MaxPol_centralized_1,'nobalance');
rank_MaxPol_centralized_1 = rank(D_MaxPol_centralized_1);
n_ord = 2;
D_MaxPol_centralized_2 = derivmtx(l_polynomial, P, n_ord, N, nod, false, false);
maxpol_eigenvalues_centralized_2 = eig(D_MaxPol_centralized_2,'nobalance');
rank_MaxPol_centralized_2 = rank(D_MaxPol_centralized_2);
%

%%  range of eigenvalues determination
eigenvalue_range = [maxpol_eigenvalues_staggered_1(:);...
    maxpol_eigenvalues_centralized_1(:);...
    maxpol_eigenvalues_staggered_2(:);...
    maxpol_eigenvalues_centralized_2(:)];

%%
y_size = 5;
x_size = y_size*range(real(eigenvalue_range))/range(imag(eigenvalue_range));

%%
figure('rend','painters','pos',[50, 50, x_size*100, 500]);
plot(real(maxpol_eigenvalues_staggered_1), ...
    imag(maxpol_eigenvalues_staggered_1), ...
    'o', ...
    'MarkerEdgeColor', col(5, :),...
    'MarkerSize', 4)
hold on
plot(real(maxpol_eigenvalues_centralized_1), ...
    imag(maxpol_eigenvalues_centralized_1), ...
    'v', ...
    'MarkerEdgeColor', col(2, :),...
    'MarkerSize', 4)
%
xlabel('Re(\lambda)')
ylabel('Im(\lambda)')
set(gca, 'FontSize', 14)
axis([min(real(eigenvalue_range)), max(real(eigenvalue_range)), min(imag(eigenvalue_range)), max(imag(eigenvalue_range))])
legend(['Staggered,', 'P=', num2str(2*l_polynomial-1) ', r(D_1)=', num2str(rank_MaxPol_staggered_1)], ...
    ['Centralized,', 'P=', num2str(2*l_polynomial) ', r(D_1)=', num2str(rank_MaxPol_centralized_1)], 'Location', 'NorthWest')

%
figure('rend','painters','pos',[650, 50, x_size*100, 500]);
plot(real(maxpol_eigenvalues_staggered_2), ...
    imag(maxpol_eigenvalues_staggered_2), ...
    'o', ...
    'MarkerEdgeColor', col(5, :),...
    'MarkerSize', 4)
hold on
plot(real(maxpol_eigenvalues_centralized_2), ...
    imag(maxpol_eigenvalues_centralized_2), ...
    'v', ...
    'MarkerEdgeColor', col(2, :),...
    'MarkerSize', 4)
%
xlabel('Re(\lambda)')
ylabel('Im(\lambda)')
set(gca, 'FontSize', 14)
axis([min(real(eigenvalue_range)), max(real(eigenvalue_range)), min(imag(eigenvalue_range)), max(imag(eigenvalue_range))])
legend(['Staggered,', 'P=', num2str(2*l_polynomial-1) ', r(D_2)=', num2str(rank_MaxPol_staggered_2)], ...
    ['Centralized,', 'P=', num2str(2*l_polynomial) ', r(D_2)=', num2str(rank_MaxPol_centralized_2)], 'Location', 'NorthWest')
