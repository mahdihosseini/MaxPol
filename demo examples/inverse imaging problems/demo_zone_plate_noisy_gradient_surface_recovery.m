clear all
close all
clc

%%
current_directory = pwd;
cd ..
cd ..
addpath([cd, filesep, 'utilities'])
cd(current_directory)

%%  setup paramters for monte-carlo simualtion
frq = 1;
sgm = 0.02;
disp(['AWGN added to zone plate gradients = ', num2str(sgm*100), '%'])
selected_colormap = 'jet';
[f, grad, hess, third_order_tensor, fourth_order_tensor, params] = sinusoidal_gratting_2D(frq, 0);
N_x = params.N_x;
N_y = params.N_y;
random_perturbation_x = randn(N_y, N_x);
random_perturbation_y = randn(N_y, N_x);
f_0 = f - mean(f(:));
error_magnification = 1;

%% MaxPol - centralized l = 1
l = 1;
P = 2*l;
n_ord = 1;
nod = 'centralized';
sparse_flag = true;
sym_flag = false;
[D{1}] = derivmtx(l, P, n_ord, N_x, nod, sparse_flag, sym_flag);
[D{2}] = derivmtx(l, P, n_ord, N_y, nod, sparse_flag, sym_flag);
%
boundary_mask = logical(ones(N_y, N_x));
boundary_mask(l+1: end-l, l+1: end-l) = false;

%  setup centralized gradients
fx = grad.centralized{1};
fy = grad.centralized{2};
fx_N = fx + range(fx(:))*sgm*random_perturbation_x;
fy_N = fy + range(fy(:))*sgm*random_perturbation_y;
%
sylvester_grad{1} = fx_N;
sylvester_grad{2} = fy_N;
%
Dx = D{1}/params.h_x;
Dy = D{2}/params.h_y;
%
D{1} = Dx;
D{2}  = Dy;
disp('Graident surface reconstruction via solving Sylvester equation')
X = lyapunov_gradient_reconstruction(sylvester_grad, D);
X = X - mean(X(:));
X = X/norm(X)*norm(f_0);
%
res_frame = f_0-X;
NMSE = norm(abs(res_frame(:)))/norm(f_0);
[bin{1}, error_magnitude{1}]=hist(res_frame(:),50);


%  prepare for image excution
h_centralized_l1_recovery = figure('rend','painters','pos',[50, 50, 300, 300]);
imagesc(X, [min(f_0(:)), max(f_0(:))])
axis image
colormap('gray')
set(gca, 'Xtick', [], 'Ytick', [])
axis off
title('reconstructed image, centralized, l=1')

%
h_centralized_l1_residual_error = figure('rend','painters','pos',[50, 450, 300, 300]);
imagesc(abs(res_frame)*error_magnification, [0, max(f_0(:))-min(f_0(:))])
axis image
colormap(selected_colormap)
set(gca, 'Xtick', [], 'Ytick', [])
axis off
title(['error res., centralized, l=1', ', NMSE=', num2str(NMSE)])

%% MaxPol - centralized l = 5
l = 5;
P = 2*l;
n_ord = 1;
nod = 'centralized';
sparse_flag = true;
sym_flag = false;
[D{1}] = derivmtx(l, P, n_ord, N_x, nod, sparse_flag, sym_flag);
[D{2}] = derivmtx(l, P, n_ord, N_y, nod, sparse_flag, sym_flag);
%
boundary_mask = logical(ones(N_y, N_x));
boundary_mask(l+1: end-l, l+1: end-l) = false;

%  setup centralized gradients
fx = grad.centralized{1};
fy = grad.centralized{2};
fx_N = fx + range(fx(:))*sgm*random_perturbation_x;
fy_N = fy + range(fy(:))*sgm*random_perturbation_y;
%
sylvester_grad{1} = fx_N;
sylvester_grad{2} = fy_N;
%
Dx = D{1}/params.h_x;
Dy = D{2}/params.h_y;
%
D{1} = Dx;
D{2}  = Dy;
X = lyapunov_gradient_reconstruction(sylvester_grad, D);
X = X - mean(X(:));
X = X/norm(X)*norm(f_0);
%
res_frame = f_0-X;
NMSE = norm(abs(res_frame(:)))/norm(f_0);
[bin{2}, error_magnitude{2}]=hist(res_frame(:),50);

%  prepare for image excution
h_centralized_l5_recovery = figure('rend','painters','pos',[350, 50, 300, 300]);
imagesc(X, [min(f_0(:)), max(f_0(:))])
axis image
colormap('gray')
set(gca, 'Xtick', [], 'Ytick', [])
axis off
title('reconstructed image, centralized, l=5')

%
h_centralized_l5_residual_error = figure('rend','painters','pos',[350, 450, 300, 300]);
imagesc(abs(res_frame)*error_magnification, [0, max(f_0(:))-min(f_0(:))])
axis image
colormap(selected_colormap)
set(gca, 'Xtick', [], 'Ytick', [])
axis off
title(['error res., centralized, l=5', ', NMSE=', num2str(NMSE)])

%% MaxPol - staggered l = 1
l = 1;
P = 2*l-1;
n_ord = 1;
nod = 'staggered';
sparse_flag = true;
sym_flag = false;
[D{1}] = derivmtx(l, P, n_ord, N_x, nod, sparse_flag, sym_flag);
[D{2}] = derivmtx(l, P, n_ord, N_y, nod, sparse_flag, sym_flag);
%
boundary_mask = logical(ones(N_y, N_x));
boundary_mask(l+1: end-l, l+1: end-l) = false;

%  setup centralized gradients
fx = grad.staggered{1};
fy = grad.staggered{2};
fx_N = fx + range(fx(:))*sgm*random_perturbation_x;
fy_N = fy + range(fy(:))*sgm*random_perturbation_y;
%
sylvester_grad{1} = fx_N;
sylvester_grad{2} = fy_N;
%
Dx = D{1}/params.h_x;
Dy = D{2}/params.h_y;
%
D{1} = Dx;
D{2}  = Dy;
X = lyapunov_gradient_reconstruction(sylvester_grad, D);
X = X - mean(X(:));
X = X/norm(X)*norm(f_0);
%
res_frame = f_0-X;
NMSE = norm(abs(res_frame(:)))/norm(f_0);
[bin{3}, error_magnitude{3}]=hist(res_frame(:),50);

%  prepare for image excution
h_staggered_l1_recovery = figure('rend','painters','pos',[650, 50, 300, 300]);
imagesc(X, [min(f_0(:)), max(f_0(:))])
axis image
colormap('gray')
set(gca, 'Xtick', [], 'Ytick', [])
axis off
title('reconstructed image, staggered, l=1')

%
h_staggered_l1_residual_error = figure('rend','painters','pos',[650, 450, 300, 300]);
imagesc(abs(res_frame)*error_magnification, [0, max(f_0(:))-min(f_0(:))])
axis image
colormap(selected_colormap)
set(gca, 'Xtick', [], 'Ytick', [])
axis off
title(['error res., staggered, l=1', ', NMSE=', num2str(NMSE)])

%% MaxPol - staggered l = 5
l = 5;
P = 2*l-1;
n_ord = 1;
nod = 'staggered';
sparse_flag = true;
sym_flag = false;
[D{1}] = derivmtx(l, P, n_ord, N_x, nod, sparse_flag, sym_flag);
[D{2}] = derivmtx(l, P, n_ord, N_y, nod, sparse_flag, sym_flag);
%
boundary_mask = logical(ones(N_y, N_x));
boundary_mask(l+1: end-l, l+1: end-l) = false;

%  setup centralized gradients
fx = grad.staggered{1};
fy = grad.staggered{2};
fx_N = fx + range(fx(:))*sgm*random_perturbation_x;
fy_N = fy + range(fy(:))*sgm*random_perturbation_y;
%
sylvester_grad{1} = fx_N;
sylvester_grad{2} = fy_N;
%
Dx = D{1}/params.h_x;
Dy = D{2}/params.h_y;
%
% X = g2sSylvester(Dy, Dx, sylvester_grad{2}, sylvester_grad{1}, ...
%     ones(params.N_y,1), ones(params.N_x,1) ) ;
D{1} = Dx;
D{2}  = Dy;
X = lyapunov_gradient_reconstruction(sylvester_grad, D);
X = X - mean(X(:));
X = X/norm(X)*norm(f_0);
%
res_frame = f_0-X;
NMSE = norm(abs(res_frame(:)))/norm(f_0);
[bin{4}, error_magnitude{4}]=hist(res_frame(:),50);

%  prepare for image excution
h_staggered_l5_recovery = figure('rend','painters','pos',[950, 50, 300, 300]);
imagesc(X, [min(f_0(:)), max(f_0(:))])
axis image
colormap('gray')
set(gca, 'Xtick', [], 'Ytick', [])
axis off
title('reconstructed image, staggered, l=5')

%
h_staggered_l5_residual_error = figure('rend','painters','pos',[950, 450, 300, 300]);
imagesc(abs(res_frame)*error_magnification, [0, max(f_0(:))-min(f_0(:))])
axis image
colormap(selected_colormap)
set(gca, 'Xtick', [], 'Ytick', [])
axis off
title(['error res., staggered, l=5', ', NMSE=', num2str(NMSE)])