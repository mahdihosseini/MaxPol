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
frq = .8;   % harmonic frequency of 2D zone plate
delta = 0;
[f, grad, hess, third_order_tensor, fourth_order_tensor, params] = sinusoidal_gratting_2D(frq, 0);
N_x = params.N_x;
N_y = params.N_y;
h_x = params.h_x;
h_y = params.h_y;

%%  recall derivative matrices
possible_l = [1: 10];
iteration = 0;
for l = possible_l
    iteration = iteration + 1
    %
    boundary_mask = logical(ones(128, 128));
    boundary_mask(l+1: end-l, l+1: end-l) = false;
    %
    n_ord = 1;
    nod = 'staggered';
    P = 2*l-1;
    D_x = derivmtx(l, P, n_ord, N_x, nod, true, false);
    D_y = derivmtx(l, P, n_ord, N_y, nod, true, false);
    grad_approximation.staggered{l}{1} = f*D_x'/h_x^n_ord;
    grad_approximation.staggered{l}{2} = D_y*f/h_x^n_ord;
    hess_approximation.staggered{l}{2} = (D_y/h_y^n_ord)*f*D_x'/h_x^n_ord;
    %
    error_frame = grad.staggered{1}-grad_approximation.staggered{l}{1};
    error_weight = grad.staggered{1};
    err.staggered.grad.boundary(l, 1) = norm(error_frame(boundary_mask))/norm(error_weight(boundary_mask));
    err.staggered.grad.interior(l, 1) = norm(error_frame(~boundary_mask))/norm(error_weight(~boundary_mask));
    %
    error_frame = grad.staggered{2}-grad_approximation.staggered{l}{2};
    error_weight = grad.staggered{2};
    err.staggered.grad.boundary(l, 2) = norm(error_frame(boundary_mask))/norm(error_weight(boundary_mask));
    err.staggered.grad.interior(l, 2) = norm(error_frame(~boundary_mask))/norm(error_weight(~boundary_mask));
    %
    nod = 'centralized';
    P = 2*l;
    D_x = derivmtx(l, P, n_ord, N_x, nod, true, false);
    D_y = derivmtx(l, P, n_ord, N_y, nod, true, false);
    grad_approximation.centralized{l}{1} = f*D_x'/h_x^n_ord;
    grad_approximation.centralized{l}{2} = D_y*f/h_x^n_ord;
    hess_approximation.centralized{l}{2} = (D_y/h_y^n_ord)*f*D_x'/h_x^n_ord;
    %
    error_frame = grad.centralized{1}-grad_approximation.centralized{l}{1};
    error_weight = grad.centralized{1};
    err.centralized.grad.boundary(l, 1) = norm(error_frame(boundary_mask))/norm(error_weight(boundary_mask));
    err.centralized.grad.interior(l, 1) = norm(error_frame(~boundary_mask))/norm(error_weight(~boundary_mask));
    %
    error_frame = grad.centralized{2}-grad_approximation.centralized{l}{2};
    error_weight = grad.centralized{2};
    err.centralized.grad.boundary(l, 2) = norm(error_frame(boundary_mask))/norm(error_weight(boundary_mask));
    err.centralized.grad.interior(l, 2) = norm(error_frame(~boundary_mask))/norm(error_weight(~boundary_mask));
    %
    n_ord = 2;
    if l>1
        nod = 'staggered';
        P = 2*l-1;
        D_x = derivmtx(l, P, n_ord, N_x, nod, true, false);
        D_y = derivmtx(l, P, n_ord, N_y, nod, true, false);
        hess_approximation.staggered{l}{1} = f*D_x'/h_x^n_ord;
        hess_approximation.staggered{l}{3} = D_y*f/h_x^n_ord;
        %
        error_frame = hess.staggered{1}-hess_approximation.staggered{l}{1};
        error_weight = hess.staggered{1};
        err.staggered.hess.boundary(l, 1) = norm(error_frame(boundary_mask))/norm(error_weight(boundary_mask));
        err.staggered.hess.interior(l, 1) = norm(error_frame(~boundary_mask))/norm(error_weight(~boundary_mask));
        %
        error_frame = hess.staggered{2}-hess_approximation.staggered{l}{2};
        error_weight = hess.staggered{2};
        err.staggered.hess.boundary(l, 2) = norm(error_frame(boundary_mask))/norm(error_weight(boundary_mask));
        err.staggered.hess.interior(l, 2) = norm(error_frame(~boundary_mask))/norm(error_weight(~boundary_mask));
        %
        error_frame = hess.staggered{3}-hess_approximation.staggered{l}{3};
        error_weight = hess.staggered{3};
        err.staggered.hess.boundary(l, 3) = norm(error_frame(boundary_mask))/norm(error_weight(boundary_mask));
        err.staggered.hess.interior(l, 3) = norm(error_frame(~boundary_mask))/norm(error_weight(~boundary_mask));
    end
    %     %
    nod = 'centralized';
    P = 2*l;
    D_x = derivmtx(l, P, n_ord, N_x, nod, true, false);
    D_y = derivmtx(l, P, n_ord, N_y, nod, true, false);
    hess_approximation.centralized{l}{1} = f*D_x'/h_x^n_ord;
    hess_approximation.centralized{l}{3} = D_y*f/h_x^n_ord;
    %
    error_frame = hess.centralized{1}-hess_approximation.centralized{l}{1};
    error_weight = hess.centralized{1};
    err.centralized.hess.boundary(l, 1) = norm(error_frame(boundary_mask))/norm(error_weight(boundary_mask));
    err.centralized.hess.interior(l, 1) = norm(error_frame(~boundary_mask))/norm(error_weight(~boundary_mask));
    %
    error_frame = hess.centralized{2}-hess_approximation.centralized{l}{2};
    error_weight = hess.centralized{2};
    err.centralized.hess.boundary(l, 2) = norm(error_frame(boundary_mask))/norm(error_weight(boundary_mask));
    err.centralized.hess.interior(l, 2) = norm(error_frame(~boundary_mask))/norm(error_weight(~boundary_mask));
    %
    error_frame = hess.centralized{3}-hess_approximation.centralized{l}{3};
    error_weight = hess.centralized{3};
    err.centralized.hess.boundary(l, 3) = norm(error_frame(boundary_mask))/norm(error_weight(boundary_mask));
    err.centralized.hess.interior(l, 3) = norm(error_frame(~boundary_mask))/norm(error_weight(~boundary_mask));
end


%%  export results
rng = [err.centralized.hess.boundary(:); ...
    err.centralized.hess.interior(:);...
    err.staggered.hess.boundary(:); ...
    err.staggered.hess.interior(:); ...
    err.centralized.grad.boundary(:); ...
    err.centralized.grad.interior(:); ...
    err.staggered.grad.boundary(:); ...
    err.staggered.grad.interior(:)];
rng(rng==0)=[];
%
figure('rend','painters','pos',[1150, 150, 250, 500]);
bar(possible_l, err.centralized.grad.boundary(:,1), 'FaceColor', [0 .6 0])
hold on
bar(possible_l, err.centralized.grad.interior(:,1), 'FaceColor', [.7 .6 0])
set(gca,'YScale','log','FontSize',10)
legend('Boundary', 'Interior', 'Location', [.65 .85 .05 .05])
axis([0, max(possible_l)+1, min(rng), max(rng)])
xlabel('Tap-length - (l)')
ylabel('NMSE')
title('Centralized x-gradient approx.')

%
figure('rend','painters','pos',[1425, 150, 250, 500]);
bar(possible_l, err.staggered.grad.boundary(:,1), 'FaceColor', [0 .6 0])
hold on
bar(possible_l, err.staggered.grad.interior(:,1), 'FaceColor', [.7 .6 0])
set(gca,'YScale','log','FontSize',10)
legend('Boundary', 'Interior', 'Location', [.65 .85 .05 .05])
axis([0, max(possible_l)+1, min(rng), max(rng)])
xlabel('Tap-length - (l)')
ylabel('NMSE')
title('Staggered x-gradient approx.')

%
figure('rend','painters','pos',[50, 150, 250, 500]);
bar(possible_l, err.centralized.hess.boundary(:,1), 'FaceColor', [0 .6 0])
hold on
bar(possible_l, err.centralized.hess.interior(:,1), 'FaceColor', [.7 .6 0])
set(gca,'YScale','log','FontSize',10)
legend('Boundary', 'Interior', 'Location', [.65 .85 .05 .05])
axis([0, max(possible_l)+1, min(rng), max(rng)])
xlabel('Tap-length - (l)')
ylabel('NMSE')
title('Centralized xx-Hessian approx.')

%
figure('rend','painters','pos',[325, 150, 250, 500]);
bar(possible_l, err.centralized.hess.boundary(:,2), 'FaceColor', [0 .6 0])
hold on
bar(possible_l, err.centralized.hess.interior(:,2), 'FaceColor', [.7 .6 0])
set(gca,'YScale','log','FontSize',10)
legend('Boundary', 'Interior', 'Location', [.65 .85 .05 .05])
axis([0, max(possible_l)+1, min(rng), max(rng)])
xlabel('Tap-length - (l)')
ylabel('NMSE')
title('Centralized xy-Hessian approx.')

%
figure('rend','painters','pos',[600, 150, 250, 500]);
bar(possible_l, err.staggered.hess.boundary(:,1), 'FaceColor', [0 .6 0])
hold on
bar(possible_l, err.staggered.hess.interior(:,1), 'FaceColor', [.7 .6 0])
set(gca,'YScale','log','FontSize',10)
legend('Boundary', 'Interior', 'Location', [.65 .85 .05 .05])
axis([0, max(possible_l)+1, min(rng), max(rng)])
xlabel('Tap-length - (l)')
ylabel('NMSE')
title('Staggered xx-Hessian approx.')

%
figure('rend','painters','pos',[875, 150, 250, 500]);
bar(possible_l, err.staggered.hess.boundary(:,2), 'FaceColor', [0 .6 0])
hold on
bar(possible_l, err.staggered.hess.interior(:,2), 'FaceColor', [.7 .6 0])
set(gca,'YScale','log','FontSize',10)
legend('Boundary', 'Interior', 'Location', [.65 .85 .05 .05])
axis([0, max(possible_l)+1, min(rng), max(rng)])
xlabel('Tap-length - (l)')
ylabel('NMSE')
title('Staggered xy-Hessian approx.')

input('press enter to continue')
close all
%%  analytical representation of 2D sinusoidal grating
%
figure('rend','painters','pos',[200, 350, 300, 300]);
img(f)
axis off
title('2D Zone Plate')

%
figure('rend','painters','pos',[525, 50, 300, 300]);
img(grad.centralized{1})
axis off
title('Analytical x-gradient')

%
figure('rend','painters','pos',[525, 350, 300, 300]);
img(hess.centralized{2})
axis off
title('Analytical xy-Hessiant')

%
figure('rend','painters','pos',[525, 675, 300, 300]);
img(hess.centralized{1})
axis off
title('Analytical xx-Hessiant')

%%  close all
abs_dx_centralized_l1 = log10(abs(grad.centralized{1}-grad_approximation.centralized{1}{1})+1);
abs_dx_staggered_l1 = log10(abs(grad.staggered{1}-grad_approximation.staggered{1}{1})+1);
abs_dx_centralized_l6 = log10(abs(grad.centralized{1}-grad_approximation.centralized{6}{1})+1);
abs_dx_staggered_l6 = log10(abs(grad.staggered{1}-grad_approximation.staggered{6}{1})+1);
%
abs_dxy_centralized_l1 = log10(abs(hess.centralized{2}-hess_approximation.centralized{1}{2})+1);
abs_dxy_staggered_l1 = log10(abs(hess.staggered{2}-hess_approximation.staggered{1}{2})+1);
abs_dxy_centralized_l6 = log10(abs(hess.centralized{2}-hess_approximation.centralized{6}{2})+1);
abs_dxy_staggered_l6 = log10(abs(hess.staggered{2}-hess_approximation.staggered{6}{2})+1);
%
abs_dxx_centralized_l1 = log10(abs(hess.centralized{1}-hess_approximation.centralized{1}{1})+1);
abs_dxx_staggered_l1 = log10(abs(hess.staggered{1}-hess_approximation.staggered{2}{1})+1);
abs_dxx_centralized_l6 = log10(abs(hess.centralized{1}-hess_approximation.centralized{6}{1})+1);
abs_dxx_staggered_l6 = log10(abs(hess.staggered{1}-hess_approximation.staggered{6}{1})+1);
%
rng = [abs_dx_centralized_l1(:);...
    abs_dx_staggered_l1(:); ...
    abs_dx_centralized_l6(:); ...
    abs_dx_staggered_l6(:); ...
    abs_dxy_centralized_l1(:);...
    abs_dxy_staggered_l1(:); ...
    abs_dxy_centralized_l6(:); ...
    abs_dxy_staggered_l6(:); ...
    abs_dxx_centralized_l1(:);...
    abs_dxx_staggered_l1(:); ...
    abs_dxx_centralized_l6(:); ...
    abs_dxx_staggered_l6(:)];
%   x-gradients
figure('rend','painters','pos',[850, 50, 300, 300]);
imagesc(abs_dx_centralized_l1, [min(rng), max(rng)])
axis image
axis off
colormap jet
set(gca, 'Xtick', [], 'Ytick', [])
title('Abs error: x-gradient centralized, l=1')

%
figure('rend','painters','pos',[1175, 50, 300, 300]);
imagesc(abs_dx_staggered_l1, [min(rng), max(rng)])
axis image
axis off
colormap jet
set(gca, 'Xtick', [], 'Ytick', [])
title('Abs error: x-gradient staggered, l=1')

%
figure('rend','painters','pos',[1500, 50, 300, 300]);
imagesc(abs_dx_centralized_l6, [min(rng), max(rng)])
axis image
axis off
colormap jet
set(gca, 'Xtick', [], 'Ytick', [])
title('Abs error: x-gradient centralized, l=6')

%
figure('rend','painters','pos',[1825, 50, 300, 300]);
imagesc(abs_dx_staggered_l6, [min(rng), max(rng)])
axis image
axis off
colormap jet
set(gca, 'Xtick', [], 'Ytick', [])
title('Abs error: x-gradient staggered, l=6')

%   xy-hessian
figure('rend','painters','pos',[850, 350, 300, 300]);
imagesc(abs_dxy_centralized_l1, [min(rng), max(rng)])
axis image
axis off
colormap jet
set(gca, 'Xtick', [], 'Ytick', [])
title('Abs error: xy-Hessian centralized, l=1')

%
figure('rend','painters','pos',[1175, 350, 300, 300]);
imagesc(abs_dxy_staggered_l1, [min(rng), max(rng)])
axis image
axis off
colormap jet
set(gca, 'Xtick', [], 'Ytick', [])
title('Abs error: xy-Hessian staggered, l=1')

%
figure('rend','painters','pos',[1500, 350, 300, 300]);
imagesc(abs_dxy_centralized_l6, [min(rng), max(rng)])
axis image
axis off
colormap jet
set(gca, 'Xtick', [], 'Ytick', [])
title('Abs error: xy-Hessian centralized, l=6')

%
figure('rend','painters','pos',[1825, 350, 300, 300]);
imagesc(abs_dxy_staggered_l6, [min(rng), max(rng)])
axis image
axis off
colormap jet
set(gca, 'Xtick', [], 'Ytick', [])
title('Abs error: xy-Hessian staggered, l=6')

%   xx-hessian
figure('rend','painters','pos',[850, 675, 300, 300]);
imagesc(abs_dxx_centralized_l1, [min(rng), max(rng)])
axis image
axis off
colormap jet
set(gca, 'Xtick', [], 'Ytick', [])
title('Abs error: xx-Hessian centralized, l=1')

%
figure('rend','painters','pos',[1175, 675, 300, 300]);
imagesc(abs_dxx_staggered_l1, [min(rng), max(rng)])
axis image
axis off
colormap jet
set(gca, 'Xtick', [], 'Ytick', [])
title('Abs error: xx-Hessian staggered, l=1')

%
figure('rend','painters','pos',[1500, 675, 300, 300]);
imagesc(abs_dxx_centralized_l6, [min(rng), max(rng)])
axis image
axis off
colormap jet
set(gca, 'Xtick', [], 'Ytick', [])
title('Abs error: xx-Hessian centralized, l=6')

%
figure('rend','painters','pos',[1825, 675, 300, 300]);
imagesc(abs_dxx_staggered_l6, [min(rng), max(rng)])
axis image
axis off
colormap jet
set(gca, 'Xtick', [], 'Ytick', [])
title('Abs error: xx-Hessian staggered, l=6')
