clear all
close all
clc

%%
current_directory = pwd;
cd ..
cd ..
cd ..
addpath([cd, filesep, 'utilities'])
addpath([cd, filesep, 'data'])
addpath([cd, filesep, 'third party codes'])
addpath([cd, filesep, 'third party codes', filesep, 'Savitzky Golay'])
cd(current_directory)

%%
fullFileName = 't1_icbm_normal_1mm_pn3_rf20.mnc';
[imVOL,scaninfo] = loadminc(fullFileName);
N = size(imVOL);

%%  construct derivative matrices for each dimension
l = 5;  % tap-length polynomial
nod = 'staggered';  % node design
n_ord = 1;  % order of differentiation
sparse_flag = true; % generate derivative matrix in sparse mode
flag_symbolic = true; % numerically compute the fullband coefficients (very fast)
dim = 3;
P = 2;
[D_maxpol] = derivmtx(l, P, n_ord, N(dim), nod, sparse_flag, flag_symbolic);
P = 4;
[D_savitzkygolay] = derivmtx_SavitzkyGolay(l, P, n_ord, N(dim), nod, sparse_flag, flag_symbolic);

%%  Decompose the MRI 3D volume data in dim-dimension using maxpol
Dic{1} = speye(N(1));
Dic{2} = speye(N(2));
Dic{3} = speye(N(3));
Dic{dim} = D_maxpol;
[cVol_maxpol] = tensordec(imVOL, Dic);

%%  Decompose the MRI 3D volume data in dim-dimension using SacitzkyGolay
Dic{1} = speye(N(1));
Dic{2} = speye(N(2));
Dic{3} = speye(N(3));
Dic{dim} = D_savitzkygolay;
[cVol_savitzkygolay] = tensordec(imVOL, Dic);

%%
start_x = 60;
end_x   = 160;
start_y = 40;
end_y   = 140;
start_z = 5;
end_z   = 110;

%%
imVOL = imVOL(start_x: end_x, start_y: end_y, start_z: end_z);
cVol_maxpol = cVol_maxpol(start_x: end_x, start_y: end_y, start_z: end_z);
cVol_savitzkygolay = cVol_savitzkygolay(start_x: end_x, start_y: end_y, start_z: end_z);

%%
close all
cut_slice = 60;
c_image_maxpol = squeeze(cVol_maxpol(:, cut_slice, :));
c_image_savitzkygolay = squeeze(cVol_savitzkygolay(:, cut_slice, :));
slice_image = squeeze(imVOL(:, cut_slice, :));

range_c = [min([c_image_maxpol(:);c_image_savitzkygolay(:)]), ...
    max([c_image_maxpol(:);c_image_savitzkygolay(:)])];

figure('rend','painters','pos',[50, 200, 400, 400]);
img(slice_image)
title('xz-slice')
xlabel('z')
ylabel('x')
figure('rend','painters','pos',[475, 200, 400, 400]);
imagesc(c_image_maxpol, range_c)
title('\partial{z} - MaxPol')
xlabel('z')
ylabel('x')
axis image
colormap jet
set(gca, 'Xtick', [], 'Ytick', [])
figure('rend','painters','pos',[900, 200, 400, 400]);
imagesc(c_image_savitzkygolay, range_c)
title('\partial{z} - Savitzky-Golay')
xlabel('z')
ylabel('x')
axis image
colormap jet
set(gca, 'Xtick', [], 'Ytick', [])




