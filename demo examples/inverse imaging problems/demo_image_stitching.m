clear all
close all
clc

%%
current_directory = pwd;
cd ..
cd ..
addpath([cd, filesep, 'utilities'])
addpath([cd, filesep, 'data'])
cd(current_directory)

%%  initialized parameters
load('digital_pathology_preview_data.mat')
l = input('Please enter plynomial degree (l) for derivative matrix generation (between 1 to 7): ');
nod_number = input('Please define the derivative scheme, 0 for centralized, 1 for staggered: ');

%%  Lyapunov reconstruction call
if nod_number
    nod = 'staggered';
else
    nod = 'centralized';
end
[m,n,q] = size(i_enhanced_grad{1});
switch nod
    case 'staggered'
        D{1} = -derivmtx(l, 2*l-1, 1, n, 'staggered', true, false);
        D{2} = -derivmtx(l, 2*l-1, 1, m, 'staggered', true, false);
    case 'centralized'
        D{1} = -derivmtx(l, 2*l, 1, n, 'centralized', true, false);
        D{2} = -derivmtx(l, 2*l, 1, m, 'centralized', true, false);
end
i_rec = lyapunov_gradient_reconstruction(i_enhanced_grad, D);

%%  histogram ehancement
[i_rec] = histenhance(i_rec);
[i_enhanced] = histenhance(i_enhanced);

%%
i_rec_enhanced = i_rec - min(i_rec(:));
i_rec_enhanced = i_rec_enhanced/max(i_rec_enhanced(:))*255;

%%
stiched_tiles = stiched_tiles/max(stiched_tiles(:))*255;
i_enhanced = i_enhanced/max(i_enhanced(:))*255;

%%  raw image export with enhances image
figure('rend','painters','pos',[50, 100, n, m]);
img(stiched_tiles/255)
title('Raw stitching preview tiles')


%%
figure('rend','painters','pos',[n+100, 100, n, m]);
img(double(ROI_gain)/255)
title('Ilumination gain')

%%
figure('rend','painters','pos',[2*n+150, 100, n, m]);
img(i_enhanced/255)
title('White balance correction')

%%
figure('rend','painters','pos',[3*n+200, 100, n, m]);
img(i_rec_enhanced/255)
title(['Gradient surface reconstruction, ', nod, ', l=', num2str(l)])