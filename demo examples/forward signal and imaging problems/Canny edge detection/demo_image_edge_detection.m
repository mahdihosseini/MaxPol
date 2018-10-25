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

%%  load image data
image_scan = imread('pirate.tif');
[m, n, q] = size(image_scan);
if q > 1
    image_scan = rgb2lab(image_scan);
    image_scan = image_scan(:,:,1);
end
image_scan = double(image_scan);

%% initialized parameters
l = 8;
%
sgm = sqrt(2);  % Gaussian
%
P_smoothing_Savitzky_Goaly = 0;
P_derivative_Savitzky_Goaly = 2;
%
P_smoothing_MaxPol = 0;
P_derivative_MaxPol = 2;

%%  Canny via Gaussian kernels (MATLAB recommendation)
[gaussKernel] = gaussian_derivatives(l, sgm);
smoothing_kernel = gaussKernel{1};
derivative_kernel = gaussKernel{2};
[segmented_Gaussian, dx_Gaussian, dy_Gaussian] = canny_edge(image_scan, smoothing_kernel, derivative_kernel);

%%  Canny via Savitzky-Golay
[smoothing_kernel] = derivcent_SavitzkyGolay(l, P_smoothing_Savitzky_Goaly, 0, 0, true);
[derivative_kernel] = derivstag_SavitzkyGolay(l, P_derivative_Savitzky_Goaly, 0, 1, true);
[segmented_Savitzky_Golay, dx_SavitzkyGolay, dy_SavitzkyGolay] = canny_edge(image_scan, smoothing_kernel, derivative_kernel);

%%  Canny via MaxPol
[smoothing_kernel] = derivcent(l, P_smoothing_MaxPol, 0, 0, true);
[derivative_kernel] = -derivstag(l, P_derivative_MaxPol, 0, 1, true);
[segmented_MaxPol, dx_MaxPol, dy_MaxPol] = canny_edge(image_scan, smoothing_kernel, derivative_kernel);

%%
figure('rend','painters','pos',[50, 50, 300, 300]);
img(image_scan/255)
title('Origianl image')

%%  excute Gaussian results
figure('rend','painters','pos',[400, 50, 300, 300]);
img(segmented_Gaussian)
title('Canny segmentation (Gaussian)')

figure('rend','painters','pos',[400, 400, 300, 300]);
img(dx_Gaussian)
title('x-gradient (Gaussian)')

figure('rend','painters','pos',[400, 750, 300, 300]);
img(dy_Gaussian)
title('y-gradient (Gaussian)')

%%  excute Savitzky-Golay results
figure('rend','painters','pos',[750, 50, 300, 300]);
img(segmented_Savitzky_Golay)
title('Canny segmentation (Savitzky-Golay)')

figure('rend','painters','pos',[750, 400, 300, 300]);
img(dx_SavitzkyGolay)
title('x-gradient (Savitzky-Golay)')

figure('rend','painters','pos',[750, 750, 300, 300]);
img(dy_SavitzkyGolay)
title('y-gradient (Savitzky-Golay)')

%%  excute MaxPol results
figure('rend','painters','pos',[1100, 50, 300, 300]);
img(segmented_MaxPol)
title('Canny segmentation (MaxPol)')

figure('rend','painters','pos',[1100, 400, 300, 300]);
img(dx_MaxPol)
title('x-gradient (MaxPol)')

figure('rend','painters','pos',[1100, 750, 300, 300]);
img(dy_MaxPol)
title('y-gradient (MaxPol)')


