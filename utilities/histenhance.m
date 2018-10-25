function [i] = histenhance(i)
%HISTENHANCE image histogram ehancement
%   Corrects the image histogram by canceling the outliers from lower and
%   higher gray values and balance the visual preview
%
%   Copyright   Mahdi S. Hosseini, BSc, MSc, MASc, PhD
%               University of Toronto, September 2015
%               The Edward S. Rogers Sr. Department of,
%               Electrical and Computer Engineering (ECE)
%               Toronto, ON, M5S3G4, Canada
%               Tel: +1 (416) 978 6845
%               email: mahdi.hosseini@mail.utoronto.ca


%% canceling high intensity outliers
[m, n, q] = size(i);
for channel = 1:q
    channel_image = i(:,:,channel);
    [N, x] = hist(channel_image(:), 100);
    max_energy_keeping = .999;
    max_point = cumsum(N)/sum(N);
    max_point = sum(max_point<max_energy_keeping);
    max_point = x(max_point);
    channel_image(channel_image>max_point) = max_point;
    i(:, :, channel) = channel_image;
end

%% canceling low intensity outliers
for channel = 1:q
    channel_image = i(:,:,channel);
    [N, x] = hist(channel_image(:), 100);
    min_energy_keeping = .999;
    min_point = cumsum(N(end:-1:1))/sum(N);
    min_point = sum(min_point<min_energy_keeping);
    min_point = x(numel(x)-min_point);
    channel_image(channel_image<min_point) = min_point;
    i(:, :, channel) = channel_image;
end