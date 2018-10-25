function [K, bases] = derivdirec(l, P_x, P_y, n, theta, sym_flag)

%DERIVDIREC 2D direction (steerable) derivative kernel
%
%   [K, bases] = derivdirec(l, P_x, P_y, n, theta, sym_flag)
%
%   returns (2l+1)x(2l+1) two dimensional (2D) nth-order derivative kernel 
%   rotated at cerntain direction 'theta'. The associated derivative 
%   kernels are MaxPol based filters
%
%   Input(s):
%   'l'             order of tap-length polynomial
%   'n'             order of differentiation n>=0
%   'P'             maximally polynomial degree controling the cutoff
%                   threshold at x- or y- axes.
%   'theta'         steering angle
%   'sym_flag'      'true' for symbolic and 'false' for numerical
%                   calculation of maxpol coefficients
%
%   Output(s):
%   'K'             output steerable fitler
%   'bases'         decomposing bases of 'K'
%
%   See also DERIVSTAG, DERIVCENT, DERIVMTX
%
%
%   Copyright (c) 2017 Mahdi S. Hosseini
%
%   University of Toronto
%   The Edward S. Rogers Sr. Department of,
%   Electrical and Computer Engineering (ECE)
%   Toronto, ON, M5S3G4, Canada
%   Tel: +1 (416) 978 6845
%   email: mahdi.hosseini@mail.utoronto.ca

K = 0;
for k = 0: n
    %%  lowpass derivative filter with cutoff level (P_x, P_y)
    if mod(n - k, 2)
        d_x = derivstag(l, P_x, -1/2, n - k, sym_flag);
        d_x = [d_x, 0];
    else
        d_x = derivcent(l, P_x, 0, n - k, sym_flag);
        d_x = d_x;
    end
    %
    if mod(k, 2)
        d_y = derivstag(l, P_y, -1/2, k, sym_flag);
        d_y = [d_y, 0]';
    else
        d_y = derivcent(l, P_y, 0, k, sym_flag);
        d_y = d_y';
    end
    %
    %%  2D basis filter construction
    bases{k+1} = [d_y * d_x];
    %
    %%  bionomial coefficients
    c(k+1) = nchoosek(n, k) * cos(theta)^(n-k) * (-sin(theta))^k;
    %
    %%  accumulate 2D basis
    K = K + c(k+1) * bases{k+1};
end

for k = 0: n
    bases{k+1} = bases{k+1};
end