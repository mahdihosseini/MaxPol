function [c] = derivstag_SavitzkyGolay(l, P, s, n, flag_symbolic)

%STAGGERED_DERIVATIVES Staggered derivative coefficients
%   [c] = staggered_derivatives(l, b, n) returns closed form solution to
%   the n-th order derivative coefficeints at staggered nodes shifted by
%   s-1/2 in global scheme. The associated coefficients are in FIR filter
%   form and comply with high-order-of-accuracy polynomial design.
%   Convolution of discrete vecotr-values function 'f' with 'c' estiamtes
%   the derivative vector d^nf/dx^n.
%
%   Input(s):
%   'l'     polynomial order l>=1 (generating 2l-tap length FIR filter)
%   's'     shifting variable in staggered form:
%           (1) zero shift: s=0
%           (2) left-side shift: s>0
%           (3) right-side shift: s<0
%   'n'     derivative order n>=1
%
%   Output(s):
%   'c'     2l tap-length derivative FIR filter
%
%   See also CENTRALIZED_DERIVATIVES, DERVMTX
%
%   Copyright   Mahdi S. Hosseini, BSc, MSc, MASc, PhD
%               University of Toronto, November 2016
%               The Edward S. Rogers Sr. Department of,
%               Electrical and Computer Engineering (ECE)
%               Toronto, ON, M5S3G4, Canada
%               Tel: +1 (416) 978 6845
%               email: mahdi.hosseini@mail.utoronto.ca

if s >= 0
    c = savitzkyGolay([-l+1:l]-1/2, P, n, -s, [], flag_symbolic);
else
    c = savitzkyGolay([-l+1:l]+1/2, P, n, -s, [], flag_symbolic);
end
if flag_symbolic
    c = eval(c(:))';
else
    c = c(:)';
end