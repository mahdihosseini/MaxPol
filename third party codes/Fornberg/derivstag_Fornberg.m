function [c] = derivstag_Fornberg(l, s, n, flag_symbolic)

%STAGGERED_DERIVATIVES_FORNBERG 
%   Calls Fornberg's fullband staggered derivative coefficients
%
%   Copyright   Mahdi S. Hosseini, BSc, MSc, MASc, PhD
%               University of Toronto, November 2016
%               The Edward S. Rogers Sr. Department of,
%               Electrical and Computer Engineering (ECE)
%               Toronto, ON, M5S3G4, Canada
%               Tel: +1 (416) 978 6845
%               email: mahdi.hosseini@mail.utoronto.ca

if s >= 0
    c = fdcoeffF(n, -s, [-l+1:l]-1/2, flag_symbolic);
    
else
    c = fdcoeffF(n, -s, [-l+1:l]-1/2, flag_symbolic);
end
c = c(:)';

%%
if flag_symbolic
    c = eval(c);
end