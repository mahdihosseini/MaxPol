function [c] = derivcent_Fornberg(l, s, n, sym_flag)

%CENTRALIZED_DERIVATIVES_FORNBERG 
%   Calls Fornberg's fullband centralized derivative coefficients
%
%   Copyright   Mahdi S. Hosseini, BSc, MSc, MASc, PhD
%               University of Toronto, November 2016
%               The Edward S. Rogers Sr. Department of,
%               Electrical and Computer Engineering (ECE)
%               Toronto, ON, M5S3G4, Canada
%               Tel: +1 (416) 978 6845
%               email: mahdi.hosseini@mail.utoronto.ca

if s >= 0
    c = fdcoeffF(n, -s, [-l:l], sym_flag);
    
else
    c = fdcoeffF(n, -s, [-l:l], sym_flag);
end
c = c(:)';

%%
if sym_flag
    c = eval(c);
end