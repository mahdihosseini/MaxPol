function [c] = derivstag(l, P, s, n, sym_flag)
%
%DERIVSTAG Staggered FIR derivative coefficients
%
%   [c] = derivstag(l, P, s, n, sym_flag)
%
%   returns closed form solution to nth-order lowpass/fullband derivative
%   coefficeints at staggered nodes shifted by 's' in global scheme. The
%   associated coefficients are FIR filter and it is based on undetermined
%   coefficient and maxflat constraint design.
%
%   Convolution of discrete vecotr-valued function 'f' with 'c' estiamtes
%   the derivative vector d^nf/dx^n.
%
%   Input(s):
%   'l'         polynomial degree l>=1 (generates 2l-tap length FIR filter)
%   'P'         maxpol degree, controls cutoff frequency: n<=P<=2l-1
%   's'         shifting variable in staggered form:
%                   (1) zero shift: s=0
%                   (2) left-side shift: s>0
%                   (3) right-side shift: s<0
%   'n'         derivative order n>=0
%   'sym_flag'  numerical calculation: 'false', symbolic calculation: 'true'
%
%   Output(s):
%   'c'     (2l)-tap length FIR filter
%
%   See also DERIVCENT, DERIVMTX, DERIVDIREC
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

%% define shifting node 'b'
b=(-l+s-1/2);

%%
if (P < n)
    P = n;
    warning(['MaxPol degree should be n<=P<=2l-1. It is fixed to minimum: P=', num2str(n)])
elseif (P>= 2*l)
    warning(['MaxPol degree should be n<=P<=2l-1. It is fixed to maximum: P=', num2str(2*l-1)])
    P = 2*l-1;
end

%%
if (P == 2*l-1) & ~sym_flag
    %%  numerical calculation for Fullband filters (closed formulation)
    L = [1: 2*l];
    nod = b + L;
    C0 = C_fun(nod);
    for k = 1: 2*l
        c(k) = hypercube_corner_sum(n, C0([1: k-1, k+1: 2*l]));
    end
    c = prod(nod)*(-1).^(L+n+1)./nod./factorial(L-1)./factorial(2*l-L).*c;
else
    if sym_flag
        b = sym(b);
    end
    %% symbolic calculation
    V_b_P = repmat(b + [1: 2*l], [P+1, 1]);
    V_b_P = V_b_P .^ repmat([0: P]', [1, 2*l]);
    %
    if P < 2*l-1
        Q = 2*l - P - 2;
        V_b_Q = repmat(b + [1: 2*l], [Q+1, 1]);
        V_b_Q = V_b_Q .^ repmat([0: Q]', [1, 2*l]);
        %
        Lambda_s = diag((-1).^[1: 2*l]);
        V_b_Q = V_b_Q*Lambda_s;
    else
        V_b_Q = [];
    end
    %
    V_s = [V_b_P; V_b_Q];
    b = zeros(2*l, 1);
    b(n+1) = factorial(n);
    c = inv(V_s) * b;
    c = c(:)';
end

%%
if sym_flag
    c = eval(c);
end
