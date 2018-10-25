function [c] = derivcent(l, P, s, n, sym_flag)

%DERIVCENT Centralized FIR derivative coefficients
%
%   [c] = derivcent(l, P, s, n, sym_flag)
%
%   returns closed form solution to nth-order lowpass/fullband derivative
%   coefficeints at centralized nodes shifted by 's' in global scheme. The
%   associated coefficients are FIR filter and it is based on undetermined
%   coefficient and maxflat constraint design.
%
%   Convolution of discrete vecotr-valued function 'f' with 'c' estiamtes
%   the derivative vector d^nf/dx^n.
%
%   Input(s):
%   'l'         polynomial degree l>=1 (generates 2l+1-tap length FIR filter)
%   'P'         maxpol degree, controls cutoff frequency: n<=P<=2l
%   's'         shifting variable in staggered form:
%                   (1) zero shift: s=0
%                   (2) left-side shift: s>0
%                   (3) right-side shift: s<0
%   'n'         derivative order n>=0
%   'sym_flag'  numerical calculation: 'false', symbolic calculation: 'true'
%
%   Output(s):
%   'c'     (2l+1)-tap length FIR filter
%
%   See also DERIVSTAG, DERIVMTX, DERIVDIREC
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
a=-l+s-1;
if sym_flag
    a = sym(a);
end

%%
if (P < n)
    P = n;
    warning(['MaxPol degree should be n<=P<=2l. It is fixed to minimum: P=', num2str(n)])
elseif (P> 2*l)
    warning(['MaxPol degree should be n<=P<=2l. It is fixed to maximum: P=', num2str(2*l)])
    P = 2*l;
end

%%
if (P == 2*l) & ~sym_flag
    %%  numerical calculation for Fullband filters (closed formulation)
    L = [1: 2*l+1];
    nod = a + L;
    C0 = C_fun(nod);
    C0(-a) = 0;
    for k = 1: 2*l+1
        c(k) = hypercube_corner_sum(n-1, C0([1: k-1, k+1: 2*l+1]));
    end
    nod(-a) = 1;
    c = n*(-1).^(L+n+1)./factorial(L-1)./factorial(2*l+1-L).*prod(nod)./nod.*c;
    c(-a) = 0;
    if n > 0
        c(-a) = -sum(c);
    else
        c(-a) = 1;
    end
else
    %%  symbolic calculation for lowpass filters
    V_a_P = repmat(a + [1: 2*l+1], [P+1, 1]);
    V_a_P = V_a_P .^ repmat([0: P]', [1, 2*l+1]);
    %
    if P < 2*l
        Q = 2*l - P - 1;
        V_a_Q = repmat(a + [1: 2*l+1], [Q+1, 1]);
        V_a_Q = V_a_Q .^ repmat([0: Q]', [1, 2*l+1]);
        %
        Lambda_s = diag((-1).^[1: 2*l+1]);
        V_a_Q = V_a_Q*Lambda_s;
    else
        V_a_Q = [];
    end
    %
    V_s = [V_a_P; V_a_Q];
    b = zeros(2*l+1, 1);
    b(n+1) = factorial(n);
    c = inv(V_s) * b;
    c = c(:)';
end

%%
if sym_flag
    c = eval(c);
end