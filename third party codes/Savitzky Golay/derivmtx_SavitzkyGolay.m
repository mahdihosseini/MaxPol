function [D, D_forward] = derivmtx_SavitzkyGolay(l, P, n, N, nod, sparse_flag, sym_flag)

%DERIVMTX_SAVITZKYGOLAY Savitzky-Golay Derivative square matrix
%
%   [D, D_forward] = derivmtx_SavitzkyGolay(l, P, n, N, nod, sparse_flag, sym_flag)
%
%   returns NxN derivative matrix incorporated by Savitzky-Golay FIR
%   derivative kernels at every row of the matrix. d^nf = D*f estimates
%   n'th order derivative of discrete vector-vlaued function 'f'.
%
%   Input(s):
%   'l'             polynomial degree
%   'n'             derivative order n>=0
%   'P'             maxpol degree, controls cutoff frequency
%                       (a) n <= P < 2l   for centralized
%                       (b) n <= P < 2l-1 for staggered
%   'N'             size of derivative matrix (NxN)
%   'nod'           node type: 'staggered' or 'cenralized'
%   'sparse_flag'   construct matrix in sparse mode
%   'sym_flag'      numerical calculation: 'false', symbolic calculation: 'true'
%
%   Output(s):
%   'D'             backward derivative matrix N-by-N
%   'D_forward'     forward  derivative matrix N-by-N
%
%   See also DERIVCENT, DERIVSTAG, DERIVDIREC
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

switch nod
    case 'staggered'
        if sparse_flag
            zero_stacks = sparse(1, N-2*l);
        else
            zero_stacks = zeros(1, N-2*l);
        end
        [c] = derivstag_SavitzkyGolay(l, P, 0, n, sym_flag);
        row = [c, zero_stacks];
        col = [row(1) fliplr(row(end-length(row)+2:end))];
        block_middle = toeplitz(col,row);
        block_middle = block_middle(1:end-(2*l-1), :);
        %
        block_left = [];
        for s = 1: l
            [c] = derivstag_SavitzkyGolay(l, P, s, n, sym_flag);
            block_left = [[c, zero_stacks]; block_left];
        end
        %
        block_right = [];
        for s = 1: l-1
            [c] = derivstag_SavitzkyGolay(l, P, -s, n, sym_flag);
            block_right = [block_right; [zero_stacks, c]];
        end
        %%
    case 'centralized'
        if sparse_flag
            zero_stacks = sparse(1, N-(2*l+1));
        else
            zero_stacks = zeros(1, N-(2*l+1));
        end
        [c] = derivcent_SavitzkyGolay(l, P, 0, n, sym_flag);
        R = [c, zero_stacks];
        C = [R(1) fliplr(R(end-length(R)+2:end))];
        block_middle = toeplitz(C,R);
        block_middle = block_middle(1:end-2*l, :);
        %
        block_left = [];
        for s = 1: l
            [c] = derivcent_SavitzkyGolay(l, P, s, n, sym_flag);
            block_left = [[c, zero_stacks]; block_left];
        end
        %
        block_right = [];
        for s = 1: l
            [c] = derivcent_SavitzkyGolay(l, P, -s, n, sym_flag);
            block_right = [block_right; [zero_stacks, c]];
        end
end

%  stitch upper, middle, and lower blocks
D = [block_left; block_middle; block_right];

%   create forward matrix
D_forward = -flip(flip(D, 1), 2);