function [C] = tensordec(X, D)
%TENSORDEC Tensor decomposition
%   Decomposing m-dimensional volumetric data X (N_1-by-N_2-by-...-by-N_m)
%   by means of basis library D composed of m bases D{j} = D_j
%   corresponding to each dimension 'j'. For instance, D_j could be
%   derivative matrix defined by 'derivmtx' function from MaxPol package.
%
%   Inputs:
%   'X'     m-dimensional data N_1-by-N_2-by-...-by-N_m
%   'D'     Basis library composed of m dictionaries for docomposing
%           each dimension accordingly
%
%   Outputs:
%   'C'     Decomposed coefficients in nD volumetric shape, same as V.
%   Mathematically speaking, vec(C) = [B^T_n o ... o B^T_2 o B^T_1]*vec(V).
%   Here 'o' stands for Kronecker Product.
%
%
%   See also DERIVSTAG, DERIVMTX
%
%   Copyright   Mahdi S. Hosseini, BSc, MSc, MASc, PhD
%               University of Toronto, November 2016
%               The Edward S. Rogers Sr. Department of,
%               Electrical and Computer Engineering (ECE)
%               Toronto, ON, M5S3G4, Canada
%               Tel: +1 (416) 978 6845
%               email: mahdi.hosseini@mail.utoronto.ca

dim = size(X);
for iteration = 1: numel(dim)
    X_shifted = shiftdim(X, iteration - 1);
    size_shift = size(X_shifted);
    N_concatinate = prod(size_shift(2:end));
    X_shifted = reshape(X_shifted, [size_shift(1), N_concatinate]);
    C = D{iteration} * X_shifted;
    C = reshape(C, size_shift);
    X = shiftdim(C, numel(dim) - (iteration-1));
end
C = X;