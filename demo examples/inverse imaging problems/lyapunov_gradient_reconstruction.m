function [X] = lyapunov_gradient_reconstruction(grad, D)

[M, N, Q] = size(grad{1});
Dx = D{1};
Dy = D{2};

%%  Lyapanov-Equation Solver (Buld-in MATLAB function)
A = (Dy'*Dy);
B = (Dx'*Dx);
%
for iteration = 1: Q
    C = grad{1}(:,:,iteration) * Dx + Dy' * grad{2}(:,:,iteration);
    X(:,  :,  iteration) = lyap(A, B, -C);
end
X = real(X);
