function [f, grad, hess, third_order_tensor, fourth_order_tensor, params] = sinusoidal_gratting_2D(frq, delta_x)

%%  synthetic image generation
%   x-axis grid points
x_a = -1;
x_b = +1;
y_a = -1;
y_b = +1;
N_x = 128;
N_y = 128;
h_x = (x_b-x_a)/(N_x-1);
x = [x_a: h_x: x_b]+delta_x;
%   y-axis grid points
h_y = (y_b-y_a)/(N_y-1);
y = [y_a: h_y: y_b];
%   2D mesh grid (x,y)
[X_centralized, Y_centralized] = meshgrid(x, y);
[X_staggered, Y_staggered] = meshgrid(x-h_x/2, y-h_y/2);

%%  sinusoidal gratting syntethic data
frq_x = frq;
frq_y = frq;
% frq_x = 0.25;
% frq_y = 0.25;
A = (2*pi*frq_x*X_centralized).^2+(2*pi*frq_y*Y_centralized).^2;
f = sin(A);
range_f = range(f(:));
f = normal(f)*255;
multiplying_factor = 255/range_f;

%%  analytical gradients
grad_centerilzied{2} = multiplying_factor*2*(2*pi*frq_y)^2*Y_centralized.*cos(A);
grad_centerilzied{1} = multiplying_factor*2*(2*pi*frq_x)^2*X_centralized.*cos(A);
%
grad_staggered{2} = multiplying_factor*2*(2*pi*frq_y)^2*Y_staggered.*cos((2*pi*frq_x*X_centralized).^2+(2*pi*frq_y*Y_staggered).^2);
grad_staggered{1} = multiplying_factor*2*(2*pi*frq_x)^2*X_staggered.*cos((2*pi*frq_x*X_staggered).^2+(2*pi*frq_y*Y_centralized).^2);

grad.centralized = grad_centerilzied;
grad.staggered = grad_staggered;

%%  analytical hessians
hess_centerilzied{1} = multiplying_factor*(2*(2*pi*frq_x)^2*cos(A) - 4*(2*pi*frq_x)^4*X_centralized.^2.*sin(A));
hess_centerilzied{3} = multiplying_factor*(2*(2*pi*frq_y)^2*cos(A) - 4*(2*pi*frq_y)^4*Y_centralized.^2.*sin(A));
hess_centerilzied{2} = -multiplying_factor*4*(2*pi*frq_x)^2*(2*pi*frq_y)^2*X_centralized.*Y_centralized.*sin(A);
%
hess_staggered{1} = multiplying_factor*(2*(2*pi*frq_x)^2*cos((2*pi*frq_x*X_staggered).^2+(2*pi*frq_y*Y_centralized).^2) - 4*(2*pi*frq_x)^4*X_staggered.^2.*sin((2*pi*frq_x*X_staggered).^2+(2*pi*frq_y*Y_centralized).^2));
hess_staggered{3} = multiplying_factor*(2*(2*pi*frq_y)^2*cos((2*pi*frq_x*X_centralized).^2+(2*pi*frq_y*Y_staggered).^2) - 4*(2*pi*frq_x)^4*Y_staggered.^2.*sin((2*pi*frq_x*X_centralized).^2+(2*pi*frq_y*Y_staggered).^2));
hess_staggered{2} = -multiplying_factor*4*(2*pi*frq_x)^2*(2*pi*frq_y)^2*X_staggered.*Y_staggered.*sin((2*pi*frq_x*X_staggered).^2+(2*pi*frq_y*Y_staggered).^2);
hess.centralized = hess_centerilzied;
hess.staggered = hess_staggered;

%%  analytical third order tensor
w_x = 2*pi*frq_x;
w_y = 2*pi*frq_y;
tensor_staggered_third{1} = multiplying_factor*(-12*w_x^4*X_staggered.*sin((w_x*X_staggered).^2+(w_y*Y_centralized).^2)-8*w_x^6*X_staggered.^3.*cos((w_x*X_staggered).^2+(w_y*Y_centralized).^2));
tensor_centralized_third{1} = multiplying_factor*(-12*w_x^4*X_centralized.*sin(A)-8*w_x^6*X_centralized.^3.*cos(A));

third_order_tensor.staggered = tensor_staggered_third;
third_order_tensor.centralized = tensor_centralized_third;

%%  analytical third order tensor
w_x = 2*pi*frq_x;
w_y = 2*pi*frq_y;
tensor_staggered_fourth{1} = multiplying_factor*((16*w_x^8*X_staggered.^4-12*w_x^4).*sin((w_x*X_staggered).^2+(w_y*Y_centralized).^2)-48*w_x^6*X_staggered.^2.*cos((w_x*X_staggered).^2+(w_y*Y_centralized).^2));
tensor_centralized_fourth{1} = multiplying_factor*((16*w_x^8*X_centralized.^4-12*w_x^4).*sin(A)-48*w_x^6*X_centralized.^2.*cos(A));

fourth_order_tensor.staggered = tensor_staggered_fourth;
fourth_order_tensor.centralized = tensor_centralized_fourth;

%%
params.N_x = N_x;
params.N_y = N_y;
params.h_x = h_x;
params.h_y = h_y;
params.X_centralized = X_centralized;
params.Y_centralized = Y_centralized;
params.X_staggered = X_staggered;
params.Y_staggered = Y_staggered;