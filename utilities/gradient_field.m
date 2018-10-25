function [grad] = gradient_field(i, params)


%%  gradient
grad{1} = imfilter(i, params.d1, 'conv', 'symmetric');
grad{2} = imfilter(i, params.d1', 'conv', 'symmetric');
