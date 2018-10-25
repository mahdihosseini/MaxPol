function [gaussKernel] = gaussian_derivatives(l, sgm)

x_c = -l: l;
x_s = (-l+1/2): (l-1/2);

c = 1/(sqrt(2*pi)*sgm);
gaussKernel_c = c * exp(-x_c.^2/2/sgm^2);
gaussKernel_c = gaussKernel_c/sum(gaussKernel_c);

gaussKernel_s = c * exp(-x_s.^2/2/sgm^2);
gaussKernel_s = gaussKernel_s/sum(gaussKernel_s);

%%  lowpass filter
gaussKernel{1} = gaussKernel_c;

%%  first-OD
gaussKernel{2} = (-x_s/sgm^2) .* gaussKernel_s;

%%  second-OD
gaussKernel{3} = (x_c.^2-sgm^2)/sgm^4 .* gaussKernel_c;

%%  third-OD
gaussKernel{4} = (x_s.^3-3*x_s*sgm^2)/sgm^6 .* gaussKernel_s;

%%  fourth-OD
gaussKernel{5} = (x_c.^4 - 6*x_c.^2*sgm^2 + 3*sgm^4)/sgm^8 .* gaussKernel_c;