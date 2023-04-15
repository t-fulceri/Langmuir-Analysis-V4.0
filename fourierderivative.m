function dfdx = fourierderivative(f,a,b)
% FOURIERDERIVATIVE Fourier derivative
%     dfdx = FOURIERDERIVATIVE(f,a,b) approximates the derivative a
%     discrete function f over the domain (a,b).  f, a vector, must be
%     uniformly sampled, periodic, and contain an even number of samples.
%     For best results, f should be periodic such that f(x + a) = f(x + b).
%     As an example,
%
%          x = linspace(0,pi);
%          f = exp(cos(x).*sin(2*x));
%          dfdx = fourierderivative(f,0,pi);
%
%     Results for nonperiodic f are dubious.

Nx = max(size(f));
k = 2*pi/(b-a)*[0:Nx/2-1 0 -Nx/2+1:-1];

dfdx = ifft(1i*k.*fft(f));

end