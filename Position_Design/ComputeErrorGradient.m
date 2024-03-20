function [E, grad_E] = ComputeErrorGradient(omega, ext, ext_t, x, p, box, xi)

% Compute the error and the gradient in the error between the current and
% target extinction spectra w.r.t. particle positions.

% INPUTS
% omega = (N_omega-by-1) frequency
% ext = (N_omega-by-1) extinction spectrum
% ext_t = (N_omega-by-1) target extinction spectrum
% x = (N-by-3) particle positions
% p = (N-by-3-by-N_k) dipole moments; dims: particle, component, freq.
% box = (1-by-3) box dimensions
% xi = (scalar) Ewald parameter
%
% OUTPUTS
% E = (scalar) error
% grad_E = (N-by-3) derivative of E w.r.t. positions

% Put frequency in the third dimension, i.e. reshape from N_omega-by-1 to 1-by-1-by-N_omega
omega = permute(omega, [3, 2, 1]);
ext = permute(ext, [3, 2, 1]);
ext_t = permute(ext_t, [3, 2, 1]);

% Compute the "force" spectrum for each particle
F = ForceSpectrum(x, p, box, xi);

% Gradient of extinction spectrum
grad_ext = 3/(4*pi)*omega.*imag(2*F);

% Compute the error in the extinction spectrum
omega = permute(omega, [3, 2, 1]); % omega has to be a column vector for trapz
E = trapz(omega, (ext-ext_t).^2, 3)/trapz(omega, ext_t.^2, 3);

% Gradient of the error
grad_E = 2*trapz(omega, (ext-ext_t).*grad_ext, 3);

end