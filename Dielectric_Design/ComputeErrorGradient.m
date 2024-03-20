function [E, grad_E] = ComputeErrorGradient(omega, ext, ext_t, p, params, type)

% Compute the error and the gradient in the error between the current and
% target extinction spectra.

% INPUTS
% omega = (N_omega-by-1) frequency
% ext = (N_omega-by-1) extinction spectrum
% ext_t = (N_omega-by-1) target extinction spectrum
% p = (N-by-3-by-N_k) dipole moments; dims: particle, component, freq.
% params = (m-by-4) Drude dielectric parameters of each particle; columns are [a, omega_p, gamma, eps_inf]
% type = (N-by-1) particle type; integers b/w [1, m] where m = number of types
%
% OUTPUTS
% E = (scalar) error
% grad_E = (m-by-4) derivative of E w.r.t. dielectric parameters

% Array sizes
N = size(p,1); % number of particles
N_omega = size(omega,1); % number of frequencies
m = max(type); % number of types 

% Extract Drude params
a = params(:,1);
omega_p = params(:,2);
gamma = params(:,3);
eps_inf = params(:,4);

% Put frequency in the third dimension, i.e. reshape from N_omega-by-1 to 1-by-1-by-N_omega
omega = permute(omega, [3, 2, 1]);
ext = permute(ext, [3, 2, 1]);
ext_t = permute(ext_t, [3, 2, 1]);

% Partial derivatives of eps_p w.r.t. Drude params
d_omega_p = -2*omega_p./(omega.^2+1i*gamma.*omega); % omega_p derivatie
d_gamma = 1i*omega_p.^2./(omega.*(omega+1i*gamma).^2); % gamma derivative
d_eps_inf = ones(m, 1, N_omega); % eps_inf derivative
grad_eps = [d_omega_p, d_gamma, d_eps_inf]; % combine into N-by-3-by-N_omega array

% Gradient of potential tensor w.r.t. Drude params
eps_p = eps_inf - omega_p.^2./(omega.^2+1i*gamma.*omega);
alpha = (eps_p - 1)./(eps_p + 2);
grad_M_lambda = -3./(4*pi*a.^3.*(eps_p-1).^2).*grad_eps;

% Gradient of potential tensor w.r.t. particle radius
grad_M_a = -3./(4*pi*a.^4.*alpha);

% Combine radius and Drude params
grad_M = [grad_M_a, grad_M_lambda];

% Loop through types and accumulate the dipole dot product 
sum_p_dot_p = zeros(m,1,N_omega); % intialize output array
for i = 1:m
    type_i = type == i;
    sum_p_dot_p(i,1,:) = sum(sum(p(type_i,:,:).^2,2),1);  % don't use dot(p,p,2), since dot complex conjugates the first argument 
end

% Gradient of extinction spectrum
grad_ext = -3/(4*pi)*omega.*imag(grad_M.*sum_p_dot_p);

% Compute the error in the extinction spectrum
omega = permute(omega, [3, 2, 1]); % have to permute back to column vector for trapz
E = trapz(omega, (ext-ext_t).^2, 3)/trapz(omega, ext_t.^2, 3);

% Gradient of the error
grad_E = 2*trapz(omega, (ext-ext_t).*grad_ext, 3);

end