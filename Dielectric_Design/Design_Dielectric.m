function [params, ext, E] = Design_Dielectric(x, type, box, omega, E_0, ext_t, params_0, step_size, error_tol, outfile, design_flag, bounds, iter_tol)

% Perform inverse design for the Drude dielectric parameters of a
% many-bodied configuration.

% INPUTS
% x = (N-by-3) particle positions
% type = (N-by-1) particle type; integers b/w [1, m] where m = number of types
% box = (1-by-3) box dimensions
% omega = (N_omega-by-1) frequency
% E_0 = (1-by-3) field polarization
% ext_t = (N_omega-by-1) target extinction spectrum
% params_0 = (m-by-4) initial values of parameters [a, omega_p, gamma, eps_inf] 
% step_size = (scalar or m-by-4) optimization step size; can be different in each direction if desired
% error_tol = (scalar) error tolerance
% outfile = (string) output file name
% design_flag = (m-by-4) 1 if parameter is a design variable, 0 is parameter is held constant
% bounds = (m-by-4-2) lower and upper bounds of each parameter 
% 
% OUTPUTS
% params = (m-by-4-by-M) design parameters at each iteration
% ext = (M-by-N_omega) extinction spectrum at each iteration
% E = (M-by-1) error at each iteration

function eps_p = eps_fun(omega, eps_params)
    
    % INPUTS
    % omega = (N_omega-by-1)
    % eps_params = (N-by-3)
    %
    % OUTPUTS
    % eps_p = (N-by-N_omega)
    
    % Extract Drude parameters
    omega_p = eps_params(:,1);
    gamma = eps_params(:,2);
    eps_inf = eps_params(:,3);
    
    % Transpose from N_omega-by-1 to 1-by-N_omega
    omega = omega.';
    
    % Compute dielectric function
    eps_p = eps_inf - omega_p.^2./(omega.^2+1i*gamma.*omega);
    
end

% Helper function to enforce parameter bounds
function params = enforce_bound(params, min_bounds, max_bounds)
    
    % INPUTS
    % params = parameters
    % min_bounds = minimum bound for each parameter; same size as params
    % max_bounds = maximum bound for each parameter; same size as params
    %
    % OUTPUTS
    % params = bounded parameters; same size as input params
    
    params = max(params, min_bounds); % pin params at the minimum bound
    params = min(params, max_bounds); % pin params at the maximum bound

end

% Add the capacitance scripts to path
addpath('/Users/zacharysherman/Documents/Scattering/MATLAB/Capacitance_Sizes');

% Numerical parameters
xi = 0.5;

% Handle the initial "iteration" separately.
tic
M = 1; % iteration counter 
params_curr = params_0; % initial parameters
params = params_curr; % store initial parameters

% Radius and dielectric function of each particle
a = params_curr(type,1);
eps_p = eps_fun(omega, params_curr(type,2:4));

% Compute dipoles and extinction
[p_ave, p] = DipoleSpectrum(x, box, a, eps_p, E_0, xi);
ext_curr = 3/(4*pi)*omega.*imag(p_ave*E_0.');

% Initial error and gradient
[E, grad_E] = ComputeErrorGradient(omega, ext_curr, ext_t, p, params_curr, type);
ext = ext_curr.'; % store initial extinction

% Print status
fprintf('Iteration %.0f: error = %.4f, t = %.1f s\n', M, E, toc);

% Iterate until the error decreases below error_tol or we reach the max
% number of iterations
while (E(M) > error_tol) && (M < iter_tol)
    
    % Increment iteration counter
    M = M+1;
    
    % Update the parameters
    params_curr = params_curr - step_size.*grad_E.*design_flag; % new parameter values
    params_curr = enforce_bound(params_curr, bounds(:,:,1), bounds(:,:,2)); % enforce min and max bounds
    params = cat(3, params, params_curr); % store new parameter values
    
    % Radius and dielectric function of each particle
    a = params_curr(type,1);
    eps_p = eps_fun(omega, params_curr(type,2:4));
    
    % Compute new dipoles and extinction
    p_guess = p; % use previous dipoles as initial guess
    [p_ave, p] = DipoleSpectrum(x, box, a, eps_p, E_0, xi, p_guess);
    ext_curr = 3/(4*pi)*omega.*imag(p_ave*E_0.');
    ext = [ext; ext_curr.']; % store new extinction
    
    % Compute the error and gradient
    [E_curr, grad_E] = ComputeErrorGradient(omega, ext_curr, ext_t, p, params_curr, type);
    E = [E; E_curr];  % store new error
    
    % Print status
    fprintf('Iteration %.0f: error = %.4f, t = %.1f s\n', M, E_curr, toc);
    
    % Save to file
    save(outfile, 'omega', 'ext_t', 'params', 'ext', 'E', 'x', 'type', 'box')
    
end

end