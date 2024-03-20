function [x, ext, E] = Design_Position(x_0, box, omega, E_0, ext_t, eps_params, step_size, T, error_tol, outfile)

% Perform inverse design for particle positions.

% INPUTS
% x_0 = (N-by-3) initial particle positions
% box = (1-by-3) box dimensions
% omega = (N_omega-by-1) frequency
% E_0 = (1-by-3) field polarization
% ext_t = (N_omega-by-1) target extinction spectrum
% eps_params = (N-by-3) Drude dielectric parameters
% step_size = (scalar) optimization step size
% T = (1-by-3) simulated annealing temperature
% error_tol = (scalar) error tolerance
% outfile = (string) output file name
%
% OUTPUTS
% x = (N-by-3-by-M) positions at each iteration
% ext = (M-by-N_omega) extinction spectrum at each iteration
% E = (M-by-1) error at each iteration

% Helper function for dielectric functions
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

% Add the capacitance scripts to path
addpath('/gscratch/zeelab/zsherm/Meta_Gen/Capacitance_Hetero');
%addpath('/Users/zacharysherman/Documents/Scattering/MATLAB/Capacitance_Hetero');

% Dielectric function
eps_p = eps_fun(omega, eps_params);

% Numerical parameters
iter_tol = 5000; % max number of iterations
xi = 0.5;

% Handle the initial "iteration" separately.
tic
M = 1; % iteration counter 
x_curr = x_0; % initial parameters
x = x_curr; % store initial parameters

% Compute dipoles and extinction
[C, p] = CapacitanceSpectrum(x, box, eps_p, xi);
ext_curr = 3/(4*pi)*omega.*imag(C);

% Initial error and gradient
[E, grad_E] = ComputeErrorGradient(omega, ext_curr, ext_t, x_curr, p, box, xi);
ext = ext_curr.'; % store initial extinction

% Print status
fprintf('Iteration %.0f: error = %.4f, t = %.1f s\n', M, E, toc);

% Iterate until the error decreases below error_tol or we reach the max
% number of iterations
while (E(M) > error_tol) && (M < iter_tol)
    
    % Increment iteration counter
    M = M+1;
    
    % Update the parameters
    %grad_E = grad_E/norm(grad_E(:)); % unit vector
    delta_x = -step_size.*grad_E;
    x_curr = x_curr + delta_x; % new parameter values
    
    % Add random noise
    x_curr = x_curr + T.*(rand(size(grad_E))-0.5);
    
    % Hard sphere repulsions
    x_curr = HardSphereDisplacement(x_curr,box);
    
    % Shift-mod-shift particles back to primary box
    x_curr = mod(x_curr+box/2, box) - box/2;
    
    % Store new positions
    x = cat(3, x, x_curr);
    
    % Compute new dipoles and extinction
    p_guess = p; % use previous dipoles as initial guess
    [C, p] = CapacitanceSpectrum(x_curr, box, eps_p, xi, p_guess);
    ext_curr = 3/(4*pi)*omega.*imag(C);
    ext = [ext; ext_curr.']; % store new extinction
    
    % Compute the error and gradient
    [E_curr, grad_E] = ComputeErrorGradient(omega, ext_curr, ext_t, x_curr, p, box, xi);
    E = [E; E_curr];  % store new error
    
    % Print status
    fprintf('Iteration %.0f: error = %.4f, t = %.1f s\n', M, E_curr, toc);
    
    % Save to file
    save(outfile, 'omega', 'ext_t', 'x', 'ext', 'E', 'box')
    
end

end
