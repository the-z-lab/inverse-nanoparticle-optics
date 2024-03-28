function [x, ext, E] = Design_Position_test(x_0, box, omega, E_0, ext_t, eps_params, step_size, T1, T2, error_tol, outfile, dim)

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
addpath('/gscratch/zeelab/rdsanghavi/mpm/mpm/inversedesign/Meta_Gen/Capacitance_Hetero');

% Dielectric function
eps_p = eps_fun(omega, eps_params);

% Numerical parameters
iter_tol = 5; % max number of iterations
xi = 0.5;

% Handle the initial "iteration" separately.
tic
M = 1; % iteration counter 
x_curr = x_0; % initial parameters
x = x_curr; % store initial parameters

% Compute dipoles and extinction
[C, p] = CapacitanceTensorSpectrum(x, box, eps_p, xi);
Cavg = 1/2*(C(1,1,:) + C(3,3,:));
Cavg_sq = squeeze(Cavg);
ext_curr = 3/(4*pi)*omega.*imag(Cavg_sq);

% x and z components of p
px = squeeze(p(:,:,1,:)); 
pz = squeeze(p(:,:,3,:));
p = [px, pz];  % concatenate px and pz along row, N-by-6-Nk, see notes  

% Initial error and gradient
[E, grad_E] = ComputeErrorGradient(omega, ext_curr, ext_t, x_curr, p, box, xi);
ext = ext_curr.'; % store initial extinction

% p is N-by-6-Nk, reshape to N-3-2-Nk
p = permute(p, [1,2,4,3]);
p = cat(3, p(:,1:3,:,:),  p(:,4:6,:,:));

% Print status
fprintf('Iteration %.0f: error = %.4f, t = %.1f s\n', M, E, toc);

% Iterate until the error decreases below error_tol or we reach the max
% number of iterations
while (E(M) > error_tol) && (M < iter_tol)
    
    % Increment iteration counter
    M = M+1;
    
    % Update the parameters
    delta_x = -step_size.*grad_E;
    x_curr = x_curr + delta_x; % new parameter values
   
    if mod(M, 2) == 0;
        T = T1;
    else
        T = T2;

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
   
    [C, p] = Capacitance2DTensorSpectrum(x_curr, box, eps_p, dim, xi, p_guess);
    
    Cavg = 1/2*(C(1,1,:) + C(2,3,:)); %Cavg = 1/2*(C(1) + C(9)), Cxx is 1st element, Czz is 9th element
    Cavg_sq = squeeze(Cavg);    
    ext_curr = 3/(4*pi)*omega.*imag(Cavg_sq); 

    % x and z components of p
    px = p(:,:,1,:,:);
    pz = p(:,:,2,:,:);
    p_err = [px, pz];
    
    ext = [ext; ext_curr.']; % store new extinction
    
    % Compute the error and gradient
    [E_curr, grad_E] = ComputeErrorGradient(omega, ext_curr, ext_t, x_curr, p_err, box, xi);
    E = [E; E_curr];  % store new error
    
    % Print status
    fprintf('Iteration %.0f: error = %.4f, t = %.1f s\n', M, E_curr, toc);
    
    % Save to file
    save(outfile, 'omega', 'ext_t', 'x', 'ext', 'E', 'box')
    
end




end

%{ 
for second T annealing

iter_tol_2 = 10;

while (E(M) > error_tol) && (M < iter_tol_2)

    % Increment iteration counter
    M = M+1;

    % Update the parameters
    delta_x = -step_size.*grad_E;
    x_curr = x_curr + delta_x; % new parameter values

    % Add random noise
    x_curr = x_curr + T2.*(rand(size(grad_E))-0.5);

    % Hard sphere repulsions
    x_curr = HardSphereDisplacement(x_curr,box);

    % Shift-mod-shift particles back to primary box
    x_curr = mod(x_curr+box/2, box) - box/2;

    % Store new positions
    x = cat(3, x, x_curr);

    % Compute new dipoles and extinction
    p_guess = p_err; % use previous dipoles as initial guess

    [C, p] = Capacitance2DTensorSpectrum(x_curr, box, eps_p, dim, xi, p_guess);

    Cavg = 1/2*(C(1,1,:) + C(2,3,:)); %Cavg = 1/2*(C(1) + C(9)), Cxx is 1st element, Czz is 9th element
    Cavg_sq = squeeze(Cavg);
    ext_curr = 3/(4*pi)*omega.*imag(Cavg_sq);

    % x and z components of p
    px = p(:,:,1,:,:);
    pz = p(:,:,2,:,:);
    p_err = [px, pz];

    ext = [ext; ext_curr.']; % store new extinction

    % Compute the error and gradient
    [E_curr, grad_E] = ComputeErrorGradient(omega, ext_curr, ext_t, x_curr, p_err, box, xi);
    E = [E; E_curr];  % store new error

    % Print status
    fprintf('Iteration %.0f: error = %.4f, t = %.1f s\n', M, E_curr, toc);

    % Save to file
    save(outfile, 'omega', 'ext_t', 'x', 'ext', 'E', 'box')

end

%}

