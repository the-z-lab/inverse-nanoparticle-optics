function [] = Design_Position_TwoBody(x_0, box, omega, E_0, eps_p, a, x_E, step_size, T, error_tol, outfile)

    % Perform inverse design for two particles aligned with E_0 in the z
    % direction. Constrain the orientation and update only the center-to-center
    % vector.
    
    % INPUTS
    % x_0 = (N-by-3) initial particle positions
    % box = (1-by-3) box dimensions
    % omega = (N_omega-by-1) frequency
    % E_0 = (1-by-3) field polarization
    % eps_p = (N-by-N_omega) particle permittiity
    % a = (N-by-1) particle radii
    % x_E = (N_E-by-1) probe positions to compute field
    % step_size = (scalar) optimization step size
    % T = (scalar) simulated annealing temperature
    % error_tol = (scalar) error tolerance
    % outfile = (string) output file name
    %
    % OUTPUTS (to file)
    % omega = (scalar) frequency
    % x = (N-by-3-by-M) positions at each iteration
    % p = (N-by-3-by-M) extinction spectrum at each iteration
    % E = (N_E-by-3-by-M) field at each probe point at each iteration
    % E_ave = (M-by-1) average field intensity at each iteration
    
    % Add MPM scripts to path
    addpath('/Users/zacharysherman/Documents/Scattering/Inverse_Design/Adjoint/MPM')
    
    % Numerical parameters
    iter_tol = 500; % max number of iterations
    xi = 0.2;
    
    % Handle the initial "iteration" separately.
    tic
    M = 1; % iteration counter 
    x_curr = x_0; % initial positions
    
    % Initial FOMs
    [~, p_curr] = DipoleSpectrum(x_curr, box, a, eps_p, E_0, xi); % dipoles
    E_curr = FieldMap(x_E, x_curr, a, p_curr, box, xi); % scattered field at target positions
    E_ave_curr = mean(vecnorm(E_curr,2,2),1); % average field intensity
    
    % Initial gradient
    [Lambda_1, Lambda_2, grad_E, ~, ~, ~] = ComputeErrorGradient(x_curr, a, eps_p, p_curr, x_E, E_curr, box, xi);
    
    % Convert gradient w.r.t. x_1 and x_2 to gradient w.r.t. r
    grad_E_r = grad_E(2,3) - grad_E(1,3); % scalar
    
    % Store initial values
    x = x_curr;
    r = norm(x_curr(2,:) - x_curr(1,:));
    p = p_curr;
    E = E_curr;
    E_ave = E_ave_curr;
    
    % Print status
    fprintf('Iteration %.0f: E_ave = %.3f, grad_E = %.2g, t = %.1f s\n', M, E_ave_curr, grad_E_r, toc);
    
    % Iterate until the error or gradient decreases below error_tol or we reach the max
    % number of iterations
    while (norm(grad_E_r) > error_tol) && (M < iter_tol)
        
        % Increment iteration counter
        M = M+1;
        
        % Convert gradient w.r.t. x_1 and x_2 to gradient w.r.t. r
        grad_E_r = grad_E(2,3) - grad_E(1,3); % scalar
    
        % Update the parameters
        %delta_r = step_size.*grad_E_r./norm(grad_E_r); % normalized isn't effective in 1D
        delta_r = step_size.*grad_E_r;
        delta_x = [0, 0, -delta_r/2; 0, 0, delta_r/2];
        x_curr = x_curr + delta_x; % new parameter values
        
        % Add random noise
        x_curr = x_curr + T.*(rand(size(grad_E))-0.5);
        
        % Hard sphere repulsions
        x_curr = HardSphereDisplacement(x_curr,box);
        x_curr = HardSphereDisplacement_Probe(x_curr,x_E,box);
        
        % Shift-mod-shift particles back to primary box
        x_curr = mod(x_curr+box/2, box) - box/2;
        r_curr = norm(x_curr(2,:) - x_curr(1,:));
    
        % Compute new FOMs
        [~, p_curr] = DipoleSpectrum(x_curr, box, a, eps_p, E_0, xi, p_curr); % dipoles; use previous dipoles as initial guess
        E_curr = FieldMap(x_E, x_curr, a, p_curr, box, xi); % scattered field at target positions
        E_ave_curr = mean(vecnorm(E_curr,2,2),1); % average field intensity
    
        % Compute new gradient; use previous adjoints as initial guess
        [Lambda_1, Lambda_2, grad_E, ~, ~, ~] = ComputeErrorGradient(x_curr, a, eps_p, p_curr, x_E, E_curr, box, xi, Lambda_1, Lambda_2);
        
        % Store new values
        x = cat(3, x, x_curr);
        p = cat(3, p, p_curr);
        E = cat(3, E, E_curr);
        E_ave = [E_ave; E_ave_curr];
        r = [r; r_curr];
    
        % Print status
        fprintf('Iteration %.0f: E_ave = %.3f, grad_E = %.2g, t = %.1f s\n', M, E_ave_curr, grad_E_r, toc);
        
        % Save to file
        save(outfile, 'omega', 'x', 'r', 'p', 'E', 'E_ave', 'box')
        
    end
    
    end
    