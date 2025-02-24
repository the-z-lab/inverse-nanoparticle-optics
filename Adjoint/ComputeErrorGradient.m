function [Lambda_1, Lambda_2, grad_E, dEdp, grad_p, F_x] = ComputeErrorGradient(x_p, a, eps_p, p, x_E, E, box, xi, varargin)

    % Compute the error and the gradient in the error between the current and
    % target extinction spectra w.r.t. particle positions.
    
    % INPUTS
    % x_p = (N_p-by-3) particle positions
    % a = (N_p-by-1) particle radii
    % eps_p = (N_p-by-1) particle permittivities
    % p = (N_p-by-3-by-N_k) dipole moments; dims: particle, component, freq.
    % x_E = (N_E-by-3) field positions
    % E = (N_E-by-3) electric field
    % box = (1-by-3) box dimensions
    % xi = (scalar) Ewald parameter
    % Lambda_1_guess, Lambda_2_guess = (N-by-3) optional; initial guess for particle adjoints
    %
    % OUTPUTS
    % Lambda_1, Lambda_2 = (N-by-3) particle adjoints; return to use as initial guess at next iteration 
    % grad_E = (N-by-3) derivative of E w.r.t. positions
    
    % Array sizes
    N_p = size(x_p,1);
    N_E = size(x_E,1);
    
    % Perform the first "adjoint" calculation, in which field points/values are
    % sources and particle positions are probe points. First we need the
    % "force" between fields (source) and dipoles (probe). This is
    % dE^2/dx at constant p' and p''.
    F_x = 2*real(ForceSpectrum(x_E, x_p, E, conj(p), box, xi));
    
    % Perform the "adjoint" calculation. In which field points/values are
    % sources and particle positions are sample points
    dEdp = FieldMap(x_p, x_E, ones(N_E, 1), E, box, xi);
    dEdp = 2*dEdp/N_E; % if target is average field
    
    % Separate the derivatives w.r.t. p' and p''
    dEdp_re = real(dEdp);
    dEdp_im = imag(dEdp);
    
    % Get initial guess from input
    if isempty(varargin)
        Lambda_1_guess = zeros(N_p,3);
        Lambda_2_guess = zeros(N_p,3);
    else
        Lambda_1_guess = varargin{1};
        Lambda_2_guess = varargin{2};
    end
    
    % Solve for the adjoint variables
    [~, Lambda_1] = DipoleSpectrum(x_p, box, a, eps_p, dEdp_re, xi, Lambda_1_guess);
    [~, Lambda_2] = DipoleSpectrum(x_p, box, a, eps_p, dEdp_im, xi, Lambda_2_guess);
    
    % Compute the "force" spectrum for each particle
    F_re = 2*real(ForceSpectrum(x_p, x_p, p, Lambda_1, box, xi)); % dE^2/dp' (constant x and p'')
    F_im = 2*imag(ForceSpectrum(x_p, x_p, p, Lambda_2, box, xi)); % dE^2/dp'' (constant x and p')
    
    % Gradient of -E^2 w.r.t. particle positions
    grad_E = F_x + F_re + F_im;
    
    % For debugging
    grad_p_re = F_re(1,3)/dEdp_re(1,3);
    grad_p_im = F_im(1,3)/dEdp_im(1,3);
    grad_p = grad_p_re + 1i*grad_p_im;
    
    end