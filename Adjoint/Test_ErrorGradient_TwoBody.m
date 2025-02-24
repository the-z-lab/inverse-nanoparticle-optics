function [] = Test_ErrorGradient_TwoBody()

    % Initialize configuration
    a = [1; 1];
    box = [50, 50, 50];
    
    % Dielectric parameters
    omega_p = [1; 1];
    gamma = [0.05; 0.05];
    eps_inf = [1; 1];
    
    % Optical parameters
    omega = 0.55; % incident frequency
    E_0 = [0, 0, 1]; % incident polarization
    x_E = [0, 0, 0]; % probe points
    
    % Particle permittivity from the Drude model
    eps_p = eps_inf - omega_p.^2./(omega.^2 + 1i*gamma.*omega);
    
    % Add MPM scripts to path
    addpath('/Users/zacharysherman/Documents/Scattering/Inverse_Design/Adjoint/MPM')
    
    
    %% Brute force parameter sweep
    
    % Initialize
    r = (5:-0.02:2)'; % separation values
    %r = 2.1;
    N_r = length(r);
    E_mag = zeros(N_r,1);
    p_mag = zeros(N_r,1);
    grad_E2_z = zeros(N_r,1);
    dEdp_z = zeros(N_r,1);
    grad_p_zz = zeros(N_r,1);
    Lambda_1_z = zeros(N_r,1);
    Lambda_2_z = zeros(N_r,1);
    partial_E2_z = zeros(N_r,1);
    E_an_mag = zeros(N_r,1);
    p_an_mag = zeros(N_r,1);
    grad_E2_an_z = zeros(N_r,1);
    dEdp_an_z = zeros(N_r,1);
    grad_p_an_zz = zeros(N_r,1);
    Lambda_1_an_z = zeros(N_r,1);
    Lambda_2_an_z = zeros(N_r,1);
    partial_E2_an_z = zeros(N_r,1);
    error = zeros(N_r,1);
    xi = 0.2;
    
    % Initial guess
    alpha = (eps_p-1)./(eps_p+2);
    alpha(isinf(eps_p)) = 1;
    p = 4*pi*a.^3.*alpha.*E_0; 
    
    % Loop through separations
    for i = 1:N_r
        
        % Particle positions
        x = [0, 0, -r(i)/2; 0, 0, r(i)/2];
        
        % Compute FOMs
        [~, p] = DipoleSpectrum(x, box, a, eps_p, E_0, xi, p); % dipoles; use previous dipoles as initial guess
        E = FieldMap(x_E, x, a, p, box, xi); % scattered field at target positions
        E_mag(i) = mean(vecnorm(E,2,2),1); % average field intensity
        p_mag(i) = mean(vecnorm(p,2,2),1); % average field intensity
        
        % Compute gradient in |E|^2
        [Lambda_1, Lambda_2, grad_E2, dEdp, grad_p_zz(i), partial_E2] = ComputeErrorGradient(x, a, eps_p, p, x_E, E, box, xi);
        grad_E2_z(i) = -grad_E2(1,3);
        dEdp_z(i) = dEdp(1,3);
        %grad_p_zz(i) = grad_E2_z(i)/dEdp_z(i);
        Lambda_1_z(i) = Lambda_1(1,3);
        Lambda_2_z(i) = Lambda_2(1,3);
        partial_E2_z(i) = partial_E2(1,3);
        
        % Analytic solution
        [p_an, E_an, grad_E2_an, dEdp_an, grad_p_an, Lambda_1_an, Lambda_2_an, partial_E2_an] = AnalyticTwoBody(x, alpha(1), box, E_0, x_E);
        
        % Store solution
        p_an_mag(i) = norm(p_an(1,:));
        E_an_mag(i) = norm(E_an);
        grad_E2_an_z(i) = grad_E2_an(3);
        dEdp_an_z(i) = dEdp_an(1,3);
        grad_p_an_zz(i) = grad_p_an(3,3);
        Lambda_1_an_z(i) = Lambda_1_an(1,3);
        Lambda_2_an_z(i) = Lambda_2_an(1,3);
        partial_E2_an_z(i) = partial_E2_an(3);
        
        error(i) = norm(p_an(1,:) - p(1,:))/norm(p_an(1,:));
        
    end
    
    % Numeric gradient
    grad_E2_num = gradient(E_an_mag.^2,r(2)-r(1));
    
    % Equivalent way to compute grad E^2
    %grad_E2_prime = -real(dEdp_z).*real(grad_p_zz) - imag(dEdp_z).*imag(grad_p_zz);
    grad_E2_an_prime = -2*real(dEdp_an_z).*real(grad_p_an_zz) - 2*imag(dEdp_an_z).*imag(grad_p_an_zz) - partial_E2_an_z;
    %grad_E2_prime = -real(dEdp_z.*grad_p_zz) - imag(dEdp_z.*grad_p_zz);
    
    % Plot dipole
    close all
    % figure; hold on
    % plot(r, p_mag, 'r', 'LineWidth', 3)
    % plot(r, p_an_mag, '--k', 'LineWidth', 3)
    % 
    % % Format axes
    % xlabel('r')
    % ylabel('|p|')
    % set(gca,'FontUnits','normalized','FontSize',0.05,'FontWeight','bold','LineWidth',1,'PlotBoxAspectRatio',[1,1,1])
    % 
    % % Plot field
    % figure; hold on
    % plot(r, E_mag, 'LineWidth', 3)
    % plot(r, E_an_mag, '--k', 'LineWidth', 3)
    % 
    % % Format axes
    % xlabel('r')
    % ylabel('|E|')
    % set(gca,'FontUnits','normalized','FontSize',0.05,'FontWeight','bold','LineWidth',1,'PlotBoxAspectRatio',[1,1,1])
    
    % Plot gradient in |E|^2
    figure; hold on
    plot(r, grad_E2_z, 'b', 'LineWidth', 3)
    %plot(r, 0.5*grad_E2_prime, ':c', 'LineWidth', 3)
    plot(r, grad_E2_an_prime, 'g', 'LineWidth', 3)
    plot(r, -grad_E2_an_z, '--r', 'LineWidth', 3)
    plot(r, grad_E2_num, ':k', 'LineWidth', 3)
    
    % Format axes
    xlabel('r')
    ylabel('\nabla |E|^2')
    set(gca,'FontUnits','normalized','FontSize',0.05,'FontWeight','bold','LineWidth',1,'PlotBoxAspectRatio',[1,1,1])
    
    % % Plot dEdp
    % figure; hold on
    % plot(r, real(dEdp_z), 'b', 'LineWidth', 3)
    % plot(r, imag(dEdp_z), 'k', 'LineWidth', 3)
    % plot(r, real(dEdp_an_z), '--r', 'LineWidth', 3)
    % plot(r, imag(dEdp_an_z), ':r', 'LineWidth', 3)
    % 
    % % Format axes
    % xlabel('r')
    % ylabel('d|E|^2/dp')
    % set(gca,'FontUnits','normalized','FontSize',0.05,'FontWeight','bold','LineWidth',1,'PlotBoxAspectRatio',[1,1,1])
    
    % % Plot grad p 
    % figure; hold on
    % plot(r, 0.5*real(grad_p_zz), 'b', 'LineWidth', 3)
    % plot(r, 0.5*imag(grad_p_zz), 'k', 'LineWidth', 3)
    % plot(r, real(grad_p_an_zz), '--r', 'LineWidth', 3)
    % plot(r, imag(grad_p_an_zz), ':r', 'LineWidth', 3)
    % 
    % % Format axes
    % xlabel('r')
    % ylabel('\nabla p')
    % set(gca,'FontUnits','normalized','FontSize',0.05,'FontWeight','bold','LineWidth',1,'PlotBoxAspectRatio',[1,1,1])
    % 
    % % Plot Lambda_1 
    % figure; hold on
    % plot(r, real(Lambda_1_z), 'b', 'LineWidth', 3)
    % plot(r, imag(Lambda_1_z), 'k', 'LineWidth', 3)
    % plot(r, real(Lambda_1_an_z), '--r', 'LineWidth', 3)
    % plot(r, imag(Lambda_1_an_z), ':r', 'LineWidth', 3)
    % 
    % % Format axes
    % xlabel('r')
    % ylabel('\Lambda_1')
    % set(gca,'FontUnits','normalized','FontSize',0.05,'FontWeight','bold','LineWidth',1,'PlotBoxAspectRatio',[1,1,1])
    % 
    % % Plot Lambda_2 
    % figure; hold on
    % plot(r, real(Lambda_2_z), 'b', 'LineWidth', 3)
    % plot(r, imag(Lambda_2_z), 'k', 'LineWidth', 3)
    % plot(r, real(Lambda_2_an_z), '--r', 'LineWidth', 3)
    % plot(r, imag(Lambda_2_an_z), ':r', 'LineWidth', 3)
    % 
    % % Format axes
    % xlabel('r')
    % ylabel('\Lambda_2')
    % set(gca,'FontUnits','normalized','FontSize',0.05,'FontWeight','bold','LineWidth',1,'PlotBoxAspectRatio',[1,1,1])
    
    % Plot partial dE^2/dx at constant p
    figure; hold on
    plot(r, partial_E2_z, 'b', 'LineWidth', 3)
    plot(r, partial_E2_an_z, '--k', 'LineWidth', 3)
    
    % Format axes
    xlabel('r')
    ylabel('\partial |E|^2/\partial r')
    set(gca,'FontUnits','normalized','FontSize',0.05,'FontWeight','bold','LineWidth',1,'PlotBoxAspectRatio',[1,1,1])
    end