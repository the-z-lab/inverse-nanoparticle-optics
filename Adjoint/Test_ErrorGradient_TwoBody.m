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
E_an_mag = zeros(N_r,1);
p_an_mag = zeros(N_r,1);
grad_E2_an_z = zeros(N_r,1);
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
    [~, ~, grad_E2] = ComputeErrorGradient(x, a, eps_p, p, x_E, E, box, xi);
    grad_E2_z(i) = -grad_E2(1,3);
    
    % Analytic solution
    [p_an, E_an, grad_E2_an] = AnalyticTwoBody(x, alpha(1), box, E_0, x_E);
    
    % Store solution
    p_an_mag(i) = norm(p_an(1,:));
    E_an_mag(i) = norm(E_an);
    grad_E2_an_z(i) = grad_E2_an(3);
    
    error(i) = norm(p_an(1,:) - p(1,:))/norm(p_an(1,:));
    
end

% Numeric gradient
grad_E2_num = gradient(E_an_mag.^2,r(2)-r(1));

% Plot dipole
close all
figure; hold on
plot(r, p_mag, 'r', 'LineWidth', 3)
plot(r, p_an_mag, '--k', 'LineWidth', 3)

% Format axes
xlabel('r')
ylabel('|p|')
set(gca,'FontUnits','normalized','FontSize',0.05,'FontWeight','bold','LineWidth',1,'PlotBoxAspectRatio',[1,1,1])

% Plot field
figure; hold on
plot(r, E_mag, 'LineWidth', 3)
plot(r, E_an_mag, '--k', 'LineWidth', 3)

% Format axes
xlabel('r')
ylabel('|E|')
set(gca,'FontUnits','normalized','FontSize',0.05,'FontWeight','bold','LineWidth',1,'PlotBoxAspectRatio',[1,1,1])

% Plot gradient in |E|^2
figure; hold on
plot(r, grad_E2_z, 'b', 'LineWidth', 3)
plot(r, -grad_E2_an_z, '--r', 'LineWidth', 3)
plot(r, grad_E2_num, ':k', 'LineWidth', 3)


% Format axes
xlabel('r')
ylabel('\nabla |E|^2')
set(gca,'FontUnits','normalized','FontSize',0.05,'FontWeight','bold','LineWidth',1,'PlotBoxAspectRatio',[1,1,1])


end