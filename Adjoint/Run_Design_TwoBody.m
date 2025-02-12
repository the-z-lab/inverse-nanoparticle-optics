function [] = Run_Design_TwoBody()

% Initialize configuration
x_0 = [0, 0, -1.5; 0, 0, 1.5];
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
    
% Initialize optimization parameters
step_size = 0.05;
T = 0;
error_tol = 10^-8;
outfile = 'twobody.mat';

% Add MPM scripts to path
addpath('/Users/zacharysherman/Documents/Scattering/Inverse_Design/Adjoint/MPM')

% Run the design
%Design_Position(x_0, box, omega, E_0, eps_p, a, x_E, step_size, T, error_tol, outfile)

%% Brute force parameter sweep

% Initialize
r = (10:-0.1:2)'; % separation values
N_r = length(r);
E_mag = zeros(N_r,1);
p_mag = zeros(N_r,1);
E_an_mag = zeros(N_r,1);
p_an_mag = zeros(N_r,1);
grad_E2_an = zeros(N_r,1);
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
    
    % Probe positions
    x_E = [0, 0, 0];
    
    % Compute FOMs
    [~, p] = DipoleSpectrum(x, box, a, eps_p, E_0, xi, p); % dipoles; use previous dipoles as initial guess
    E = FieldMap(x_E, x, a, p, box, xi); % scattered field at target positions
    E_mag(i) = mean(vecnorm(E,2,2),1); % average field intensity
    p_mag(i) = mean(vecnorm(p,2,2),1); % average field intensity
    
    % Analytic solution
    [p_an, E_an, grad_E2_an(i,:)] = AnalyticTwoBody(x, alpha(1), box, E_0, x_E);
    
    % Store solution
    p_an_mag(i) = norm(p_an(1,:));
    E_an_mag(i) = norm(E_an);
   
    error(i) = norm(p_an(1,:) - p(1,:))/norm(p_an(1,:));
    
end

% Analytic dipole
%p_an = 4*pi*alpha(1)./(1-2*alpha(1)./r.^3);
%p_an_ave = abs(p_an);

% Plot brute force sweep
close all
figure; hold on
plot(r, E_mag, 'b', 'LineWidth', 2)
plot(r, p_mag, 'r', 'LineWidth', 2)
plot(r, p_an_mag, '--k', 'LineWidth', 2)

% Format axes
xlabel('r')
ylabel('|E|')
set(gca,'FontUnits','normalized','FontSize',0.05,'FontWeight','bold','LineWidth',1,'PlotBoxAspectRatio',[1,1,1])

% Plot error
figure; hold on
plot(r, error, 'k', 'LineWidth', 2)

% Format axes
xlabel('r')
ylabel('|p - p_an|')
set(gca,'FontUnits','normalized','FontSize',0.05,'FontWeight','bold','LineWidth',1,'PlotBoxAspectRatio',[1,1,1])


end