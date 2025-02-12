function [] = Tests()

% Add files to path
addpath('MPM')

% Initialize configuration
x_0 = [0, 0, -1.1; 0, 0, 1.1]+2*pi;
%x_0 = [-2, 0, 0; 2, 0, 0];
a = [1; 1];
box = [50, 50, 50];

% Dielectric parameters
omega_p = [1; 1];
gamma = [0.05; 0.05];
eps_inf = [1; 1];

% Optical parameters
omega = 0.55; % incident frequency
E_0 = [0, 0, 1/(4*pi)]; % incident polarization
x_E = [0, 0, 0]; % probe points

% Particle permittivity from the Drude model
eps_p = eps_inf - omega_p.^2./(omega.^2 + 1i*gamma.*omega);
alpha = (eps_p-1)./(eps_p+2);
alpha(isinf(eps_p)) = 1;

% Analytic dipole
r = x_0(1,:) - x_0(2,:);
r_mag = norm(r);
r_hat = r/r_mag;
p_z = 4*pi*alpha./(1-2*alpha./r_mag.^3);
p = [zeros(2,2), p_z];
%p = [1, 2, 3; 2.5, 1.5, 0.3];
%p = 8*pi*[0, 0, 1; 0, 0, 1];

% Analytic force
p_i = p(1,:);
p_j = p(2,:);
F_i = 3/(4*pi*r_mag^4)*( (p_i*p_j.')*r_hat - 5*(p_i*r_hat')*(p_j*r_hat')*r_hat + (p_i*r_hat')*p_j + (p_j*r_hat')*p_i );
F_an = [F_i; -F_i]

xi = 0.1;
F = ForceSpectrum(x_0, p, p, box, xi);
F(1,:)
norm(F(1,:) - F_an(1,:))

xi = 0.2;
F = ForceSpectrum(x_0, p, p, box, xi);
F(1,:)
norm(F(1,:) - F_an(1,:))

xi = 0.5;
F = ForceSpectrum(x_0, p, p, box, xi);
F(1,:)
norm(F(1,:) - F_an(1,:))

xi = 1.0;
F = ForceSpectrum(x_0, p, p, box, xi);
F(1,:)
norm(F(1,:) - F_an(1,:))

end