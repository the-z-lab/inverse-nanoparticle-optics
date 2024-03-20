function [] = RunDesign_PhaseSeparation()

% Add needed directories to the file path
addpath('/gscratch/zeelab/zsherm/Meta_Gen/Lattices')
addpath('/gscratch/zeelab/zsherm/Meta_Gen/Capacitance_Hetero')

%% Target extinction spectrum

% Construct a lattice
N_x = 9; % number of unit cells in x dimension
N_y = 5; % number of unit cells in y dimension
phi = pi/(2*sqrt(3)); % area fraction
[x_t, L] = HexLattice(N_x, N_y, phi); % place particles on a lattice
N = size(x_t, 1); % number of particles

% Put particles in the xz plane
x_t = [x_t(:,1), zeros(N,1)+10^-6, x_t(:,2)];
L_y = 20; % box length in the out-of-plane (y) dimenison
box = [L(1), L_y, L(2)]; % box dimensions

% Choose dielectric parameters
omega = (0.01:0.01:1).'; % frequency
gamma = 0.1; % damping 
eps_inf = 1; % high-freq. permittivity

% Particle permittivity
eps_p = eps_inf - 1./(omega.^2+1i*gamma.*omega);

% Assign all particles the same dielectric function
eps_p = repmat(eps_p.', N, 1); % dims: particle, freq.
eps_params = [1, gamma, eps_inf].*ones(N,1); % format dielectric parameters for later

E_0 = [0, 0, 1]; % field polarization
xi = 0.50; % Ewald parameter

% Compute the target extinction spectrum
[C, ~] = CapacitanceSpectrum(x_t, box, eps_p, xi);
ext_t = 3/(4*pi)*omega.*imag(C); % target extinction

%% Optimization initialization

% Optimization geometry
N = N; % number of particles
phi = 0.40; % area fraction

% Initialize particles radomly with no overlaps
[x_0, L] = GenerateConfig_2D(N, phi, 2);
x_0 = [x_0(:,1), zeros(N,1), x_0(:,2)]; % put particles in the xz plane
box = [L(1), L_y, L(2)]; % box dimensions

%% Run optimization

% Numerical parameters
step_size = 0.01;
error_tol = 0.001;
T = 0.10*[1, 0, 1];

% File name
outfile = sprintf('phasesep_T%.2f.mat', T(1));

% Run the design
[x, ext, E] = Design_Position(x_0, box, omega, E_0, ext_t, eps_params, step_size, T, error_tol, outfile);

end