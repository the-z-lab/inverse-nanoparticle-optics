function [] = RunDesignLattice_Dielectric()

% Target parameters
omega_p_t = 1;
gamma_t = 0.10;
eps_inf_t = 1;

% Geometry
N = 25; % number of particles
N_layer = 3; % number of layers
h = 3; % spacing between layers
phi = 0.20; % area fraction
addpath('/Users/zacharysherman/Documents/Code/MATLAB/Lattices')
[x_0,L] = SQLattice(N, phi); % place particles on a 2D lattice
x_0 = [x_0, zeros(N,1)];
x = [];
type = [];

% Stack layers
for i = 1:N_layer
    x = [x; x_0 + (i-1)*[0, 0, h]];
    type = [type; i*ones(N,1)];
end

x = x + 1e-6; % slight perturbation to kick particles off the PBs
box = [L,L,50]; % box dimensions
E_0 = [1, 0, 0]; % field polarization

% Initial guess
a_0 = 1;
omega_p_0 = 1;
gamma_0 = 0.1;
eps_inf_0 = 1;
params_0 = [a_0, omega_p_0, gamma_0, eps_inf_0];
params_0 = params_0 + 0.05*(rand(N_layer,4)-1/2).*params_0;

% Specify which parameters are going to be varied by the optimizer (1);
% others (0) are held constant
design_flag = ones(N_layer, 4);

% Parameter bounds
min_bounds = [0.1, 0.1, 0.02, 0.1].*ones(N_layer,1);
max_bounds = [2, inf, inf, inf].*ones(N_layer,1);
bounds = cat(3, min_bounds, max_bounds);

% File name
outfile = 'test_2.mat';

% Numerical parameters
rel_steps = [1, 1, 0.125, 2.5];
step_size = 0.0001/N;
error_tol = 0.001;
iter_tol = 5000;

% Compute the target extinction spectrum
omega = (0.2:0.01:0.8);
% omega_res = omega_p_t/sqrt(2+eps_inf_t);
% h = 9*omega_p_t^2/(gamma_t*(2+eps_inf_t)^2);
% W = gamma_t;
% ext_t = h./(1+(omega_res/W)^2*(omega.'/omega_res - omega_res./omega.').^2);
ext_t = 3*(omega >= 0.3 & omega <= 0.7).';

% Run the design
[params, ext, E] = Design_Dielectric(x, type, box, omega.', E_0, ext_t, params_0, step_size, error_tol, outfile, design_flag, bounds, iter_tol);

end