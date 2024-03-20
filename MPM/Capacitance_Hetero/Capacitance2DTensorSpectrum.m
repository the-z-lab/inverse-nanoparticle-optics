function [C, p] = CapacitanceTensorSpectrum(x, box, eps_p, dim, xi, varargin)

% Compute the capacitance spectrum.
%
% INPUTS
% x = (N-by-3-by-N_frames) particle positions at each frame
% box = (1-by-3) box dimensions
% eps_p = (N-by-N_k) particle dielectric functions at each wavevector
% dim = (1-by-2) dimensions of incident polarization; 1 = x, 2 = y, 3 = z
% xi = (scalar) Ewald parameter
% p_guess = (N-by-3-by-2-by-N_k-by-N_frames) optional; initial guess for particle dipoles
%
% OUTPUTS
% C = (2-by-3-by-N_k) capacitance at each wavevector
% p = (N-by-3-by-2-by-N_k-by-N_frames) p(:,:,i,:,:) are the particle dipoles at each wavevector and frame for dim(i) polarization

% Set numerical parameters
errortol = 10^-3; % error tolerance

% Config parameters
N = size(x,1); % number of particles
N_frames = size(x,3); % number of frames
N_k = size(eps_p,2); % number of wavevectors
eta = 4*pi*N/(3*prod(box)); % volume fraction

% Build the real space table
r_table = (1:0.001:10)'; % separation values for real space table
[~, ~, field_dip_1, field_dip_2, ~, ~] = RealSpaceTable(r_table, xi);
r_table = [0; r_table]; % prepend 0 to r_table

% Pre-calculations
kcut = 2*xi*sqrt(-log(errortol)); % wave space cutoff
Ngrid = ceil(1+box*kcut/pi); % number of grid nodes in each dimension
h = box./Ngrid; % grid spacing in each dimension
P = ceil(-2*log(errortol)/pi); % number of grid nodes over which to support the Gaussians
[offset, offsetxyz] = PreCalculations(P, h);

% Initializations. Wavenumbers need to get the last index because they are
% in the parallel loop
C = zeros(2, 3, N_frames, N_k); % capacitance tensor
p = zeros(N, 3, 2, N_frames, N_k); % particle dipoles

% Loop through configurations
for j = 1:N_frames
    
    % Loop through wavevectors
    parfor i = 1:N_k

        if isempty(varargin) 
            
            % Construct initial guess
            beta = (eps_p(:,i)-1)./(eps_p(:,i)+2);
            p_guess = zeros(N,3,2);
            p_guess(:,dim(1),1) = 4*pi*beta./(1-beta*eta);
            p_guess(:,dim(2),2) = 4*pi*beta./(1-beta*eta);
            
        else
            p_guess = varargin{1}(:,:,:,i,j);
        end
        
        % Compute the capacitance
        [C(:,:,j,i),p(:,:,:,j,i)] = Capacitance_Tensor(x(:,:,j), eps_p(:,i), box, dim, p_guess, xi, errortol, field_dip_1, field_dip_2, r_table, offset, offsetxyz);
        
    end
end

% Average the capcitance over configurations
C = squeeze(mean(C,3));

% Switch the 4th and 5th indices; now frames are the last index
p = permute(p,[1,2,3,5,4]);

end