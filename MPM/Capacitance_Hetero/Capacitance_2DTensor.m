function [C,S] = Capacitance_2DTensor(x,lambdap,box,dim,Sguess,xi,errortol,dipoletable1,dipoletable2,rtable,offset,offsetxyz)

% Compute the two-dimensional capacitance tensor of a configuration of particles.
%
% INPUTS
% x = (N-by-3) particle positions
% lambdap = (N-by-1) particle conductivities
% box = (1-by-3) box dimensions
% dim = (1-by-2) dimensions of incident polarization; 1 = x, 2 = y, 3 = z
% Sguess = (N-by-3-by-2) initial guess for particle dipoles for field in ith direction
% xi = (scalar) Ewald parameter
% errortol = (scalar) error tolerance
% dipoletable1 = (Ntable vec) I-rr component of the field/dipole real space coupling
% dipoletable2 = (Ntable vec) rr component of the field/dipole real space coupling
% rtable = (Ntable vec) separation values for table entries
% offset = (P^3-by-3) indices of the P^3 surrounding nodes relative to a center node at [0,0,0]
% offsetxyz = (P^3-by-3) coordinates of the P^3 surrounding nodes relative to a center node at [0,0,0]
%
% OUTPUTS
% C = (2-by-3) effective capacitance
% S = (N-by-3-by-2) S(:,:,i) are the induced dipoles for field in ith direction

% Initializations
H0 = zeros(2,3);
H0(1, dim(1)) = 1; % unit magnetic field in the dim(1) direction
H0(1, dim(2)) = 1; % unit magnetic field in the dim(2) direction
C = zeros(2,3); % capacitance tensor
S = zeros(size(x,1),3,2); % particle dipoles

% Parameters for the spectral Ewald method
rc = sqrt(-log(errortol))/xi; % real space cutoff radius
kcut = 2*xi*sqrt(-log(errortol)); % wave space cutoff
Ngrid = ceil(1+box*kcut/pi); % number of grid nodes in each dimension
h = box./Ngrid; % grid spacing in each dimension
P = ceil(-2*log(errortol)/pi); % number of grid nodes over which to support the Gaussians
eta = P*(h*xi).^2/pi; % spectral splitting parameter

% Throw a warning if the real space cutoff is larger than half the box size
% in any dimension
if any(rc > box/2)
    error('Real space cutoff (rc = %.3f) larger than half the box size (%.3f, %.3f, %.3f). Accuracy to specified error tolerance not guaranteed.',rc,box(1),box(2),box(3))
end

% Compute the neighbor list.
[cell,Ncell] = CellList(x,box,rc);
[p1,p2] = NeighborList(cell,Ncell);

% Loop through two orthogonal field directions
for i = 1:2
    
    % Compute the magnetic dipoles
    S_i  = MagneticDipole(x,lambdap,H0(i,:),box,p1,p2,Ngrid,h,P,xi,eta,rc,...
                   Sguess(:,:,i),offset,offsetxyz,dipoletable1,dipoletable2,rtable,errortol);
    S(:,:,i) = S_i;

    % Compute the capcitance from the average dipole
    C(i,:) = mean(S_i,1);    
    
end

end