function [Fr] = RealSpaceForce(x,m,box,p1,p2,rc,F1,F2,rvals)

% INPUTS
% x = (N-by-3 array) particle positions
% m = (N-by-3 array) particle magnetic dipoles
% box = (3 row vector) box dimensions
% p1 = (column vector) first particles in the neighbor pair list
% p2 = (column vector) second particles in the neighbor pair list
% rc = real space cutoff radius
% F1 = (column vector) coefficient multiplying -(mi*mj)r and -( (mj*r)mi +
%      (mi*r)mj - 2(mi*r)(mj*r)r ) for real space force from derivatives
% F2 = (column vector) coefficient multiplying (mi*r)(mj*r)r for real space
%      force from derivatives
% rvals = (column vector) separation values for which to compute the real
%         space contribution
%
% OUTPUTS
% Fr = (N-by-3 array) real space contribution to the magnetic field

% Other quantities
N = size(x,1); % number of particles

% Initializations
Fr = zeros(N,3);

for i = 1:length(p1) 
    
    % Closest image distance between the particles
    r = x(p1(i),:) - x(p2(i),:);
    r = r - box.*fix(2*r./box);
    d = sqrt(r*r');
    r = r ./ d; % convert r to a unit vector
   
    % Compute potential and force if the distance is less than the cutoff distance
    if d < rc
        
        % Particle dipoles
        mi = m(p1(i),:);
        mj = m(p2(i),:);
        
        % Find the entries in the real space tables between which to
        % interpolate.
        interpind = find(floor(rvals/d),1);
        rhigh = rvals(interpind);
        rlow = rvals(interpind-1);
        
        % Interpolate between the values in the table.
        A = F1(interpind-1)*(rhigh - d)/(rhigh-rlow) + ...
            F1(interpind)*(d - rlow)/(rhigh-rlow);
        B = F2(interpind-1)*(rhigh - d)/(rhigh-rlow) + ...
            F2(interpind)*(d - rlow)/(rhigh-rlow);
        
        % Accumulate the p2(i)th particle's contribution to the p1(i)th
        % particle's potential
        Fr(p1(i),:) = Fr(p1(i),:) - A*( (mi*mj.')*r+(mj*r')*mi+(mi*r')*mj ) + (2*A-B)*(mi*r')*(mj*r')*r;
    end 
     
end

end