function [p, E, grad_E2] = AnalyticTwoBody(x, alpha, box, E_0, x_E)

% Separation between the two particles
d_vec = x(1,:) - x(2,:); % points from x_2 to x_1
d_vec = d_vec - box.*fix(2*d_vec./box); % nearest image distance
d = norm(d_vec);
d_hat = d_vec/d;

% Dipole moments
A = 1/(1-2*alpha/d^3);
B = 1/(1+alpha/d^3);

p = 4*pi*alpha*( (A-B)*d_hat*(d_hat*E_0.') + B*E_0 );
p = [p; p]; % 2-by-3

% Dipole gradient w.r.t. x_1
A = 1/d^2*1/(1-2*alpha/d^3);
B = -1/d^2*1/(1+alpha/d^3);
A_prime = -2*(d^3+alpha)/(d^3-2*alpha)^2;
B_prime = (2*d^3-alpha)/(d^3+alpha)^2;
C_prime = 3*d^2*alpha/(d^3+alpha)^2;

grad_p = 4*pi*alpha*( d^2*(A_prime+B_prime)*d_hat'.*d_hat*(d_hat*E_0.') + d*(A+B)*(eye(3)*(d_hat*E_0.') + E_0.'*d_hat) + C_prime*d_hat'.*E_0 );

% Initialize
E = 0;
grad_E = 0;

% Loop through particles
for i = 1:size(x,1)
    
    % Separation between probe point and particles
    r_vec = x_E - x(i,:); % points from x_i to probe point
    r_vec = r_vec - box.*fix(2*r_vec./box); % nearest image distance
    r = norm(r_vec);
    r_hat = r_vec/r;
    
    % Grab the ith dipole
    p_i = p(i,:);
    
    % Scattered field
    M_sc = 1/(4*pi*r^3)*(3*r_hat'.*r_hat - eye(3)); % 3-by-3
    E_i = (M_sc*p_i.').'; % 1-by-3
    
    % Field gradient: contribution from grad M_sc; only grad_i M_sc,i is nonzero
    if i == 1
        grad_E_i = -3/(4*pi*r^4)* ( r_hat'.*p_i + p_i.'.*r_hat + eye(3)*(r_hat*p_i.') - 5*r_hat'.*r_hat*(r_hat*p_i.') ); % 3-by-3; negative sign b/c r = x_E - x_i and grad_x_i = -grad_r
    else
        grad_E_i = 0;
    end
    
    % Field gradient: contribution from grad p
    grad_E_i = grad_E_i + grad_p*M_sc;
    
    % Accumulate
    E = E + E_i;
    grad_E = grad_E + grad_E_i;
    
end

% Gradient of |E|^2 w.r.t. x_1
grad_E2 = grad_E*conj(E).' + conj(grad_E)*E.';


end