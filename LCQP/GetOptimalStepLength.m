function [ alpha ] = GetOptimalStepLength(C, g, Rk, Q, pk, xk)

% Quadratic matrix
Qk = Q + Rk*C;

% Quadratic term
qk = pk'*Qk*pk;

% Linear term
lk = pk'*(Qk*xk + g);

% Step size function along step
meritAlongStep = @(a) a^2*qk + a*lk;

if (qk <= eps)
    % Handle not convex case
    if (meritAlongStep(1) < meritAlongStep(0))
        alpha = 1;
    else
        alpha = 0;
    end
elseif (lk > 0)
    % Handle no descent case
    alpha = 0;
else
    % At most take step equal to one.
    alpha = min(-lk/qk, 1);
end

% Handle no step case
if (alpha == 0)
    alpha = 0;
end

end