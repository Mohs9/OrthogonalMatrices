function [normA4,xBest] = norm_4(A,nStart)

if nargin < 2
    nStart = 1000;
end

K = size(A,2);

% Maximize sum((A*x).^4)
% fmincon minimizes, hence the minus sign
objective = @(x) -sum((A*x).^4);

% L4 unit sphere constraint
nonlcon = @(x) norm_4_constraint(x);

options = optimoptions('fmincon', ...
    'Display', 'off', ...
    'Algorithm', 'sqp', ...
    'MaxIterations', 5000, ...
    'OptimalityTolerance', 1e-10, ...
    'StepTolerance', 1e-12);

bestValue = -Inf;
xBest = [];

for s = 1:nStart

    % Random starting point on ||x||_4 = 1
    z = randn(K,1);
    x0 = z / norm(z,4);

    [x,fval] = fmincon( ...
        objective, ...
        x0, ...
        [],[],[],[],[],[], ...
        nonlcon, ...
        options);

    value = -fval;

    if value > bestValue
        bestValue = value;
        xBest = x;
    end

end

% Recover the induced norm
normA4 = bestValue^(1/4);

end


