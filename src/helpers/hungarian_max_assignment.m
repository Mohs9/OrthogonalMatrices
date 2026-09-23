function assignment = hungarian_max_assignment(score)
%HUNGARIAN_MAX_ASSIGNMENT Maximize sum_i score(i, assignment(i)).

if any(~isfinite(score(:)))
    error('The assignment score matrix must have finite entries.');
end

K = size(score, 1);

% Convert maximization to minimization. Shifting by max(score) keeps costs
% nonnegative and preserves the optimizer.
cost = max(score(:)) - score;

u = zeros(K + 1, 1);
v = zeros(K + 1, 1);
p = zeros(K + 1, 1);
way = zeros(K + 1, 1);

for i = 1:K
    p(1) = i;
    j0 = 1;
    minv = Inf(K + 1, 1);
    used = false(K + 1, 1);
    way(:) = 0;

    while true
        used(j0) = true;
        i0 = p(j0);
        delta = Inf;
        j1 = 1;

        for j = 2:K + 1
            if used(j)
                continue
            end

            cur = cost(i0, j - 1) - u(i0) - v(j);

            if cur < minv(j)
                minv(j) = cur;
                way(j) = j0;
            end

            if minv(j) < delta
                delta = minv(j);
                j1 = j;
            end
        end

        for j = 1:K + 1
            if used(j)
                u(p(j)) = u(p(j)) + delta;
                v(j) = v(j) - delta;
            else
                minv(j) = minv(j) - delta;
            end
        end

        j0 = j1;

        if p(j0) == 0
            break
        end
    end

    while true
        j1 = way(j0);
        p(j0) = p(j1);
        j0 = j1;

        if j0 == 1
            break
        end
    end
end

assignment = zeros(1, K);

for j = 2:K + 1
    assignment(p(j)) = j - 1;
end
end