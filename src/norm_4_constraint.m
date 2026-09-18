function [c,ceq] = norm_4_constraint(x)

c = [];

ceq = sum(x.^4) - 1;

end