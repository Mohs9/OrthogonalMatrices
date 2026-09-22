%% Local function to wrap angles
function theta = wrap_to_pi_local(theta)
theta = mod(theta + pi, 2*pi) - pi;
end


