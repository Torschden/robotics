function F_handle = force_traj(profile)
% Returns a function F(t) implementing smooth force ramps
% based on time/value pairs using logistic transitions

    times = profile.times;
    values = profile.values;
    k = profile.k;
    t0 = profile.t0;

    % Validate lengths
    if length(times) ~= length(values)
        error('Times and values must have the same length');
    end

    % Build function handle
    F_handle = @(t) logistic_sum(t, times, values, k, t0);
end

function F = logistic_sum(t, times, values, k, t0)
    % Compute force value as a sum of weighted logistic ramps
    t = t(:);  % support vector inputs
    F = zeros(size(t));
    F_current = 0;

    for i = 1:length(times)
        F_next = values(i);
        delta = F_next - F_current;
        F = F + delta ./ (1 + exp(-k * (t - times(i))));
        F_current = F_next;
    end

    % Force is zero before t0
    F(t < t0) = 0;
end
