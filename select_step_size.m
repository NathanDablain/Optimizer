function alpha = select_step_size(z, delta_z, cons, cost_func)
  alpha = 1.0;

  % Loop through slack variables, find which will become the most negative, modify
  % delta_z so that it does not become negative
  tau = 0.95;
  z_new = z + delta_z;
  M = length(z);
  sv_start = M - cons.N_sv + 1;
  lowest_z = 1;
  lowest_index = 0;
  for i = sv_start : M
    if z_new(i) <= 0 && z_new(i) < lowest_z && delta_z(i) ~= 0.0
      lowest_index = i;
      lowest_z = z_new(i);
    end
  end

  if lowest_index ~= 0
    shaveoff = lowest_z;
    alpha = (delta_z(lowest_index) - shaveoff) / delta_z(lowest_index);
  end
  alpha = alpha * tau;

  % alpha is the max step we can take that will not violate our bounds
  % now perform a ternary search from 0 -> alpha to find the lowest cost
  alow = 0;
  ahigh = alpha;
  for i = 1:10
    amid1 = alow + (1/3)*(ahigh - alow);
    amid2 = alow + (2/3)*(ahigh - alow);
    [Jmid1, ~, ~] = cost_func(z + amid1*delta_z);
    [Jmid2, ~, ~] = cost_func(z + amid2*delta_z);
    %
    if Jmid1 <= Jmid2
      ahigh = amid2;
    elseif Jmid1 > Jmid2
      alow = amid1;
    end
  end
  alpha = amid1;
end

