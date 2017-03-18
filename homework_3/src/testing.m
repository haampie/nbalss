function testing
  thetas = init(6, 100, 1e-6);

  k = 6;

  first = thetas(:, k);

  from = (k * pi)^2 * 1.1;
  to = from + 500;
  steps = 50;

  cont(from, first, to, steps, 1e-6)

  % plot(thetas)
end