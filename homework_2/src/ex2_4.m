function ex2_4(n, mu)
  % Mesh width
  h = 1 / (n + 1);

  % Internal grid points
  xs = linspace(h, 1 - h, n).';

  % Eigensolution of linearized system
  theta = sin(sqrt(mu) * xs);

  % Random unit vector
  v = rand(n, 1) - 0.5;
  v = v / norm(v);

  epsilons = 10 .^ -(1 : 12);
  errors = [];

  % Jacobi matrix in theta
  J = JacBuck(theta, n, mu);

  for epsilon = epsilons
    func = fBuck(theta + epsilon * v, n, mu);
    approx = fBuck(theta, n, mu) + epsilon * J * v;
    errors(end + 1) = norm(func - approx, Inf);
  end

  figure;
  loglog(1 ./ epsilons, epsilons .^ 2, '-+'); hold on;
  loglog(1 ./ epsilons, errors, '-+'); hold off;
  grid on
  xlabel('e^{-1}')
  legend('e^2', 'r(\epsilon)')
end
