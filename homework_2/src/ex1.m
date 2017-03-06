function ex1(n, mu)
  close all
  % Let s_i = (i - 1/2)h for i = 1, 2, ..., n.
  % and h := 1 / n the mesh width.

  h = 1 / n;
  xs = linspace(h / 2, 1 - h / 2, n).';

  % mu_sqrt = 9 * pi;
  % mu = mu_sqrt ^ 2;
  theta = sin(sqrt(mu) * xs);

  A = snd_derivative(n, h);
  v = sin(10 * pi * xs);
  v = rand(n, 1);
  v = v / norm(v);

  denominators = 10 .^ (1 : 12);
  errors = [];
  jacobi = JacBuck(theta, n, mu);
  
  for denominator = denominators
    epsilon = 1 / denominator;
    dif = fBuck(theta + epsilon * v, n, mu) - fBuck(theta, n, mu);
    lin = epsilon * jacobi * v;
    errors(end + 1) = norm(dif - lin);
  end

  figure('Name','Error');
  loglog(denominators, 1 ./ denominators .^ 2, '-+'); hold on;
  loglog(denominators, errors, '-+'); hold off;
  grid on
  xlabel('e^{-1}')
  ylabel('Residual')
  legend('e^2', 'error')


  % x = newton(@(x, a, b) x^2 - 1, @(x, a, b) 2 * x, 0.1, 1e-4, n, mu)

  pert = rand(n, 1) - 0.5;
  pert = pert / norm(pert) * norm(theta);
  theta_0 = theta + 0.0 * pert;
  [improvement, it] = newton(@fBuck, @JacBuck, theta_0, 1e-9, n, mu);

  figure
  plot(xs, improvement);
  title(sprintf('Solution after it = %d', it));

  norm(fBuck(improvement, n, mu))
end

function fB = fBuck(theta, n, mu)
  A = snd_derivative(n, 1 / n);
  fB = A * theta + mu * sin(theta);
end

function J = JacBuck(theta, n, mu)
  A = snd_derivative(n, 1 / n);
  J = A + spdiags(mu * cos(theta), 0, n, n);
end

function A = snd_derivative(n, h)
  e = ones(n, 1) / h ^ 2;
  A = spdiags([e, -2 * e, e], -1 : 1, n, n);
  A(1, 1) = -4 / h ^ 2;
  A(1, 2) = 4 / (3 * h ^ 2);
  A(n, n - 1) = 4 / (3 * h ^ 2);
  A(n, n) = -4 / h ^ 2;
end

function [x, it] = newton(f, df, x, tol, p1, p2)
  it = 0;
  while norm(f(x, p1, p2), 2) > tol
    x = x - df(x, p1, p2) \ f(x, p1, p2);
    it = it + 1;
  end
end