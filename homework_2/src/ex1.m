function ex1(n)
  % Let s_i = (i - 1/2)h for i = 1, 2, ..., n.
  % and h := 1 / n the mesh width.

  h = 1 / n;
  xs = linspace(h / 2, 1 - h / 2, n).';

  mu_sqrt = 2 * pi;
  mu = mu_sqrt ^ 2;
  theta = sin(mu_sqrt * xs);

  A = snd_derivative(n, h);
  v = sin(10 * pi * xs);
  v = rand(n, 1);
  v = v / norm(v);
  jacobi = JacBuck(theta, n, mu);


  denominators = 10 .^ (1 : 12);
  errors = [];
  
  for denominator = denominators
    epsilon = 1 / denominator;
    dif = fBuck(theta + epsilon * v, n, mu) - fBuck(theta, n, mu);

    lin = epsilon * jacobi .* [v, A * v] * ones(2, 1);

    errors(end + 1) = norm(dif - lin);
  end

  loglog(denominators, 1 ./ denominators.^2, '-+'); hold on;
  loglog(denominators, errors, '-+');


  x = newton(@(x) x^2 - 1, @(x) 2 * x, 0.1, 1e-7, n, mu)

end

function fB = fBuck(theta, n, mu)
  A = snd_derivative(n, 1 / n);
  fB = A * theta + mu * sin(theta);
end

function J = JacBuck(theta, n, mu)
  J = [mu * cos(theta), ones(n, 1)];
end

function A = snd_derivative(n, h)
  e = ones(n, 1) / h ^ 2;
  A = spdiags([e, -2 * e, e], -1 : 1, n, n);
  A(1, 1) = -4 / h ^ 2;
  A(1, 2) = 4 / (3 * h ^ 2);
  A(n, n - 1) = 4 / (3 * h ^ 2);
  A(n, n) = -4 / h ^ 2;
end

function x = newton(f, df, x, tol, p1, p2)
  while norm(f(x)) > tol
    x = x - f(x) / df(x)
  end
end