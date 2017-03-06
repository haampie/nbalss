function ex1(n)
  % Only solve for n - 1 interior
  % nodes.

  h = 1 / n;
  xs = linspace(0, 1, n + 1).';
  mu_sqrt = pi;
  mu = mu_sqrt ^ 2;

  theta = sin(mu_sqrt * xs);
  perturb = rand(n + 1, 1);
  perturb(1) = 0;
  perturb(end) = 0;
  epsilon = 0.0001;

  fst = fBuck(theta, n, mu);
  snd = fBuck(theta + epsilon * perturb, n, mu);

  norm((fst - snd) - epsilon * perturb .* JacBuck(theta, n, mu))
end

function fB = fBuck(theta, n, mu)
  A = central_diff(n - 1, 1 / n);
  inner = theta(2 : end - 1);
  fB = [0; A * inner + mu * sin(inner); 0];
end

function J = JacBuck(theta, n, mu)
  J = mu * cos(theta);
end