function A = linBuck(n)
  % Discrete Poisson matrix
  e = ones(n, 1);
  A = spdiags([e, -2 * e, e], -1 : 1, n, n);
end