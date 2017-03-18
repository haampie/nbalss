function J = JacBuck(theta, n, mu)
  A = linBuck(n);
  % h ^ 2 = 1 / (n + 1) ^ 2
  diag_vec = mu * cos(theta) / (n + 1) ^ 2;
  J = A + spdiags(diag_vec, 0, n, n);
end