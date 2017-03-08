function J = JacBuck(theta, n, mu)
  A = linBuck(n);
  h_inv = (n + 1)^2;
  J = A + spdiags(mu * cos(theta) / h_inv, 0, n, n);
end