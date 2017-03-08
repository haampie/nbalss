function f = fBuck(theta, n, mu)
  A = linBuck(n);
  h_inv = (n + 1)^2;
  f = A * theta + mu * sin(theta) / h_inv;
end