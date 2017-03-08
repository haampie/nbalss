function f = fBuck(theta, n, mu)
  A = linBuck(n);

  % Here h = 1 / (n + 1)^2
  f = A * theta + mu * sin(theta) / (n + 1)^2;
end