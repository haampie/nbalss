function ex2_5
  % Solve f(x) = 0 for x.
  f = @(x) x^2 - 1;
  df = @(x) 2 * x;  % Derivative of f.

  % Print 16 digits of the approximation to the zero.
  fprintf('%.16f\n', newton(f, df, 0.1, 1e-4))
end