function [A, b, x_grid] = discretize(rhs, n)
  % Discretizes the problem using n rather than n + 1
  % grid points, omitting the dirichlet grid point at x_0.

  h = 1 / n;
  x_grid = linspace(0, 1, n + 1)';
  ys = rhs(x_grid);

  % Central difference matrix A + r.h.s. b
  A = central_diff(n);
  b = h * h * ys(2 : end);
  
  % Dirichlet left
  b(1) = b(1) - 1;

  % Neumann right
  b(end) = b(end) - 2 * h * exp(1);
  A(n, n - 1) = 2;
end