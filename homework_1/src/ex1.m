function ex1_4
  es = [];
  ns = 2 .^ (3 : 11);
  exact = @(x) exp(x);

  % Solve for one value of N to plot the error
  [A, b, x_grid] = discretize(@(x) exp(x), 32);
  numerical = [1; A \ b];

  figure;
  plot(numerical - exact(x_grid));
  title('')

  % Show experimentally 2nd-order convergence
  for n = ns
    
    [A, b, x_grid] = discretize(@(x) exp(x), n);

    % Solve Ax = b and prepend the Dirichlet grid point.
    numerical = [1; A \ b];

    % Compute the normed error (max norm)
    es(end + 1) = norm(numerical - exact(x_grid), Inf);
  end

  figure;
  loglog(ns, es, 'b'); hold on;
  loglog(ns, ns .^ -2, '-.'); hold off;
  grid
  title('Convergence rate central differences');
  xlabel('N');
  ylabel('Error');
  legend('Error', '1 / N^2')
end

function A = central_diff(n)
  e = ones(n, 1);
  A = spdiags([e, -2 * e, e], -1 : 1, n, n);
end

function [A, b, x_grid] = discretize(rhs, n)
  % Discretizes the problem using n
  % grid points, omitting the dirichlet
  % grid point at x_0.

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