function ex1
  es = [];
  ns = 2 .^ (3 : 11);

  for n = ns
    h = 1 / n;
    x_grid = linspace(0, 1, n + 1)';
    ys = exp(x_grid);

    % Central difference matrix A + r.h.s. b
    A = central_diff(n);
    b = h * h * ys(2 : end);
    
    % Dirichlet left
    b(1) = b(1) - 1;

    % Neumann right
    b(end) = b(end) - 2 * h * exp(1);
    A(n, n - 1) = 2;

    % Solve Ax = b.
    x = [1; A \ b];

    % Compute the error (approximate L2-norm)
    es(end + 1) = L2_norm(x - ys);
  end

  loglog(ns, es, 'b'); hold on;
  loglog(ns, ns .^ -2, '-.'); hold off;
  grid
  title('Convergence rate central differences');
  xlabel('N');
  ylabel('Error');
  legend('Error', '1 / N^2')
end

function value = L2_norm(f)
  value = norm(f(1 : end - 1), 2) / sqrt(length(f) - 1);
end

function A = central_diff(n)
  e = ones(n, 1);
  A = spdiags([e, -2 * e, e], -1 : 1, n, n);
end