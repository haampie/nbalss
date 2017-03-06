function ex1_5
  es = [];
  ns = 2 .^ (3 : 11);
  exact = @(x) exp(x);

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
  ylabel('Global error in inf norm');
  legend('Error', '1 / N^2')
end