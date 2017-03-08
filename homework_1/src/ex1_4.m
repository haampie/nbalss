function ex1_4
  % Exact solution
  exact = @(x) exp(x);

  % Solve for one value of N to plot the error
  [A, b, x_grid] = discretize(@(x) exp(x), 32);
  
  % Solve Ax = b and prepend the Dirichlet grid point.
  numerical = [1; A \ b];

  % Local error.
  err = numerical - exact(x_grid);

  figure;
  plot(x_grid, abs(err), 'r.-');
  title('Local error')
  xlabel('x')
  ylabel('error')

  fprintf('Global error (inf norm) = %f\n', norm(err, Inf));
end

