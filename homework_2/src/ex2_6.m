function ex2_6(n, mu, tol)
  % Mesh width
  h = 1 / (n + 1);

  % Internal grid points
  xs = linspace(h, 1 - h, n).';

  % Do some Newton iterations
  theta_start = sin(sqrt(mu) * xs);
  my_zero = newton(@fBuck, @JacBuck, theta_start, tol, n, mu);

  figure
  plot(xs, my_zero); hold on;
  plot(xs, theta_start); hold off;
  legend('Approximate solution', 'Initial guess')
end

