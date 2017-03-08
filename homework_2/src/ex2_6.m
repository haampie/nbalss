function ex2_6(n, mu, tol)
  h = 1 / (n + 1);
  xs = linspace(0, 1, n + 2).';
  xs = xs(2 : end - 1);

  theta_start = sin(sqrt(mu) * xs);
  theta_zero = newton(@fBuck, @JacBuck, theta_start, tol, n, mu);

  figure
  plot(xs, theta_zero); hold on;
  plot(xs, theta_start); hold off;
  legend('Approximate solution', 'Initial guess')
end

