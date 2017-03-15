function my_zero = zeroOfBuck(n, mu, tol, theta_0)
  % Do some Newton iterations
  my_zero = newton(@fBuck, @JacBuck, theta_0, tol, n, mu);
end

