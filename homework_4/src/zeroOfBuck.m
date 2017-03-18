function [my_zero, history] = zeroOfBuck(n, mu, tol, theta_0)
  % Do some Newton iterations
  [my_zero, history] = newton(@fBuck, @JacBuck, theta_0, tol, n, mu);
end

