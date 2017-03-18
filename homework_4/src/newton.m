function [x, incr] = newton(f, df, x, tol, varargin)
  % Performs Newton-Raphson iterations until a prescribed
  % tolerance is met.

  % Matlab note: varargin collects the remaining arguments
  % by calling newton(f, df, x, tol, 1, 2, 3) varargin
  % contains the arguments {1, 2, 3}. But these additional
  % arguments are optional :).

  % Increment of previous iteration
  incr = [inf];

  % Loop as long as the increment is larger than the tolerance
  while incr(end) > tol
    % Pass the additional parameters to f and df.
    dx = df(x, varargin{:}) \ f(x, varargin{:});
    x = x - dx;
    incr(end + 1) = norm(dx, 2);
  end

  incr = incr(2 : end);
end