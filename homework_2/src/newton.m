function x = newton(f, df, x, tol, varargin)
  % Performs Newton-Raphson iterations until a prescribed
  % tolerance is met.

  % Matlab note: varargin collects the remaining arguments
  % by calling newton(f, df, x, tol, 1, 2, 3) varargin
  % contains the arguments {1, 2, 3}. But these additional
  % arguments are optional :).

  % Increment of previous iteration
  incr = inf;

  % Loop as long as the increment is larger than the tolerance
  while norm(incr, 2) > tol
    % Pass the additional parameters to f and df.
    incr = df(x, varargin{:}) \ f(x, varargin{:});
    x = x - incr;
  end
end