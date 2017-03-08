function x = newton(f, df, x, tol, varargin)
  incr = inf;

  while norm(incr, 2) > tol
    incr = df(x, varargin{:}) \ f(x, varargin{:});
    x = x - incr;
  end
end