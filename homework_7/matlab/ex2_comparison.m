function ex2_comparison(n, k)
  A = spdiags(ones(n, 1) * (n + 1)^2 * [-1 2 -1], -1 : 1, n, n);

  opts.disp = 1;
  opts.p = 1;
  opts.maxit = 1;

  eigs(A, k, 'sm', opts)
end
