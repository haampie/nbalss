function res = ex1(n, k)
  A = spdiags(ones(n, 1) * (n + 1)^2 * [-1 2 -1], -1 : 1, n, n);
  
  C = chol(A);
  mv = @(x) C \ (C' \ x);
  % eigs(mv, n, k, 'sm', opts)

  B = spdiags(ones(n, 1), 0, n, n);
  B(n, n) = 0;

  res = eigs(A, B, k, 'sm')
  res = eigs(mv, n, B, k, 'sm')
end