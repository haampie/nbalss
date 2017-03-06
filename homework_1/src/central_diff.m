function A = central_diff(n)
  e = ones(n, 1);
  A = spdiags([e, -2 * e, e], -1 : 1, n, n);
end