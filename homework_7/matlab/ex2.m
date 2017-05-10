function ex2(n, k, l)
  A = spdiags(ones(n, 1) * (n + 1)^2 * [-1 2 -1], -1 : 1, n, n);

  tol = 1e-10;

  xs = cheby_nodes(l, 0, 50);
  shifts = cell(l, 1);

  for i = 1 : l
    [L, U] = lu(A - xs(i) * speye(n));
    shifts{i} = @(x) U \ (L \ x);
  end

  mv_A = @(x) applyA(shifts, x);

  [D, V] = orthogonal_iteration(mv_A, n, k, 50);

  reduction = abs(D(2 : end, :) - D(1 : end - 1, :));

  mvs = (nnz(reduction > tol) + k) * l;

  mvs

  for i = 1 : k
    history = reduction(:, i);
    history(history < tol) = 1e-16;
    semilogy(history); hold on
  end
  hold off;
end

function y = applyA(shifts, x)
  y = x;
  for idx = 1 : length(shifts)
    y = shifts{idx}(y);
  end
end

function xs = cheby_nodes(k, a, b)
  xs = (a + b) / 2.0 + (b - a) / 2 * cos((2 * (k : -1 : 1) - 1) * pi / (2 * k));
end

function [approx, V] = orthogonal_iteration(mv_A, n, k, iter)
  V = rand(n, k);
  [V, ~] = qr(V, 0);

  approx = zeros(iter, k);

  for i = 1 : iter
    B = mv_A(V);
    approx(i, :) = sort(eig(V' * B));
    [V, ~] = qr(B, 0);
  end
end