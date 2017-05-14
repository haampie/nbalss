function ex8_1()
  n = 100;
  A = central_diff(n);

  for mu = [0, -5, -10, -20]
    Abar = A - mu * speye(n);
    % LU = my_lu(Abar);
    [L, U, P] = lu(Abar);

    P' * P - speye(n)

    % spy(L)

    % P

    % disp(norm(L * U - P * Abar, 'fro'))

    % spy(L + U)

    pause



    % [L, U, P] = lu(Abar);
    % % L = chol(-Abar);
    % spy(L)
    % full(L)
    % % disp(norm(Abar - L * U * P, 'fro'))
    % pause
  end
end

function A = central_diff(n)
  e = ones(n, 1) * (n + 1) ^ 2;
  A = spdiags([e, -2 * e, e], -1 : 1, n, n);
end


function A = my_lu(A)
  n = size(A, 1);

  for k = 1 : n - 1
    if A(k, k) == 0;
      error('Null pivot element')
    end

    if abs(A(k, k)) < abs(A(k + 1, k))
      A([k, k+1], :) = A([k+1, k], :);
      disp(sprintf('Swapping %d and %d\n', k, k + 1));
    end

    A(k + 1 : n, k) = A(k + 1 : n, k) / A(k, k);
    i = [k + 1 : n];
    for j = i
      A(i, j) = A(i, j) - A(i, k) * A(k, j);
    end
  end
end