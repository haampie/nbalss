function thetas = init(k, n, tol)
  % larger n -> more accurate eigenpairs.
  h = 1 / (n + 1);

  % Find eigenpairs.
  A = linBuck(n);
  [V, D] = eig(full(A));
  eigen = diag(D);
  [~, I] = sort(abs(eigen));
  eigen = eigen(I);
  V = V(:, I);
  mus = eigen / -h ^ 2;

  % Given epsilon, compute gamma
  % and the initial theta's with
  % corresponding mu's.
  correction = sqrt((n + 1) / 2);
  gamma = 0.1
  epsilon = 8 * gamma / (1 + gamma)
  theta_normalized = correction * V;
  thetas_0 = epsilon * theta_normalized;
  mus_0 = mus * (1 + gamma);

  % Do the Newton iterations.
  % Let the user add the Dirichlet boundaries :).
  thetas = zeros(n, k);
  for idx = 1 : k
    thetas(:, idx) = zeroOfBuck(n, mus_0(idx), tol, thetas_0(:, idx));
  end
end