function testdeflat()
  n = 100;
  [B, rhs] = problem(n, 10);
  [L, U] = preconditioner(B);
  v = approximate_bad_eigenvector(B, 10);

  % Projects z on R^n \ span{v}
  P1 = @(z) z - dot(v, z) * v;

  % Restrict B to R^n \ span{v}
  B1 = @(z) P1(B * P1(z));

  % Solves P1 B P1 x = b for x.
  solver = @(b) pcg(B1, b, 1e-12, 400, L, U);

  % Other projection and its transpose
  P2 = @(z) z - P1(solver(P1(B * z)));
  P2_t = @(z) z - B * P1(solver(P1(z)));

  % Restrict -B to image of P2
  B2 = @(z) -P2_t(B * P2(z));
  
  [x1, ~, ~, ~, history_1] = solver(P1(rhs));
  [x2, ~, ~, ~, history_2] = pcg(B2, -P2_t(rhs), 1e-12, 400);

  x = P1(x1) + P2(x2); 
  fprintf('|Residual| = %e\n', norm(B * x - rhs));
  
  semilogy(history_1); hold on
  semilogy(history_2);
end

function [L, U] = preconditioner(B)
  opts.type = 'ilutp';
  opts.droptol = 1.0;
  [L, U] = ilu(B, opts); %makes incomplete LU, in fact a diagonal matrix
end

function [B, rhs] = problem(n, mu)
  n = 100;  
  B = -linbuck(n) - mu * speye(n);  
  rhs = sin(0.01 * (1 : n)'); % An arbitrary smooth rhs.
end

function v = approximate_bad_eigenvector(B, m)
  m = 10;
  n = size(B, 1);
  P = kron(eye(m), ones(n / m, 1));
  [V, ~] = eig(P' * B * P);
  v = P * V(:, 1);
  v = v / norm(v);
end

function A = linbuck(n)
  e = ones(n, 1) * (n + 1) ^ 2;
  A = spdiags([e, -2 * e, e], -1 : 1, n, n);
end