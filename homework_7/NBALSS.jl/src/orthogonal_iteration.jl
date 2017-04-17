function orthogonal_iteration(A, k::Int, max_iter = 35)
  n = size(A, 1)

  # Random orthonormal basis vectors
  V, = qr(rand(n, k))

  # Initialize Z and B
  Z = similar(V)
  B = zeros(k, k)

  history = zeros(max_iter, k)

  for j = 1 : max_iter
    # Z = A * V
    A_mul_B!(Z, A, V)

    # B = V' * Z = V' * A * V
    Ac_mul_B!(B, V, Z)

    # History
    history[j, :] = sort!(eigvals(B))

    V, = qr(Z)
  end

  V, history[max_iter, :], history
end
