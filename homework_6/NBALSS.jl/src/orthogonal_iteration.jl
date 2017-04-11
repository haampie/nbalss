function apply_A!(method::Multiply, Z, A, V)
  A_mul_B!(Z, A, V)
end

function apply_A!(method::Invert, Z, A, V)
  A_ldiv_B!(Z, A, V)
end

function approx_eigenvals(method::Multiply, B)
  sort!(eigvals(B), by = abs, rev = true)
end

function approx_eigenvals(method::Invert, B)
  1 ./ sort!(eigvals(B), by = abs, rev = true)
end

function orthogonal_iteration{T<:Method}(A, k::Int, method::T = Multiply(), max_iter = 35)
  n = size(A, 1)

  # Random orthonormal basis vectors
  V, = qr(rand(n, k))

  # Initialize Z and B
  Z = similar(V)
  B = zeros(k, k)

  history = zeros(max_iter, k)

  for j = 1 : max_iter
    # Z = A * V or A \ V
    apply_A!(method, Z, A, V)

    # B = V' * Z = V' * A * V
    Ac_mul_B!(B, V, Z)

    # History
    history[j, :] = approx_eigenvals(method, B)

    V, = qr(Z)
  end

  V, history[max_iter, :], history
end
