function init(n::Int, k::Int = 3; γ = 0.1)
  f, ∂f = problem(n)

  # Find the eigenvalues & vectors of the linearized system
  F = factorize(∂f(zeros(n), 0.0))
  V, D, = orthogonal_iteration(F, k, Invert());

  correction = √((n + 1) / 2);
  ɛ = 8 * γ / (1 + γ);
  θs = ɛ * correction * V;
  μs = -(1 + γ) * D;

  # Starting points
  starting_points = Tuple{Float64, Array{Float64}}[]

  for idx = 1 : k
    θ = θs[:, idx]
    newton!(x -> f(x, μs[idx]), x -> ∂f(x, μs[idx]), θ)
    push!(starting_points, (μs[idx], θ))
  end

  return starting_points
end