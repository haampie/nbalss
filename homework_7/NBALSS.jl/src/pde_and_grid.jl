function grid_interior(n::Int)
  xs = linspace(0, 1, n + 2)
  xs[2 : end - 1]
end

function poisson(n::Int)
  t = (n + 1)^2
  spdiagm((fill(t, n - 1), fill(-2 * t, n), fill(t, n - 1)), (-1, 0, 1))
end

# Returns a tuple with function and its Jacobian (f, ∂f).
function problem(n::Int)
  A = poisson(n)
  f = (θ::Vector{Float64}, μ::Float64) -> A * θ + μ * sin(θ)
  ∂f = (θ::Vector{Float64}, μ::Float64) -> A + Diagonal(μ * cos(θ))
  f, ∂f  
end

# Return f - ɛ s (1 - s) as function with Jacobian unchanged
function quadratic(n::Int, ɛ::Float64 = 0.01)
  ss = grid_interior(n)
  rhs = ɛ * ss .* (1 - ss)
  f, ∂f = problem(n)
  g = (θ::Vector{Float64}, μ::Float64) -> f(θ, μ) - rhs
  g, ∂f
end

