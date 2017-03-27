function newton!(f::Function, ∂f::Function, x::Vector{Float64}; ɛ = 1e-10)
  Δx_norm = Inf
  history = Float64[]

  while Δx_norm > ɛ
    Δx = ∂f(x) \ f(x)
    x[:] = x - Δx
    Δx_norm = norm(Δx)
    push!(history, Δx_norm)
  end

  history
end

function nonlinear_backward_euler!(f::Function, ∂f::Function, θ::Vector{Float64}, T::Float64 = 10.; Δt::Float64 = 0.1)
  I = Diagonal(ones(θ))

  for t = Δt : Δt : T
    newton!(x -> x - θ - Δt * f(x), x -> I - Δt * ∂f(x), θ)
  end
end

function grid_interior(n)
  xs = linspace(0, 1, n + 2);
  xs[2 : end - 1];
end

function poisson(n)
  t = (n + 1)^2
  SymTridiagonal(-2 * t * ones(n), t * ones(n - 1))
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

function ex5_1(;μ = 27.0)
  n = 100
  ɛ = 0.01
  f, ∂f = quadratic(n, ɛ)
  ss = grid_interior(n)
  
  θ = ɛ * sin(π * ss)
  
  Δt = 0.00001
  T = 0.0001
  for t = Δt : Δt : T
    println(θ)
    nonlinear_backward_euler!(x -> f(x, μ), x -> ∂f(x, μ), θ, 1.1 * Δt, Δt = Δt)
    Plots.plot!(θ)
  end


  
end
