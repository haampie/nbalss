function newton!(f::Function, ∂f::Function, x::Vector{Float64}; ɛ = 1e-10)
  Δx_norm = Inf
  history = Float64[]

  while Δx_norm > ɛ
    Δx = ∂f(x) \ f(x)
    x .= x - Δx
    Δx_norm = norm(Δx)
    push!(history, Δx_norm)
  end

  history
end

function nonlinear_backward_euler!(f::Function, ∂f::Function, θ::Vector{Float64}, T = 3.; Δt = 0.01)
  I = Diagonal(ones(θ))
  θₙ = copy(θ)

  for t = Δt : Δt : T
    g = θₙ₊₁ -> θₙ₊₁ - θₙ - Δt * f(θₙ₊₁);
    ∂g = θₙ₊₁ -> I - Δt * ∂f(θₙ₊₁);
    newton!(g, ∂g, θ)
  end
end

function continuation!(f::Function, ∂f::Function, θ::Vector{Float64}, μs::LinSpace{Float64})
  for μ = μs
    b = x -> f(x, μ)
    ∂b = x -> ∂f(x, μ)
    newton!(b, ∂b, θ)
  end
end

function grid_interior(n::Int)
  xs = linspace(0, 1, n + 2)
  xs[2 : end - 1]
end

function poisson(n::Int)
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

function ex5_1(;μ = 20.0, n = 100, ɛ = 0.01)
  f, ∂f = quadratic(n, ɛ)
  ss = grid_interior(n)
  θ = ɛ * sin(π * ss)

  fμ = x -> f(x, μ)
  ∂fμ = x -> ∂f(x, μ)
  nonlinear_backward_euler!(fμ, ∂fμ, θ, 10.0, Δt = 0.1)

  Plots.plot(ss, θ)
end

function ex5_3(; n = 100, ɛ = 0.01, μ_end = 20.)
  # Initialize the problem.
  θ = zeros(n)
  f, ∂f = quadratic(n, ɛ)

  # First do continuation from 0 to 20
  continuation!(f, ∂f, θ, linspace(0, μ_end, 400))
  
  # Minus the "easy" solution
  θ₂ = -copy(θ)

  # Do time integration for fixed μ = μ_end.
  fμ = x -> f(x, μ_end);
  ∂fμ = x -> ∂f(x, μ_end);
  nonlinear_backward_euler!(fμ, ∂fμ, θ₂, 10., Δt = 0.1)
  
  # Plot the solutions
  ss = grid_interior(n);
  Plots.plot(ss, θ)
  Plots.plot!(ss, θ₂)
end

function ex5_4(; n = 100, ɛ = 0.01)
  f, ∂f = quadratic(n, ɛ)
  ss = grid_interior(n)
  θ_norms_time = Float64[]
  θ_norms_cont = Float64[]
  μ_end = 15.0;
  μs = linspace(0, μ_end, 100)

  # Time integration
  for μ = μs
    fμ = x -> f(x, μ)
    ∂fμ = x -> ∂f(x, μ)
    θ = ɛ * sin(π * ss)
    nonlinear_backward_euler!(fμ, ∂fμ, θ, 1000.0, Δt = 5.0)
    push!(θ_norms_time, norm(θ))
  end

  # Continuation
  θ₂ = zeros(n)
  μs_cont = linspace(0, μ_end, 1000)
  for μ = μs_cont
    b = x -> f(x, μ)
    ∂b = x -> ∂f(x, μ)
    newton!(b, ∂b, θ₂)
    push!(θ_norms_cont, norm(θ₂))
  end

  # Let's try some perturbation
  # θ₃ = copy(θ₂) + 0.1 * rand(n)
  # backup = copy(θ₃)
  # fμ = x -> f(x, μ_end)
  # ∂fμ = x -> ∂f(x, μ_end)
  # nonlinear_backward_euler!(fμ, ∂fμ, θ₃, 100.0, Δt = 5.0)
  # Plots.plot(θ₂)
  # Plots.plot!(θ₃)
  # Plots.plot!(backup)

  Plots.plot(μs, θ_norms_time)
  Plots.plot!(μs_cont, θ_norms_cont)
end