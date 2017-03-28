include("newton.jl")
include("continuation.jl")
include("nonlinear_backward_euler.jl")
include("pde_and_grid.jl")

using Plots, GR

gr()

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
  nonlinear_backward_euler!(fμ, ∂fμ, θ₂, 100., Δt = 5.)
  
  # Plot the solutions
  ss = grid_interior(n);
  Plots.plot(ss, θ, label = "Easy solution")
  Plots.plot!(ss, θ₂, label = "Isolated solution")
end

function ex5_4(; n = 100, ɛ = 0.01)
  f, ∂f = quadratic(n, ɛ)
  ss = grid_interior(n)
  θ_norms_time = Float64[]
  θ_norms_cont = Float64[]
  μ_end = 14.0;
  μs = linspace(6.0, μ_end, 100)

  # Time integration
  for μ = μs
    fμ = x -> f(x, μ)
    ∂fμ = x -> ∂f(x, μ)
    θ = ɛ * sin(μ * ss)
    nonlinear_backward_euler!(fμ, ∂fμ, θ, 100.0, Δt = 5.0)
    push!(θ_norms_time, norm(θ))
  end

  # Continuation
  θ₂ = zeros(n)
  μs_cont = linspace(6.0, μ_end, 300)
  for μ = μs_cont
    b = x -> f(x, μ)
    ∂b = x -> ∂f(x, μ)
    newton!(b, ∂b, θ₂)
    push!(θ_norms_cont, norm(θ₂))
  end

  Plots.plot(μs, θ_norms_time, xlabel = "\$\\mu\$", ylabel = "norm \\theta", label = "Time integration")
  Plots.plot!(μs_cont, θ_norms_cont, label = "Time-independent continuation")
end