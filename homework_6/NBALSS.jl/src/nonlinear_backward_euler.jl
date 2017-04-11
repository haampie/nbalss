function nonlinear_backward_euler!(f::Function, ∂f::Function, θ::Vector{Float64}, T = 3.; Δt = 0.01)
  I = Diagonal(ones(θ))

  for t = Δt : Δt : T
    θₙ = copy(θ)
    g = θₙ₊₁ -> θₙ₊₁ - θₙ - Δt * f(θₙ₊₁);
    ∂g = θₙ₊₁ -> I - Δt * ∂f(θₙ₊₁);
    newton!(g, ∂g, θ)
  end
end