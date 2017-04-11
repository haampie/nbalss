function continuation!(f::Function, ∂f::Function, θ::Vector{Float64}, μs::LinSpace{Float64})
  for μ = μs
    b = x -> f(x, μ)
    ∂b = x -> ∂f(x, μ)
    newton!(b, ∂b, θ)
  end
end