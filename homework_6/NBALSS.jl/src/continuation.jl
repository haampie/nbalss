function continuation!(f::Function, ∂f::Function, θ::Vector{Float64}, μs)
  for μ = μs
    b = x -> f(x, μ)
    ∂b = x -> ∂f(x, μ)
    newton!(b, ∂b, θ)
  end
end

function continuation_with_stability_check!(f::Function, ∂f::Function, θ::Vector{Float64}, μs)
  for μ = μs
    # Do a continuation step
    b = x -> f(x, μ)
    ∂b = x -> ∂f(x, μ)
    newton!(b, ∂b, θ)
    
    # Find three eigenvalues close to the origin
    F = factorize(∂f(θ, μ))
    _, D, _ = orthogonal_iteration(F, 3, Invert())

    num_positive = count(x -> x > 0, D)
    println("$μ")
    println(num_positive > 0 ? "Unstable" : "Stable")
  end
end