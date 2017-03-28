# Does in-place update of θ
function newton!(f::Function, ∂f::Function, x::Vector{Float64}; ɛ = 1e-5)
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