module Exercises
  
  importall NBALSS
  
  import Plots, PGFPlots

  Plots.pgfplots()

  export ex_1_convergence_history, ex1_init, ex3

  # Creates a plot of the convergence history
  # of the orthogonal_iteration method applied
  # on a 1d Poisson matrix.
  function ex1_convergence_history(n = 300, k = 3)
    # Apply inverted iteration with a factorized poisson matrix
    F = factorize(poisson(n))
    _, _, history = orthogonal_iteration(F, k, Invert())

    # Compute the difference |λᵢ - λᵢ₊₁| for all iterations i
    history_diff = abs(history[1 : end - 1, :] - history[2 : end, :]) + eps()

    p = Plots.plot()

    # Make sure we dont take log(0) anywhere
    for j = 1 : k
      Plots.plot!(filter(x -> x > 1e-12, history_diff[:, j]), yscale = :log10, mark = :auto, label = "$j")
    end

    p
  end

  # Shows that the init method works when
  # using orthogonal_iteration
  function ex1_init(n = 300, k = 3)
    # Get k eigenpairs
    starting_points = init(n, k)

    # Problem with quadratic rhs
    f, ∂f = quadratic(n)

    # Pick an initial μ and θ
    μ, θ = starting_points[3]

    # Apply continuation
    continuation!(f, ∂f, θ, linspace(μ, μ + 10, 100))

    # Plot the solution at μ + 10
    Plots.plot(θ)
    
  end

  function ex3(n = 300)
    f, ∂f = quadratic(n)
    
    # Nearly the trivial solution of the problem with the quadratic rhs
    θ = zeros(n)

    p = Plots.plot()

    # Let μ go from 30 to 0 in 1000 steps
    continuation_with_stability_check!(f, ∂f, θ, linspace(30, 0, 1000))
    Plots.plot!(θ)

    # Let μ got from 0 to 30 in 1000 steps
    continuation_with_stability_check!(f, ∂f, θ, linspace(0, 30, 1000))
    Plots.plot!(θ)
    
    # Flip the sign of θ
    θ *= -1;

    # Go from 30 to 20
    continuation_with_stability_check!(f, ∂f, θ, linspace(30, 20, 333))
    Plots.plot!(θ)

    # Go from 20 to 0.
    continuation_with_stability_check!(f, ∂f, θ, linspace(20, 0, 666))

    p
  end
end