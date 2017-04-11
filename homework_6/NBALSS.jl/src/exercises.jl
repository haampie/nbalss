module Exercises
  
  using NBALSS: poisson, orthogonal_iteration, Invert
  
  import Plots, PGFPlots

  Plots.pgfplots()

  export ex_1

  function ex_1()
    A = poisson(300)
    F = factorize(A)
    k = 3
    history = orthogonal_iteration(F, k, Invert())

    @show history

    history_diff = abs(history[1 : end - 1, :] - history[2 : end, :]) + eps()

    p = Plots.plot()

    for j = 1 : k
        Plots.plot!(filter(x -> x > 1e-12, history_diff[:, j]), yscale = :log10, mark = :auto, label = "$j")
    end

    p
  end
end