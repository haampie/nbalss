module Exercises
  
  importall NBALSS
  
  import Plots, PGFPlots

  import Base: size

  import LinearMaps: LinearMap

  Plots.pgfplots()

  export cheby_test, A_mul_B!

  immutable InvertedMatrixPoly
    A
    Fs::Vector

    function InvertedMatrixPoly(A::AbstractMatrix, xs::AbstractVector)
      Fs = []
      for x = xs
        push!(Fs, factorize(A - UniformScaling(x)))
      end
      new(A, Fs)
    end
  end

  Base.size(S::InvertedMatrixPoly, n) = size(S.A, n)

  function Base.A_mul_B!(S::InvertedMatrixPoly, B)
    for F = S.Fs
      B .= F \ B
    end
    B
  end

  function Base.A_mul_B!(Y, S::InvertedMatrixPoly, B)
    copy!(Y, B) # Copy can be avoided in the first multiplication.
    A_mul_B!(S, Y)
    Y
  end

  function cheby_nodes(k::Int, a, b)
    (a + b) / 2.0 + (b - a) / 2 * cos((2 * (1 : k) - 1) * pi / (2k))
  end

  function cheby_test(n = 300, k = 5; max_iter = 15)
    A = -poisson(n)
    nodes = cheby_nodes(2, 0.0, 20.0)
    println(nodes)
    J = InvertedMatrixPoly(A, nodes)

    V, λ, history = orthogonal_iteration(J, k)

    history_diff = abs(history[1 : end - 1, :] - history[2 : end, :]) + eps()

    p = Plots.plot()

    # Make sure we dont take log(0) anywhere
    for j = 1 : k
      Plots.plot!(filter(x -> x > 1e-12, history_diff[:, j]), yscale = :log10, mark = :auto, label = "$j")
    end

    p

    # a = nodes[1]
    # b = nodes[2]
    # (b + a - √((b + a)^2 - 4*a*b + 4 ./ λ)) / 2
  end
end