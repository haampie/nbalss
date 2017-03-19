function newton(f::Function, ∂f::Function, x, ɛ::Float64)
  Δx_norm = Inf;
  history = Float64[];

  while Δx_norm > ɛ
    Δx = ∂f(x) \ f(x);
    x = x - Δx;
    Δx_norm = norm(Δx);
    push!(history, Δx_norm);
  end

  return (x, history);
end

function Poisson1D(n)
  return SymTridiagonal(-2 * ones(n), ones(n - 1));
end

function buckling(θ::Vector{Float64}, μ::Float64)
  n = length(θ);
  A = Poisson1D(n);
  return A * θ + μ * sin(θ) / (n + 1) ^ 2;
end

function ∂buckling(θ::Vector{Float64}, μ::Float64)
  n = length(θ);
  A = Poisson1D(n);
  return A + Diagonal(μ * cos(θ) / (n + 1) ^ 2);
end

function branch_off(k::Int64, n::Int64, tol::Float64 = 1e-8; γ::Float64 = 0.1)
  # First k eigenvalues & eigenvectors
  eigen = eigfact(Poisson1D(n) * -(n + 1) ^ 2, 1 : k);

  correction = √((n + 1) / 2);
  ɛ = 8 * γ / (1 + γ);
  θs = ɛ * correction * eigen[:vectors];
  μs = (1 + γ) * eigen[:values];

  for idx = 1 : k
    (θs[:, idx], _) = newton(x -> buckling(x, μs[idx]), x -> ∂buckling(x, μs[idx]), θs[:, idx], tol);
  end

  return (θs, μs);
end

function continue_on_branch(θ::Vector{Float64}, μ_start::Float64, μ_end::Float64, steps::Int64, tol::Float64 = 1e-8)
  for μ = linspace(μ_start, μ_end, steps)
    (θ, _) = newton(x -> buckling(x, μ), x -> ∂buckling(x, μ), θ, tol);
  end

  return θ;
end

function two_norm_of_solution(;k = 1, n = 200)
  (θs, bifurcation_points) = branch_off(k, n);
  θ = θs[:, k];
  
  μs = linspace(bifurcation_points[k], bifurcation_points[k] * 20, 600);
  two_norms = Float64[];

  for μ = μs
    (θ, _) = newton(x -> buckling(x, μ), x -> ∂buckling(x, μ), θ, 1e-8);
    push!(two_norms, norm(θ) / √(n + 1));
  end

  Plots.plot!(μs, two_norms, label = "$k");
end

function plot_solutions(;n::Int64 = 1000, k::Int64 = 6, μ_end::Float64 = 300., μ_steps::Int64 = 100)
  (θs, μs) = branch_off(k, n);

  p = Plots.plot();

  for idx = 1 : k
    θ_end = [0; continue_on_branch(θs[:, idx], μs[idx], μ_end, μ_steps); 0];
    x = cumsum(sin(θ_end)) / (n + 1);
    y = cumsum(cos(θ_end)) / (n + 1);
    Plots.plot!(x, y, label = "$idx");
  end

  return p;
end


## Exercise 4.3

function buckling_plus_ɛ(θ::Vector{Float64}, μ::Float64, ɛ::Float64, k::Int64)
  n = length(θ);
  h = 1 / (n + 1);
  A = Poisson1D(n);
  rhs = ɛ * sin(k * pi * linspace(h, 1 - h, n));
  return A * θ + (μ * sin(θ) - rhs) * h ^ 2;
end

function ex4_3(; n::Int64 = 1000, k::Int64 = 1, ɛ::Float64 = 0.1, tol::Float64 = 1e-8)
  θ = zeros(n);

  iterations = Int64[];
  μs = linspace(0, 50, 100);

  for μ = μs
    (θ, its) = newton(x -> buckling_plus_ɛ(x, μ, ɛ, k), x -> ∂buckling(x, μ), θ, tol);
    push!(iterations, length(its));
  end

  return (θ, iterations, μs);
end

## Exercise 4.4

function buckling_plus_poly(θ::Vector{Float64}, μ::Float64, ɛ::Float64, k::Int64)
  n = length(θ);
  h = 1 / (n + 1);
  A = Poisson1D(n);
  s = linspace(h, 1 - h, n);
  rhs = ɛ * s .* (1 - s);
  return A * θ + (μ * sin(θ) - rhs) * h ^ 2;
end

function ex4_4(; n::Int64 = 1000, k::Int64 = 1, ɛ::Float64 = 0.01, tol::Float64 = 1e-8)
  θ = zeros(n);

  iterations = Int64[];
  μs = linspace(0, 50, 50);

  for μ = μs
    (θ, its) = newton(x -> buckling_plus_poly(x, μ, ɛ, k), x -> ∂buckling(x, μ), θ, tol);
    push!(iterations, length(its));
  end

  return (θ, iterations, μs);
end