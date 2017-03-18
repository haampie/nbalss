function newton(f::Function, ∂f::Function, x, ɛ::Float64, args...)

  Δx_norm = Inf;
  history = Float64[];

  while Δx_norm > ɛ
    Δx = ∂f(x, args...) \ f(x, args...);
    x = x - Δx;
    Δx_norm = norm(Δx);
    push!(history, Δx_norm);
  end

  return (x, history);
end

function poisson1d(n)
  return SymTridiagonal(-2 * ones(n), ones(n - 1));
end

function buckling(θ::Vector{Float64}, μ::Float64)
  n = length(θ);
  A = poisson1d(n);
  return A * θ + μ * sin(θ) / (n + 1) ^ 2;
end

function ∂buckling(θ::Vector{Float64}, μ::Float64)
  n = length(θ);
  A = poisson1d(n);
  return A + Diagonal(μ * cos(θ) / (n + 1) ^ 2);
end

function branch_off(k::Int64, n::Int64, tol::Float64 = 1e-8)
  # Last k eigenvalues & eigenvectors
  eigen = eigfact(poisson1d(n) * -(n + 1) ^ 2, 1 : k);

  # Reverse order of eigen pairs (from small magnitude to large)
  θs = eigen[:vectors];
  μs = eigen[:values];

  correction = √((n + 1) / 2);
  γ = 0.1;
  ɛ = 8 * γ / (1 + γ);

  θs *= ɛ * correction;
  μs *= 1 + γ;

  for idx = 1 : k
    (θs[:, idx], _) = newton(buckling, ∂buckling, θs[:, idx], tol, μs[idx]);
  end

  return (θs, μs);
end

function step_along_branch(θ::Vector{Float64}, μ_start::Float64, μ_end::Float64, steps::Int64, tol::Float64 = 1e-8)
  for μ = linspace(μ_start, μ_end, steps)
    (θ, _) = newton(buckling, ∂buckling, θ, tol, μ);
  end

  return θ;
end

function run_example()
  idx = 5;
  k = 6;
  n = 200;
  (θs, μs) = branch_off(k, n);

  my_θ = θs[:, idx];

  Plots.plot();

  anim = @animate for μ = linspace(μs[idx], μs[idx] + 100, 20)
    (my_θ, _) = newton(buckling, ∂buckling, my_θ, 1e-8, μ);
    x = cumsum(sin([0;my_θ;0])) / (n + 1);
    y = cumsum(cos([0;my_θ;0])) / (n + 1);
    Plots.plot(x, y, title = "$μ");
  end

  gif(anim, "./anim_fps15.gif", fps = 15);
end