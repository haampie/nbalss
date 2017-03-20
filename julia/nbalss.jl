function newton(f::Function, ∂f::Function, x, ɛ = 1e-10)
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
  t = (n + 1)^2;
  return SymTridiagonal(-2 * t * ones(n), t * ones(n - 1));
end

function fBuck(θ::Vector{Float64}, μ)
  A = Poisson1D(length(θ));
  return A * θ + μ * sin(θ);
end

function ∂fBuck(θ::Vector{Float64}, μ)
  A = Poisson1D(length(θ));
  return A + Diagonal(μ * cos(θ));
end

function branch_off(k::Int64, n::Int64; γ = 0.1)
  # First k eigenvalues & eigenvectors
  eigen = eigfact(Poisson1D(n) * -1, 1 : k);
  correction = √((n + 1) / 2);
  ɛ = 8 * γ / (1 + γ);
  θs = ɛ * correction * eigen[:vectors];
  μs = (1 + γ) * eigen[:values];

  for idx = 1 : k
    (θs[:, idx], _) = newton(
      x -> fBuck(x, μs[idx]),
      x -> ∂fBuck(x, μs[idx]), 
      θs[:, idx]
    );
  end

  return (θs, μs);
end

function continue_on_branch(θ::Vector{Float64}, μ_start, μ_end, steps)
  for μ = linspace(μ_start, μ_end, steps)
    (θ, _) = newton(x -> fBuck(x, μ), x -> ∂fBuck(x, μ), θ);
  end

  return θ;
end

## Exercise 4.1
function ex4_1_two_norm(;k = 3, n = 1000, μ_end = 300.)
  (θs, μs) = branch_off(k, n);
  
  p = Plots.plot();

  for idx = 1 : k
    two_norms = Float64[];
    μ_range = μs[idx] : 0.4 : μ_end;
    θ = θs[:, idx];

    for μ = μ_range
      (θ, _) = newton(
        x -> fBuck(x, μ),
        x -> ∂fBuck(x, μ),
        θ
      );
      push!(two_norms, norm(θ) / √(n + 1));
    end

    Plots.plot!(μ_range, two_norms, label = "\$k = $idx\$");
  end

  return p;
end

function ex4_1_beam(;n = 1000, k = 3, μ_end = 200., μ_steps = 100)
  (θs, μs) = branch_off(k, n);

  p = Plots.plot();

  for idx = 1 : k
    θ = [0; continue_on_branch(θs[:, idx], μs[idx], μ_end, μ_steps); 0];
    x = cumsum(sin(θ)) / (n + 1);
    y = cumsum(cos(θ)) / (n + 1);
    Plots.plot!(x, y, label = "$idx");
  end

  return p;
end


## Exercise 4.3

function fBuckSinRHS(θ::Vector{Float64}, μ, ɛ, k)
  n = length(θ);
  h = 1 / (n + 1);
  rhs = ɛ * sin(k * pi * linspace(h, 1 - h, n));
  return fBuck(θ, μ) - rhs;
end

function ex4_3(; n = 1000, k = 1, ɛ = 0.1)
  θ = zeros(n);

  iterations = Int64[];
  μs = linspace(0, 50, 100);

  for μ = μs
    (θ, its) = newton(
      x -> fBuckSinRHS(x, μ, ɛ, k),
      x -> ∂fBuck(x, μ),
      θ
    );
    push!(iterations, length(its));
  end

  return (θ, iterations, μs);
end

## Exercise 4.4

function fBuckQuad(θ::Vector{Float64}, μ, ɛ)
  n = length(θ);
  h = 1 / (n + 1);
  s = linspace(h, 1 - h, n);
  rhs = ɛ * s .* (1 - s);
  return fBuck(θ, μ) - rhs;
end

function ex4_4a(; n = 1000, ɛ = 0.01)
  θ = zeros(n);
  iterations = Int64[];
  norms = Float64[];
  μs = linspace(0, 30, 1000);

  for μ = μs
    (θ, its) = newton(
      x -> fBuckQuad(x, μ, ɛ),
      x -> ∂fBuck(x, μ),
      θ
    );
    push!(iterations, length(its));
    push!(norms, norm(θ));
  end

  return (
    Plots.plot(μs, iterations, xlabel = "\$\\mu\$", label = ""),
    Plots.plot(μs, norms, xlabel = "\$\\mu\$", label = "")
  );
end

function ex4_4b(; n = 100, ɛ = 0.01)
  θ = zeros(n);
  iterations = Int64[];
  μ_switch = 20.;

  # Continue untill μ = μ_switch.
  for μ = linspace(0, μ_switch, 1000)
    (θ, _) = newton(
      x -> fBuckQuad(x, μ, ɛ),
      x -> ∂fBuck(x, μ),
      θ
    );
  end

  # Flip signs.
  θ_easy = copy(θ);
  θ *= -1;

  # Apply Newton.
  (θ, _) = newton(
    x -> fBuckQuad(x, μ_switch, ɛ), 
    x -> ∂fBuck(x, μ_switch),
    θ
  );

  # Return plots
  xs = linspace(0, 1, n + 2);
  Plots.plot(xs, [0; θ_easy; 0], label = "Easy solution");
  return Plots.plot!(xs, [0; θ; 0], label = "Branch switch");
end

function ex4_4c(; n = 100, ɛ = 0.01, μ_end = 20.)
  θ = zeros(n);

  # Continuation in μ
  for μ = linspace(0, μ_end, 1000)
    (θ, _) = newton(
      x -> fBuckQuad(x, μ, ɛ),
      x -> ∂fBuck(x, μ),
      θ
    );
  end

  θ_easy = copy(θ);

  # Continuation in ɛ
  for my_ɛ = linspace(ɛ, 0, 5)
    (θ, _) = newton(
      x -> fBuckQuad(x, μ_end, my_ɛ),
      x -> ∂fBuck(x, μ_end),
      θ
    );
  end

  # Flip signs.
  θ *= -1;

  # Continuation in ɛ
  for my_ɛ = linspace(0, ɛ, 5)
    (θ, _) = newton(
      x -> fBuckQuad(x, μ_end, my_ɛ),
      x -> ∂fBuck(x, μ_end),
      θ
    );
  end

  xs = linspace(0, 1, n + 2);
  Plots.plot(xs, [0; θ_easy; 0], label = "Easy solution");
  return Plots.plot!(xs, [0; θ; 0], label = "Branch switch");
end

function ex4_4d(; n = 100, ɛ = 0.01, μ_end = 20.)
  θ = zeros(n);

  # Continuation in μ
  for μ = linspace(0, μ_end, 1000)
    (θ, _) = newton(
      x -> fBuckQuad(x, μ, ɛ),
      x -> ∂fBuck(x, μ),
      θ
    );
  end

  θ_easy = copy(θ);

  # Flip signs.
  θ *= -1;

  r = fBuckQuad(θ, μ_end, ɛ);

  # Continuation in α
  for α = linspace(0, 1, 10)
    (θ, _) = newton(
      x -> fBuckQuad(x, μ_end, ɛ) - (1 - α) * r,
      x -> ∂fBuck(x, μ_end),
      θ
    );
  end

  xs = linspace(0, 1, n + 2);
  Plots.plot(xs, [0; θ_easy; 0], label = "Easy solution");
  return Plots.plot!(xs, [0; θ; 0], label = "Branch switch");
end
