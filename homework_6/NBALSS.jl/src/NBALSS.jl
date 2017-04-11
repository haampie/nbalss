module NBALSS

  export continuation, 
    newton, 
    nonlinear_backward_euler,
    grid_interior,
    poisson,
    problem,
    quadratic,
    orthogonal_iteration

  abstract Method
  immutable Invert <: Method end
  immutable Multiply <: Method end

  include("continuation.jl")
  include("nonlinear_backward_euler.jl")
  include("newton.jl")
  include("pde_and_grid.jl")
  include("orthogonal_iteration.jl")

end