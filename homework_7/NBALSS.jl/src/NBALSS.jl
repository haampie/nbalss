module NBALSS

  export
    grid_interior,
    poisson,
    problem,
    quadratic,
    orthogonal_iteration

  include("pde_and_grid.jl")
  include("orthogonal_iteration.jl")

end