__precompile__(true)

module LevySupSim

  using Distributions
  
  export
    # methods
    rand_levy_sup,
    VG,
    NIG,
    WS

  # Source Files
  include("levysupsim.jl")

end
