using Test

using FinEtoolsMultimaterialVibEP: solve_ep

try rm("twoblocks.mat"); catch end
# The name of the parameter file is up to you
parameterfile = joinpath(@__DIR__, "twoblocks.json")
solve_ep(parameterfile)
@test isfile("twoblocks.mat")

try rm("twoblocks.mat"); catch end
# The name of the parameter file is up to you
parameterfile = joinpath(@__DIR__, "twoblocks_ssit.json")
solve_ep(parameterfile)
@test isfile("twoblocks.mat")
