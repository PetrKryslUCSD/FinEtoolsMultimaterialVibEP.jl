using Test

using FinEtools
using FinEtoolsMultimaterialVibEP: solve_ep

using FinEtools.MeshImportModule
using FinEtools.MeshExportModule: VTK
using MAT
using DataDrop
using LinearAlgebra

try rm("fin81mm.mat"); catch end
# # The name of the parameter file is up to you
parameterfile = "fin81mm.json"
solve_ep(parameterfile)
@test isfile("Finned_81mm_Max_12kHz.mat")

vars = matread("Finned_81mm_Max_12kHz.mat")
mX = vars["X"]
mconn = vars["conn"]
interiorxyz = vars["interiorxyz"]

vtkexportmesh("fin81mm-boundary.vtk", mconn, mX, VTK.T3; scalars = [("label", 1:size(mconn, 1))]);

VTK.vtkexportmesh("fin81mm-interiorxyz.vtk", reshape(collect(1:size(interiorxyz, 1)), size(interiorxyz, 1), 1), interiorxyz, VTK.P1)

