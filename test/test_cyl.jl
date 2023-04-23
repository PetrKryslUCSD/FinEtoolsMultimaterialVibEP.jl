using Test

using FinEtools
using FinEtoolsMultimaterialVibEP: solve_ep

using FinEtools.MeshImportModule
using FinEtools.MeshExportModule: VTK
using MAT

R = 1.0/2 * phun("ft")
L = 2 * phun("ft")
nR = 5
nL = 20

fens, fes = T4cylindern(R, L, nR, nL)
@show boundingbox(fens.xyz)

File = "cyl-volume.vtk"
VTK.vtkexportmesh(File, fens, fes)
let
    e = NASTRANExporter("cyl.nas")
    BEGIN_BULK(e)
    for i in eachindex(fens)
        GRID(e, i, fens.xyz[i,:])
    end
    j = 1
    for i in eachindex(fes)
        CTETRA(e, j, 1, collect(fes.conn[i]))
        j = j + 1
    end
    PSOLID(e, 1, 1)
    ENDDATA(e)
    close(e)
end

try rm("cyl.mat"); catch end
# # The name of the parameter file is up to you
parameterfile = "cyl-mode_based.json"
solve_ep(parameterfile)
@test isfile("cyl.mat")

vars = matread("cyl.mat")
mX = vars["X"]
mconn = vars["conn"]
interiorxyz = vars["interiorxyz"]

vtkexportmesh("cyl-boundary.vtk", mconn, mX, VTK.T3);

VTK.vtkexportmesh("cyl-interiorxyz.vtk", reshape(collect(1:size(interiorxyz, 1)), size(interiorxyz, 1), 1), interiorxyz, VTK.P1)


