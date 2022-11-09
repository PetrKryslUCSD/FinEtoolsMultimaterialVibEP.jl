using Test

using FinEtools
using FinEtoolsMultimaterialVibEP: solve_ep

using FinEtools.MeshImportModule
using FinEtools.MeshExportModule: VTK
using MAT

output = MeshImportModule.import_ABAQUS("./two-mat-cyl-w-voids.inp"; allocationchunk = 11)

fens, fes = output["fens"], output["fesets"][1]

fens, fes = T10toT4(fens, fes)

@test length(output["elsets"]["SET-1"]) == 5737
@test length(output["elsets"]["SET-2"])  == 5774
File = "cyl-w-voids-two-mats-elset-1.vtk"
VTK.vtkexportmesh(File, fens, subset(fes, output["elsets"]["SET-1"]))
File = "cyl-w-voids-two-mats-elset-2.vtk"
VTK.vtkexportmesh(File, fens, subset(fes, output["elsets"]["SET-2"]))

let
    e = NASTRANExporter("two-mat-cyl-w-voids.nas")
    BEGIN_BULK(e)
    for i in 1:count(fens)
        GRID(e, i, fens.xyz[i,:])
    end
    j = 1
    sfes = subset(fes, output["elsets"]["SET-1"])
    conn = connasarray(sfes)
    for i in 1:count(sfes)
        CTETRA(e, j, 1, conn[i, :])
        j = j + 1
    end
    sfes = subset(fes, output["elsets"]["SET-2"])
    conn = connasarray(sfes)
    for i in 1:count(sfes)
        CTETRA(e, j, 2, conn[i, :])
        j = j + 1
    end
    PSOLID(e, 1, 1)
    PSOLID(e, 2, 2)
    ENDDATA(e)
    close(e)
end

# # The name of the parameter file is up to you
parameterfile = "two-mat-cyl-w-voids.json"
solve_ep(parameterfile)

vars = matread("two-mat-cyl-w-voids.mat")
mX = vars["X"]
mconn = vars["conn"]
vtkexportmesh("two-mat-cyl-w-voids-boundary.vtk", mconn, mX, VTK.T3);

# @test isfile("twoblocks.mat")