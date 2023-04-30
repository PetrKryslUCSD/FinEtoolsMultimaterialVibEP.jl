using Test

using FinEtools
using FinEtoolsMultimaterialVibEP: solve_ep

using FinEtools.MeshImportModule
using FinEtools.MeshExportModule: VTK
using MAT
using DataDrop
using LinearAlgebra
using PlotlyLight

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


G = DataDrop.retrieve_matrix("sim~modal~Finned_81mm_Max_12kHz_nas~50_50_20000~ne=200-no-chief-7_json~G.mdat");
@show size(G)

mG = vars["G"]
@show size(mG)

@show G - mG


Phi = DataDrop.retrieve_matrix("sim~modal~Finned_81mm_Max_12kHz_nas~50_50_20000~ne=200-no-chief-7_json~Phi.mdat");
@show size(Phi)


mE = vars["E"]
@show size(mE)

Phi = abs.(Phi)
mE = abs.(mE)
@show norm(Phi[:, 7:end] - mE[:, 7:end]) / norm(mE[:, 7:end])

p = PlotlyLight.Plot()
j = 7
@show norm(Phi[:, j] - mE[:, j]) / norm(mE[:, j])
p(x = 1:size(Phi, 1), y = Phi[:, j])
p(x = 1:size(mE, 1), y = mE[:, j])
display(p)

mOmega = vars["Omega"];
@show size(mOmega)


Omega = DataDrop.retrieve_matrix("sim~modal~Finned_81mm_Max_12kHz_nas~50_50_20000~ne=200-no-chief-7_json~Omega.mdat");

@show norm(Omega - mOmega) / norm(Omega)


bconn = DataDrop.retrieve_matrix("sim~modal~Finned_81mm_Max_12kHz_nas~50_50_20000~ne=200-no-chief-7_json~bconn.mdat");
xyz = DataDrop.retrieve_matrix("sim~modal~Finned_81mm_Max_12kHz_nas~50_50_20000~ne=200-no-chief-7_json~xyz.mdat");

vtkexportmesh("fin81mm-boundary-2.vtk", bconn, xyz, VTK.T3; scalars = [("label", 1:size(bconn, 1))]);


areas = DataDrop.retrieve_matrix("sim~modal~Finned_81mm_Max_12kHz_nas~50_50_20000~ne=200-no-chief-7_json~areas.mdat");

@show norm(areas - vars["areas"]) / norm(areas)


xyz = DataDrop.retrieve_matrix("sim~modal~Finned_81mm_Max_12kHz_nas~50_50_20000~ne=200-no-chief-7_json~xyz.mdat");

@show norm(xyz - vars["X"]) / norm(xyz)


bconn = DataDrop.retrieve_matrix("sim~modal~Finned_81mm_Max_12kHz_nas~50_50_20000~ne=200-no-chief-7_json~bconn.mdat");

@show norm(bconn - vars["conn"]) / norm(bconn)

