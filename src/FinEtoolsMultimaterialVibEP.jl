module FinEtoolsMultimaterialVibEP

using JSON
using Printf
using DelimitedFiles
using FinEtools
using FinEtools.MeshExportModule: MESH
using FinEtoolsDeforLinear.AlgoDeforLinearModule
using FinEtoolsDeforLinear
using FinEtoolsAcoustics
using FinEtools.MeshExportModule
using LinearAlgebra
using SparseArrays
using Arpack
using MAT

function solve_ep(parameterfile)

    parameters = open(parameterfile, "r") do file
        JSON.parse(file)
    end

    materials = parameters["materials"]
    meshfile = parameters["meshfile"]
    neigvs = parameters["neigvs"]
    frequencyshift = parameters["frequencyshift"]

    meshfilebase, ext = splitext(meshfile)

    OmegaShift = (2*pi*frequencyshift) ^ 2; # to resolve rigid body modes
    
    # No need to change anything below this line ##########
    MR = DeforModelRed3D
    output = import_NASTRAN(meshfile)
    fens, fesets, pids = output["fens"], output["fesets"], output["property_ids"] 

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field

    numberdofs!(u)

    K = spzeros(u.nfreedofs, u.nfreedofs)
    M = spzeros(u.nfreedofs, u.nfreedofs)
    allfes = nothing
    for i in 1:length(fesets)
        pid = string(pids[i])
        E, nu, rho = materials[pid]["E"], materials[pid]["nu"], materials[pid]["rho"]
        material = MatDeforElastIso(MR, rho, E, nu, 0.0)
        femm = FEMMDeforLinearESNICET4(MR, IntegDomain(fesets[i], NodalSimplexRule(3)), material)
        femm = associategeometry!(femm, geom)
        K += stiffness(femm, geom, u)
        M += mass(femm, geom, u)
        if allfes == nothing
            allfes = fesets[i]
        else
            allfes = cat(allfes, fesets[i])
        end
    end
    
    # eigs returns the nev requested eigenvalues in d, the corresponding Ritz vectors
    # v (only if ritzvec=true), the number of converged eigenvalues nconv, the number
    # of iterations niter and the number of matrix vector multiplications nmult, as
    # well as the final residual vector resid.

    d,v,nev,nconv = eigs(K+OmegaShift*M, M; nev=neigvs, which=:SM)
    d = d .- OmegaShift;
    fs = real(sqrt.(complex(d)))/(2*pi)
    # println("Eigenvalues: $fs [Hz]")
    # open(meshfilebase * "-eval" * ".mat", "w") do file
    #     writedlm(file, d)
    #     # write(file, matrix[:])
    # end
    # mode = 7
    # scattersysvec!(u, v[:,mode])
    # File =  "multimaterial_nas.vtk"
    # vtkexportmesh(File, fens, allfes; vectors=[("mode$mode", u.values)])
    # @async run(`"paraview.exe" $File`)
    # for i in 1:length(d)
    #     open(meshfilebase * "-evec$(i)" * ".mat", "w") do file
    #         writedlm(file, v[:, i])
    #         # write(file, matrix[:])
    #     end
    # end

    # Extract the boundary
    bfes = meshboundary(allfes)
    bconn = connasarray(bfes)

    # MESH.write_MESH("$(meshfilebase)-surface.mesh", fens, bfes) 

    # Create a machine for the surface integrals, and calculate the coupling structural-acoustic matrix. Note that the acoustic fluid properties are not used in this method and they are supplied as dummy values.
    bfemm = FEMMAcoustSurf(IntegDomain(bfes, SimplexRule(2, 1)), MatAcoustFluid(1.0, 1.0))
    G = acousticcouplingpanels(bfemm, geom, u);

    # Compute the areas of all the boundary triangles. 
    areas = fill(0.0, count(bfes))
    for panel = 1:count(bfes)
        femm1  =  FEMMBase(IntegDomain(subset(bfes, [panel]), SimplexRule(2, 1)))
        areas[panel] = integratefunction(femm1, geom, (x) ->  1.0)
    end

    # open(meshfilebase * "-areas" * ".mat", "w") do file
    #     writedlm(file, areas)
    #     # write(file, matrix[:])
    # end

    I, J, V =  SparseArrays.findnz(G)
    # open(meshfilebase * "-G" * ".mat", "w") do file
    #     for i in 1:length(I)
    #         @printf(file, "%d, %d, %g\n", I[i], J[i], V[i])
    #     end
    #     # writedlm(file, hcat(I, J, V))
    #     # write(file, matrix[:])
    # end

    Omega = d
    E = v
    X = fens.xyz
    conn = bconn
    file = matopen(meshfilebase * ".mat", "w")
    write(file, "Omega", d)
    write(file, "E", v)
    write(file, "X", fens.xyz)
    write(file, "conn", bconn)
    write(file, "G", G)
    close(file)

    true

end # solve_ep

end # module FinEtoolsMultimaterialVibEP
