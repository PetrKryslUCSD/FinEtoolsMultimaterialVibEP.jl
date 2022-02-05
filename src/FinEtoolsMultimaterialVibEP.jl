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
using Infiltrator

function solve_ep(parameterfile)

    parameters = open(parameterfile, "r") do file
        JSON.parse(file)
    end
    @info "Loaded $parameterfile"

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

    # Check that all property IDs are defined in the input file
    for p in pids
        if !(string(p) in keys(materials))
            @error "Material $(string(p)) was not defined in $(parameterfile)"
        end
    end

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
    @info "Solving eigenvalue problem for $neigvs frequencies"
    d, v, conv = eigs(K+OmegaShift*M, M; nev=neigvs, which=:SM, explicittransform=:none)
    @assert  conv == length(d)
    d = d .- OmegaShift;
    fs = real(sqrt.(complex(d)))/(2*pi)

    # @show size(d), size(v), size(M)
    # @show d
    # @info "Checking orthogonality"
    # tol = 1.0e-6
    # for i in 1:length(d), j in 1:length(d)
    #     p = v[:, i]' * M * v[:, j]
    #     if i == j && abs(p - 1) > tol
    #         @show i, p
    #     end
    #     if i != j && abs(p) > tol
    #         @show i, j, p
    #     end
    # end
    # for i in 1:length(d), j in 1:length(d)
    #     p = v[:, i]' * K * v[:, j]
    #     if i == j && abs(p - d[i]) > max(tol, tol*abs(d[i]))
    #         @show i, p
    #     end
    #     if i != j && abs(p) > max(tol, tol*abs(d[i]))
    #         @show i, j, p, d[i]
    #     end
    # end
    println("Eigenvalues: $fs [Hz]")
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

    # Define an auxiliary field: it has 1 degree of freedom at all the nodes of the mesh
    P = NodalField(collect(1:count(fens)));
    numberdofs!(P);

    # Find the nodes on the surface
    sn = connectednodes(bfes);
    on = setdiff(1:count(fens), sn);
    # Put a constraint "is on surface" == 0 (false) at the interior nodes
    setebc!(P, on, true, 1, 0);
    applyebc!(P);
    # Only the surface nodes will carry a nonzero equation number
    permutation = zeros(Int, count(fens));
    permutation[sn] = 1:length(sn);
    permutation[on] = length(sn)+1:length(permutation);
    numberdofs!(P, permutation);

    mX = fill(Inf, P.nfreedofs, 3);
    for j in 1:count(fens)
        if P.dofnums[j] > 0
            mX[P.dofnums[j], :] = fens.xyz[j, :]
        end
    end
    # Create a machine for the surface integrals, and calculate the coupling
    # structural-acoustic matrix. Note that the acoustic fluid properties are
    # not used in this method and they are supplied as dummy values.
    bfemm = FEMMAcoustSurf(IntegDomain(bfes, SimplexRule(2, 1)), MatAcoustFluid(1.0, 1.0))
    G = acousticcouplingpanels(bfemm, geom, u);

    # Construct the matrix of the connectivities of all the surface finite elements (panels)
    mConn = zeros(Int, count(bfes), nodesperelem(bfes))
    for j in 1:size(mConn, 1)
        mConn[j,:] .= P.dofnums[bconn[j, :]]; # need to renumber to match mX
    end

    # @infiltrate

    # MESH.write_MESH("$(meshfilebase)-surface.mesh", fens, bfes) 

    sgeom = NodalField(mX)
    sfes = fromarray!(bfes, mConn)
    # Compute the areas of all the boundary triangles. 
    areas = fill(0.0, count(sfes))
    for panel = 1:count(sfes)
        femm1  =  FEMMBase(IntegDomain(subset(sfes, [panel]), SimplexRule(2, 1)))
        areas[panel] = integratefunction(femm1, sgeom, (x) ->  1.0)
    end

    file = matopen(meshfilebase * ".mat", "w")
    write(file, "Omega", d)
    write(file, "E", v)
    write(file, "X", mX)
    write(file, "conn", mConn)
    write(file, "G", G)
    write(file, "areas", areas)
    close(file)

    true

end # solve_ep

end # module FinEtoolsMultimaterialVibEP
