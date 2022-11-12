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
using Random
using Statistics
using Arpack
using MAT

function solve_ep(parameterfile)

    parameters = open(parameterfile, "r") do file
        JSON.parse(file)
    end
    @info "Loaded $parameterfile"

    materials = parameters["materials"]
    meshfile = parameters["meshfile"]
    neigvs = parameters["neigvs"]
    frequencyshift = parameters["frequencyshift"]
    maxiter = haskey(parameters, "maxiter") ? parameters["maxiter"] : 300
    method = haskey(parameters, "method") ? parameters["method"] : "eigs"
    withrr = haskey(parameters, "withrr") ? parameters["withrr"] : false
    tol = haskey(parameters, "tol") ? parameters["tol"] : 1.0e-3

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
        setlabel!(fesets[i], pids[i])
        E, nu, rho = materials[pid]["E"], materials[pid]["nu"], materials[pid]["rho"]
        material = MatDeforElastIso(MR, rho, E, nu, 0.0)
        femm = FEMMDeforLinearESNICET4(MR, IntegDomain(fesets[i], NodalSimplexRule(3)), material)
        femm = associategeometry!(femm, geom)
        K += stiffness(femm, geom, u)
        M += mass(femm, geom, u)
        if allfes == nothing
            allfes = deepcopy(fesets[i])
        else
            allfes = cat(allfes, deepcopy(fesets[i]))
        end
    end
    
    # eigs returns the nev requested eigenvalues in d, the corresponding Ritz vectors
    # v (only if ritzvec=true), the number of converged eigenvalues nconv, the number
    # of iterations niter and the number of matrix vector multiplications nmult, as
    # well as the final residual vector resid.
    @info "Solving eigenvalue problem for $neigvs frequencies"
    if method == "eigs"
        d, v, conv = eigs(Symmetric(K+OmegaShift*M), Symmetric(M); nev=neigvs, which=:SM, maxiter = maxiter, explicittransform=:none, check = 1)
        # @show "eigs", d
        d = d .- OmegaShift;
    else
        lamb, v, nconv, niter, lamberr = AlgoDeforLinearModule.ssit(K+OmegaShift*M, M; nev=neigvs, tol = tol, maxiter = maxiter, verbose=false)
        d = lamb
         # @show "ssit", d
        d = d .- OmegaShift;
        conv = nconv
    end
    @info "$conv eigenvalues converged"
    
    
    fs = real(sqrt.(complex(d)))/(2*pi)

    # size(d), size(v), size(M)
    # @show d
    @info "Checking orthogonality"
    tol = 1.0e-6
    max_vMv_diag_error = 0.0
    max_vMv_offdiag_error = 0.0
    for i in 1:length(d), j in 1:length(d)
        p = v[:, i]' * M * v[:, j]
        if i == j && abs(p - 1) > tol
            max_vMv_diag_error = max(max_vMv_diag_error, abs(p - 1))
        end
        if i != j && abs(p) > tol
            max_vMv_offdiag_error = max(max_vMv_offdiag_error, abs(p))
        end
    end
    max_vKv_diag_error = 0.0
    max_vKv_offdiag_error = 0.0
    for i in 1:length(d), j in 1:length(d)
        p = v[:, i]' * K * v[:, j]
        if i == j && abs(p - d[i]) > max(tol, tol*abs(d[i]))
            max_vKv_diag_error = max(max_vKv_diag_error, abs(p - d[i]))
        end
        if i != j && abs(p) > max(tol, tol*abs(d[i]))
            max_vKv_offdiag_error = max(max_vKv_offdiag_error, abs(p))
        end
    end
    @info "Mass: diagonal error = $(max_vMv_diag_error), off-diagonal error = $(max_vMv_offdiag_error) "
    @info "Stiffness: diagonal error = $(max_vKv_diag_error), off-diagonal error = $(max_vKv_offdiag_error) "
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
    bfes = outer_surface_of_solid(fens, bfes)
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
    npanels = count(sfes)
    # Compute the areas of all the boundary triangles. 
    areas = fill(0.0, npanels)
    for panel = 1:npanels
        femm1  =  FEMMBase(IntegDomain(subset(sfes, [panel]), SimplexRule(2, 1)))
        areas[panel] = integratefunction(femm1, sgeom, (x) ->  1.0)
    end

    numinteriorpoints_fraction = 0.05
    if "numinteriorpoints_fraction" in keys(parameters)
        numinteriorpoints_fraction = parameters["numinteriorpoints_fraction"]
    end
    numinteriorpoints = Int(round(numinteriorpoints_fraction * npanels))
    interiorxyz = interiorpoints(fens, allfes, numinteriorpoints; use_labels = true)

    file = matopen(meshfilebase * ".mat", "w")
    write(file, "Omega", d)
    write(file, "E", v)
    write(file, "X", mX)
    write(file, "conn", mConn)
    write(file, "G", G)
    write(file, "areas", areas)
    write(file, "interiorxyz", interiorxyz)
    close(file)

    true

end # solve_ep


# The interior points are added as randomly selected centroids of (a fraction
# of) the elements
function interiorpoints(fens, fes, numpoints; use_labels = false)::Matrix{Float64}
    xyz = fens.xyz
    # Fraction of finite elements that should have an interior point
    r = Float64(numpoints) / count(fes)
    if use_labels
        ulab  = unique(fes.label)
        # First count how many interior points we will generate
        numpoints = 0
        for l in ulab
            lst = selectelem(fens, fes, label = l)
            numpoints += max(1, Int(round(r * length(lst))))
        end
        numpoints = max(1, min(numpoints, count(fes)))
        interiorxyz = fill(zero(eltype(xyz)), numpoints, 3)
        p = 1
        for l in ulab
            lst = selectelem(fens, fes, label = l)
            sfes = subset(fes, lst)
            rix = randperm(count(sfes));
            conn = connasarray(sfes)
            np = max(1, Int(round(r * count(sfes))))
            for i in 1:np
                k = rix[i]
                interiorxyz[p, :] = mean(xyz[view(conn, k, :), :], dims=1)
                p += 1
            end
        end
        @assert p == numpoints + 1
    else
        conn = connasarray(fes)
        rix = randperm(count(fes));
        numpoints = min(numpoints, count(fes))
        interiorxyz = fill(zero(eltype(xyz)), numpoints, 3)
        for i in 1:numpoints
            k = rix[i]
            interiorxyz[i, :] = mean(xyz[view(conn, k, :), :], dims=1)
        end
    end
    
    return interiorxyz
end # function 



end # module FinEtoolsMultimaterialVibEP
