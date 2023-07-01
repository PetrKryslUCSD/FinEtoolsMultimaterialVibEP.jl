module FinEtoolsMultimaterialVibEP

using JSON
using Printf
using DelimitedFiles
using FinEtools
using FinEtools.MeshExportModule: MESH
using SubSIt: ssit
using FinEtoolsDeforLinear
using FinEtoolsAcoustics
using FinEtools.MeshExportModule
using LinearAlgebra
using SparseArrays
using Random
using Statistics
using Arpack
using MAT
using DataDrop
using Statistics
using FinEtoolsTetsFromTris

function solve_ep(parameterfile)

    Random.seed!(1234);

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
    verbose = haskey(parameters, "verbose") ? parameters["verbose"] : false
    tol = haskey(parameters, "tol") ? parameters["tol"] : 1.0e-3
    check_orthogonality = haskey(parameters, "check_orthogonality") ? parameters["check_orthogonality"] : true

    meshfilebase, _ = splitext(meshfile)

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
    @info "$(count(fens)) nodes"

    K = spzeros(nalldofs(u), nalldofs(u))
    M = spzeros(nalldofs(u), nalldofs(u))
    allfes = nothing
    for i in 1:length(fesets)
        pid = string(pids[i])
        setlabel!(fesets[i], pids[i])
        E, nu, rho = Float64(materials[pid]["E"]), Float64(materials[pid]["nu"]), Float64(materials[pid]["rho"])
        material = MatDeforElastIso(MR, rho, E, nu, 0.0)
        femm = FEMMDeforLinearESNICET4(MR, IntegDomain(fesets[i], NodalSimplexRule(3)), material)
        associategeometry!(femm, geom; stabilization_parameters=(2.0, 3.0))
        K += stiffness(femm, geom, u)
        M += mass(femm, geom, u)
        if allfes === nothing
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
        d, v, nconv = eigs(Symmetric(K+OmegaShift*M), Symmetric(M); nev=neigvs, which=:SM, maxiter = maxiter, explicittransform=:none, check = 1)
        d = d .- OmegaShift;
    else
        d, v, nconv, niter, lamberr = ssit(K+OmegaShift*M, M; nev=neigvs, tol = tol, maxiter = maxiter, verbose=verbose)
        d = d .- OmegaShift;
    end
    @info "$nconv eigenvalues converged"

    Omega = deepcopy(d)
    Phi = deepcopy(v)

    fs = real(sqrt.(complex(d)))/(2*pi)

    if check_orthogonality
        @info "Checking orthogonality"
        tol = 1.0e-6
        max_vMv_diag_error = 0.0
        max_vMv_offdiag_error = 0.0
        Mred = v' * M * v
        for i in 1:length(d), j in i:length(d)
            p = Mred[i, j]
            if i == j && abs(p - 1) > tol
                max_vMv_diag_error = max(max_vMv_diag_error, abs(p - 1))
            end
            if i != j && abs(p) > tol
                max_vMv_offdiag_error = max(max_vMv_offdiag_error, abs(p))
            end
        end
        max_vKv_diag_error = 0.0
        max_vKv_offdiag_error = 0.0
        Kred = v' * K * v
        for i in 1:length(d), j in i:length(d)
            p = Kred[i, j]
            if i == j && abs(p - d[i]) > max(tol, tol*abs(d[i]))
                max_vKv_diag_error = max(max_vKv_diag_error, abs(p - d[i]))
            end
            if i != j && abs(p) > max(tol, tol*abs(d[i]))
                max_vKv_offdiag_error = max(max_vKv_offdiag_error, abs(p))
            end
        end
        @info "Mass: diagonal error = $(max_vMv_diag_error), off-diagonal error = $(max_vMv_offdiag_error) "
        @info "Stiffness: diagonal error = $(max_vKv_diag_error), off-diagonal error = $(max_vKv_offdiag_error) "
    end

    @info("Natural frequencies:\n $(round.(fs, digits=5)) [Hz]")

    # File =  "modes.vtk"
    # vectors = []
    # for mode  in  1:20
    #     scattersysvec!(u, v[:,mode])
    #     push!(vectors, ("mode-$(mode)", deepcopy(u.values)))
    # end
    # vtkexportmesh(File, fens, allfes; vectors=vectors)

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

    # mX = fill(Inf, P.nfreedofs, 3);
    # for j in 1:count(fens)
    #     if P.dofnums[j] > 0
    #         mX[P.dofnums[j], :] = fens.xyz[j, :]
    #     end
    # end
    # Create a machine for the surface integrals, and calculate the coupling
    # structural-acoustic matrix. Note that the acoustic fluid properties are
    # not used in this method and they are supplied as dummy values.
    bfemm = FEMMAcoustSurf(IntegDomain(bfes, SimplexRule(2, 1)), MatAcoustFluid(1.0, 1.0))
    G = acousticcouplingpanels(bfemm, geom, u);

    # Construct the matrix of the connectivities of all the surface finite elements (panels)
    # mConn = zeros(Int, count(bfes), nodesperelem(bfes))
    # for j in 1:size(mConn, 1)
        # mConn[j,:] .= P.dofnums[bconn[j, :]]; # need to renumber to match mX
    # end

    # @infiltrate

    # MESH.write_MESH("$(meshfilebase)-surface.mesh", fens, bfes) 

    # sgeom = NodalField(mX)
    # sfes = fromarray!(bfes, mConn)
    sfes = bfes
    npanels = count(sfes)
    # Compute the areas of all the boundary triangles. 
    areas = fill(0.0, npanels)
    for panel = 1:npanels
        femm1  =  FEMMBase(IntegDomain(subset(sfes, [panel]), SimplexRule(2, 1)))
        areas[panel] = integratefunction(femm1, geom, (x) ->  1.0)
    end

    numinteriorpoints_fraction = 0.05
    if "numinteriorpoints_fraction" in keys(parameters)
        numinteriorpoints_fraction = parameters["numinteriorpoints_fraction"]
    end
    if "numinteriorpoints" in keys(parameters)
        numinteriorpoints = parameters["numinteriorpoints"]
    else
        numinteriorpoints = Int(round(numinteriorpoints_fraction * npanels))
    end
    interiorpoint_method = "random"
    if "interiorpoint_method" in keys(parameters)
        interiorpoint_method = parameters["interiorpoint_method"]
    end
    nummodes = 1
    single_point_per_mode = false
    rhow = 1000.0
    cw = 1500.0
    if interiorpoint_method == "mode_based"
        if "single_point_per_mode" in keys(parameters)
            single_point_per_mode = parameters["single_point_per_mode"]
        end
        if "nummodes" in keys(parameters)
            nummodes = parameters["nummodes"]
        else
            nummodes = (single_point_per_mode ? numinteriorpoints : Int(round(numinteriorpoints / 5)))
        end
        if "rhow" in keys(parameters)
            rhow = parameters["rhow"]
        end
        if "cw" in keys(parameters)
            cw = parameters["cw"]
        end
    end
    interiorxyz = interiorpoints(fens, allfes; interiorpoint_method, numinteriorpoints, nummodes, single_point_per_mode, rhow, cw)

    file = matopen(meshfilebase * ".mat", "w")
    write(file, "Omega", Omega)
    write(file, "E", Phi)
    write(file, "X", fens.xyz)
    write(file, "conn", bconn)
    write(file, "G", G)
    write(file, "areas", areas)
    write(file, "interiorxyz", interiorxyz)
    close(file)

    f = DataDrop.with_extension(DataDrop.clean_file_name(meshfilebase * "-model"), ".h5")
    DataDrop.empty_hdf5_file(f)
    DataDrop.store_matrix(f, "/xyz", fens.xyz)
    DataDrop.store_matrix(f, "/dofnums", u.dofnums)
    DataDrop.store_matrix(f, "/conn", connasarray(allfes))

    true
end # solve_ep


# # The interior points are added as randomly selected centroids of (a fraction
# # of) the elements
# function interiorpoints(fens, fes, numpoints; use_labels = false)::Matrix{Float64}
#     xyz = fens.xyz
#     # Fraction of finite elements that should have an interior point
#     r = Float64(numpoints) / count(fes)
#     if use_labels
#         ulab  = unique(fes.label)
#         # First count how many interior points we will generate
#         numpoints = 0
#         for l in ulab
#             lst = selectelem(fens, fes, label = l)
#             numpoints += max(1, Int(round(r * length(lst))))
#         end
#         numpoints = max(1, min(numpoints, count(fes)))
#         interiorxyz = fill(zero(eltype(xyz)), numpoints, 3)
#         p = 1
#         for l in ulab
#             lst = selectelem(fens, fes, label = l)
#             sfes = subset(fes, lst)
#             rix = randperm(count(sfes));
#             conn = connasarray(sfes)
#             np = max(1, Int(round(r * count(sfes))))
#             for i in 1:np
#                 k = rix[i]
#                 interiorxyz[p, :] = mean(xyz[view(conn, k, :), :], dims=1)
#                 p += 1
#             end
#         end
#         @assert p == numpoints + 1
#     else
#         conn = connasarray(fes)
#         rix = randperm(count(fes));
#         numpoints = min(numpoints, count(fes))
#         interiorxyz = fill(zero(eltype(xyz)), numpoints, 3)
#         for i in 1:numpoints
#             k = rix[i]
#             interiorxyz[i, :] = mean(xyz[view(conn, k, :), :], dims=1)
#         end
#     end

#     return interiorxyz
# end # function

function _random_points(fens, fes, numpoints)
    xyz = fens.xyz
    conn = connasarray(fes)
    rix = randperm(count(fes));
    interiorxyz = fill(0.0, numpoints, 3)
    for i in 1:numpoints
        k = rix[i]
        interiorxyz[i, :] = mean(xyz[view(conn, k, :), :], dims=1)
    end
    return interiorxyz
end


function _do_mode_based_points(fens, fes, nummodes, rho, c, compute_threshold)
    bulk =  c^2*rho;

    geom = NodalField(fens.xyz)
    P = NodalField(zeros(size(fens.xyz,1),1))
    bfes = meshboundary(fes)
    setebc!(P, connectednodes(bfes))
    numberdofs!(P)

    femm = FEMMAcoust(IntegDomain(fes, TetRule(1)),
        MatAcoustFluid(bulk,rho))

    S = acousticstiffness(femm, geom, P);
    C = acousticmass(femm, geom, P);

    d, v, _ = eigs(C, S; nev=nummodes, which=:SM, explicittransform=:none)
    v = real.(v)
    fs=real(sqrt.(complex(d)))./(2*pi)
    @info("Interior frequencies: $fs [Hz]")
    # ks = (2*pi).*fs./c./phun("m")
    # println("Wavenumbers: $(ks) [m]")

    magP = deepcopy(P)
    interiorxyz = fill(0.0, 0, 3)
    for s in 1:nummodes
        scattersysvec!(magP, abs.(v[:, s])) # rectified
        thresholdP = compute_threshold(magP)
        x = _identify_peaks(fens, fes, magP, thresholdP, s)
        interiorxyz = vcat(interiorxyz, x)
    end

    return interiorxyz
end # _do_mode_based_points

function _mode_based_points(fens, fes, nummodes, rho, c, compute_threshold)
    bfes = meshboundary(fes)
    obfes = outer_surface_of_solid(fens, bfes)
    if count(bfes) == count(obfes)
        # the outer boundary is all the boundary there is
        return _do_mode_based_points(fens, fes, nummodes, rho, c, compute_threshold)
    else
        # the outer boundary needs to be filled in with tetrahedral mesh to
        # represent the complement of the fluid domain
        _fens = deepcopy(fens)
        connected = findunconnnodes(_fens, obfes);
        _fens, new_numbering = compactnodes(_fens, connected);
        bfes = renumberconn!(obfes, new_numbering);
        __fens, __fes = FinEtoolsTetsFromTris.mesh(_fens, obfes)
        return _do_mode_based_points(__fens, __fes, nummodes, rho, c, compute_threshold)
    end
end # _mode_based_points

function _first_included(inclusion_list)
    for j in eachindex(inclusion_list)
        if inclusion_list[j]
            return j
        end
    end
    return 0
end

function _identify_peaks(fens, fes, magP, thresholdP, s)
    # First prune all elements whose node(s)'s value lies below a threshold
    # vtkexportmesh("all.vtk", connasarray(fes), fens.xyz, FinEtools.MeshExportModule.VTK.T4)
    inclusion_list = fill(false, count(fes))
    for i in 1:count(fes)
        isin = true
        for k in fes.conn[i]
            if magP.values[k] < thresholdP
                isin = false; break
            end
        end
        inclusion_list[i] = isin
    end
    pv = fill(0.0, length(magP.values))
    peak_locations = []
    # Now identify bunches of connected elements, and from each connected bunch
    # pick the node with the largest absolute value of the pressure
    while true
        f = _first_included(inclusion_list)
        if f == 0
            break
        end
        il = findall(inclusion_list)
        sfes = subset(fes, il)
        # vtkexportmesh("il-$s-$f.vtk", connasarray(sfes), fens.xyz, FinEtools.MeshExportModule.VTK.T4)
        el = selectelem(fens, sfes, flood = true, startnode = fes.conn[f][1])
        @assert !isempty(el)
        # vtkexportmesh("el-$s-$f.vtk", connasarray(subset(sfes, el)), fens.xyz, FinEtools.MeshExportModule.VTK.T4)
        # @show length(el)
        inclusion_list[il[el]] .= false
        nl = connectednodes(subset(sfes, el))
        pv .= 0
        for k in nl
            pv[k] = magP.values[k]
        end
        npeak = argmax(pv)
        push!(peak_locations, fens.xyz[npeak, :])
    end
    # Now we have the peak locations
    interiorxyz = fill(0.0, length(peak_locations), 3)
    for i in eachindex(peak_locations)
        interiorxyz[i, :] .= peak_locations[i]
    end
    return interiorxyz
end

function interiorpoints(fens, fes; interiorpoint_method = "random", numinteriorpoints = 10, nummodes = 1, single_point_per_mode = false, rhow = 1000.0, cw = 1500.0)::Matrix{Float64}
    compute_threshold(magP) = (single_point_per_mode ? 0.0
        : min(maximum(magP.values) * 0.8, mean(magP.values)  * 2)
    )

    if interiorpoint_method == "mode_based"
        interiorxyz = _mode_based_points(fens, fes, nummodes, rhow, cw, compute_threshold)
    elseif interiorpoint_method == "random"
        interiorxyz = _random_points(fens, fes, numinteriorpoints)
    else
        error("Unknown interior-point method")
        interiorxyz = nothing
    end

    return interiorxyz
end # function



end # module FinEtoolsMultimaterialVibEP
