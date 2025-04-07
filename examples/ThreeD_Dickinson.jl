using WaterLily,StaticArrays,LinearAlgebra
using CUDA,ForcePartition,WriteVTK,Plots

# strating from rest
WaterLily.CFL(a::Flow) = WaterLily.CFL(a;Δt_max=0.5)

# make simulation following Dickinson's setup
function Dickinson(p=5;Re=5e2,mem=Array,T=Float32)
    # Define simulation size, geometry dimensions, & viscosity
    L = 2^p; U=1; ν = U*L/Re
    ϵ = 0.5; thk = 2ϵ+√3
    AR = 2.0
    function sdf(x,t)
        # offset it to put origin at tip, and stretch it by AR. remove radius for sdf
        √sum(abs2, SA[x[1], (x[2]-1.5L/2)/AR, 0]) - L/2/AR
    end
    function wing(xyz,t)
        # put in symmetry, correct scaling
        x = abs.(SA[xyz[1],xyz[2]-1.5L/2]); ab=SA[L/2AR,L/2]

        # find root with Newton solver
        q = ab.*(x-ab);
        w = (q[1]<q[2]) ? π/2 : 0.0;
        for i ∈ 1:5
            u = ab.*SA[ cos(w),sin(w)]
            v = ab.*SA[-sin(w),cos(w)]
            w += dot(x-u,v)/(dot(x-u,u)+dot(v,v));
        end
        # compute final point and distance
        d = norm(x-ab.*SA[cos(w),sin(w)]);
        
        # return signed distance
        return (dot(x./ab,x./ab)>1.0) ? d : -d
    end
    
    # the mapping
    function map(x,t)
        # Dickinson kinemtics
        _α = π/2 - π/4*sin(π*t/L)  # positive pitch increase AoA
        _ϕ = 0.35π*cos(π*t/L)      # positive to the rear of the mosquito 
        # rotation mmatrix
        Ry = SA[cos(_α) 0 sin(_α); 0 1 0; -sin(_α) 0 cos(_α)] # alpha
        Rz = SA[cos(_ϕ) -sin(_ϕ) 0; sin(_ϕ) cos(_ϕ) 0; 0 0 1] # phi
        return Ry*Rz*(x .- SA[2L,0,L]) # the order matters
    end
   
    # Build the mosquito from a mapped elipsoid and two plane that trim it to the correct thickness
    elipsoid = AutoBody(sdf, map)      
    upper_lower = AutoBody((x,t)->(abs(x[3])-thk/2), map)
    # this creates the final body
    body1 = elipsoid ∩ upper_lower # intersection of sets
    body2 = AutoBody(sdf, (x,t)->map(x.+SA[L,0,0],t+π/2)) ∩ AutoBody((x,t)->(abs(x[3])-thk/2), (x,t)->map(x.+SA[L,0,0],t+π/2))
    body = body1 ∪ body2 # union of sets

    # Return initialized simulation
    Simulation((4L,3L,4L),(0,0,0),L;ν,U,body,mem,T),body1,body2
end

# make a simulation
sim,body1,body2 = Dickinson(4;mem=CUDA.CuArray,T=Float32)

# make a writer with some attributes, need to output to CPU array to save file (|> Array)
vtk_velocity(a::AbstractSimulation) = a.flow.u |> Array;
vtk_pressure(a::AbstractSimulation) = a.flow.p |> Array;
vtk_body(a::AbstractSimulation) = (measure_sdf!(a.flow.σ, a.body, WaterLily.time(a)); a.flow.σ |> Array;)
vtk_lambda(a::AbstractSimulation) = (@inside a.flow.σ[I] = WaterLily.λ₂(I, a.flow.u); a.flow.σ |> Array;)
vtk_Q(a::AbstractSimulation) = (@inside a.flow.σ[I] = ForcePartition.Qcriterion(I, a.flow.u); a.flow.σ |> Array)
vtk_ϕ1(a::AbstractSimulation) = (potential!(fpm,body1;axis=2,tᵢ=sum(a.flow.Δt)); fpm.ϕ |> Array)
vtk_ϕ2(a::AbstractSimulation) = (potential!(fpm,body2;axis=2,tᵢ=sum(a.flow.Δt)); fpm.ϕ |> Array)

custom_attrib = Dict("u" => vtk_velocity,"p" => vtk_pressure,"d" => vtk_body,
                     "λ2" => vtk_lambda,"Q" => vtk_Q,"ϕ1" => vtk_ϕ1,"ϕ2"=> vtk_ϕ2)# this maps what to write to the name in the file
# make the writer
writer = vtkWriter("Dickinson"; attrib=custom_attrib)

# force moment
fpm = ForcePartitionMethod(sim)
potential!(fpm,fpm.body;axis=2) # change axis to 2 for lift

# a simulation
t₀ = sim_time(sim); duration = 12; tstep = 0.05 # print time
forces,fp,fm,fv = [],[],[],[] # empty list to store forces

# simulate
@time for tᵢ in range(t₀,t₀+duration;step=tstep)
    # update until time tᵢ in the background
    t = sum(sim.flow.Δt[1:end-1])
    # this is normaly hiden into sim_step!(sim,tᵢ)
    while t < tᵢ*sim.L/sim.U
        #update the body
        sim_step!(sim;remasure=true)
        # update time
        t += sim.flow.Δt[end]
        # compute and save pressure forces, this is actually a force coefficient, 
        # to find the actual force we have Fi = Ci/(1/2*ρ*U^2*L^2)
        force = -2WaterLily.pressure_force(sim)
        push!(forces,[t/sim.L,force...])
        # the totential is the same, we just recompute it once and then use it
        push!(fp,-∫2QϕdV!(fpm,sim.flow,recompute=true,axis=2))
        push!(fm,-∮UϕdS!(fpm,sim.flow,recompute=false,axis=2))
        push!(fv, ∮ReωdS!(fpm,sim.flow,recompute=false,axis=2))
    end
    # this writes every tstep to paraview files, not every time step
    write!(writer,sim);
    println("tU/L=",round(tᵢ,digits=4),", Δt=",round(sim.flow.Δt[end],digits=3))
end
forces = reduce(vcat,forces')
println("Done...")
close(writer)

# plot results
time = cumsum(sim.flow.Δt)./sim.L # time vector
plot(time[1:end-1],[forces[:,2]./sim.L^2, fp/sim.L^2, fm/sim.L^2, fv/sim.L^2, (fp+fm+fv)/sim.L^2],
     label=["Total force" "vorticity force" "added-mass force" "viscous force" "total partition"],
     xlabel="tU/L",ylabel="2F/ρU²L²") #,ylims=(1,3),xlims=(0,12))
savefig("force_partition_mosquito.png")