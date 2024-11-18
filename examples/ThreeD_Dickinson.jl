using WaterLily
using StaticArrays
using LinearAlgebra
using CUDA
using ForcePartition
using WriteVTK
using Plots

# make simulation foloowing Dickinson's setup
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
    body =  elipsoid ∩ upper_lower # intersection of sets

    # Return initialized simulation
    Simulation((4L,3L,4L),(0,0,0),L;ν,U,body,mem,T)
end

# make a simulation
sim = Dickinson(4;mem=CUDA.CuArray,T=Float32)

# make a writer with some attributes, need to output to CPU array to save file (|> Array)
velocity(a::Simulation) = a.flow.u |> Array;
pressure(a::Simulation) = a.flow.p |> Array;
_body(a::Simulation) = (measure_sdf!(a.flow.σ, a.body, WaterLily.time(a));
                                     a.flow.σ |> Array;)
lamda(a::Simulation) = (@inside a.flow.σ[I] = WaterLily.λ₂(I, a.flow.u);
                        a.flow.σ |> Array;)

custom_attrib = Dict(
    "Velocity" => velocity,
    "Pressure" => pressure,
    "Body" => _body,
    "Lambda" => lamda
)# this maps what to write to the name in the file
# make the writer
# writer = vtkWriter("Dickinson"; attrib=custom_attrib)

# force moment
fpm = ForcePartitionMethod(sim)
potential!(fpm;axis=2) # change axis to 2 for lift

# a simulation
t₀ = sim_time(sim); duration = 12; tstep = 0.05 # print time
forces,fp,fm,fv = [],[],[],[] # empty list to store forces

# simulate
@time for tᵢ in range(t₀,t₀+duration;step=tstep)
    @show tᵢ
    # update until time tᵢ in the background
    t = sum(sim.flow.Δt[1:end-1])
    # this is normaly hiden into sim_step!(sim,tᵢ)
    while t < tᵢ*sim.L/sim.U
        #update the body
        measure!(sim,t)
        # update flow
        mom_step!(sim.flow,sim.pois) 
        # update time
        t += sim.flow.Δt[end]
        # compute and save pressure forces, this is actually a force coefficient, 
        # to find the actual force we have Fi = Ci/(1/2*ρ*U^2*L^2)
        force = -2WaterLily.pressure_force(sim)/sim.L^2
        push!(forces,[t/sim.L,force...])
        push!(fp,-∫2QϕdV!(fpm,sim.flow,recompute=true,axis=4))
        push!(fm,-∮UϕdS!(fpm,sim.flow,recompute=false,axis=4))
        push!(fv, ∮ReωdS!(fpm,sim.flow,recompute=false,axis=4))
    end
    # this writes every tstep to paraview files, not every time step
    # write!(writer,sim);
end
forces = reduce(vcat,forces')
println("Done...")
# close(writer)

# plot results
plot(cumsum(@views(sim.flow.Δt[1:end-1]))./sim.L,forces[:,2],label="Total force",
     xlabel="tU/L",ylabel="2F/ρU²L²") #,ylims=(1,3),xlims=(0,12))
plot!(cumsum(@views(sim.flow.Δt[1:end-1]))./sim.L,fp/sim.L^2,label="vorticity force")
plot!(cumsum(@views(sim.flow.Δt[1:end-1]))./sim.L,fm/sim.L^2,label="added-mass force")
plot!(cumsum(@views(sim.flow.Δt[1:end-1]))./sim.L,fv/sim.L^2,label="viscous force")
plot!(cumsum(@views(sim.flow.Δt[1:end-1]))./sim.L,(fp+fm+fv)/sim.L^2,label="total partition")
savefig("force_partition_mosquito.png")