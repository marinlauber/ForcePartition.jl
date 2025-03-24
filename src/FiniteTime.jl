using WaterLily
using StaticArrays
using LinearAlgebra: eigvals

function f(x, t, omega, epsilon)
    return epsilon*sin(omega*t)*x^2 + (1-2*epsilon*sin(omega*t))*x
end
function phi(x, y, t, A=0.1, omega=2pi/10, epsilon=0.2)
    return @. A*sin(pi*f(x, t, omega, epsilon)) * sin(pi*y)
end

function update(r, t, Delta=0.0001)
    x = r[:,1]; y = r[:,2]
    vx = (phi(x,y.+Delta,t).-phi(x,y.-Delta,t))./(2*Delta)
    vy = (phi(x.-Delta,y,t).-phi(x.+Delta,y,t))./(2*Delta)
    return [-vx -vy]
end

function Jacobian(x, y, delta)
    
    # pre-processing
    nx, ny = size(x)
    J = Array{Float64}(undef, 2, 2)
    FTLE = Array{Float64}(undef, nx-2, ny-2)

    for j ∈ 1:ny-2, i ∈ 1:nx-2
        
        # jacobian of that particle
        J[1,1] = (x[i+2,j]-x[i,j])/(2*delta)
        J[1,2] = (x[i,j+2]-x[i,j])/(2*delta)
        J[2,1] = (y[i+2,j]-y[i,j])/(2*delta)
        J[2,2] = (y[i,j+2]-y[i,j])/(2*delta)

        # Green-Cauchy tensor
        D = transpose(J).*J
        
        # its largest eigenvalue
        lamda = eigvals(D)
        FTLE[i,j] = maximum(lamda)
    end
    return FTLE
end

function J(I::CartesianIndex{2},x)
    J = @SMatrix [∂(i,j,I,x) for i ∈ 1:2, j ∈ 1:2]
    maximum(eigvals(J'*J))
end


function test()
    N = 256
    r = zeros(N,N,2)
    apply!(x->x[1],r[:,:,1])
    apply!(x->x[2],r[:,:,2])

    @time for t ∈ T
        @inside r[I] = update!(r,I,t)
    end
end


function main()
    # number of positions
    N = 256

    # particle position (2D grid)
    x = collect(LinRange(0.01, 1.99, 2*N))' .* ones(N)
    y = collect(LinRange(0.01, 0.99,   N))  .* ones(2*N)';

    # particle positions
    r = [vec(x) vec(y)]

    # time step
    dt = 0.1
    T = collect(range(0, 20, step=dt))

    # intergate in time with RK4
    for t ∈ T
        k1 = dt*update(r, t)
        k2 = dt*update(r+0.5*k1, t+0.5*dt)
        k3 = dt*update(r+0.5*k2, t+0.5*dt)
        k4 = dt*update(r+k3, t+dt)
        r += (k1 + 2*k2 + 2*k3 + k4) / 6
    end

    # extract position
    x = reshape(r[:,1], N, 2*N);
    y = reshape(r[:,2], N, 2*N);

    # compute FTLE
    FTLE = Jacobian(x, y, r[1,2]-r[1,1])
    FTLE = log.(FTLE);

    contour(FTLE, levels=31, cmap="viridis")
end
main()