# I use the equations from evolve12.c for my hexagonal-to-orthorhombic model.
# I conserve composition and evolve the three variants using Allen–Cahn equations.
using FFTW, LinearAlgebra, Random, Statistics, DelimitedFiles
Base.@kwdef struct Parameters
    A1::Float64=2.0
    A2::Float64=2.0
    A41::Float64=-3.0
    A42::Float64=36.0
    A61::Float64=2.0
    A62::Float64=4.0
    kappa_c::Float64=0.001
    kappa_eta::NTuple{3,Float64}=(1.0,1.0,1.0)
    mobility::Float64=1.0
    relaxation::NTuple{3,Float64}=(1.0,1.0,1.0)
    eps_c::Float64=0.0
    eps_eta::Float64=0.03
    tetra::Float64=-0.93
    mu::Float64=3000.0
    nu::Float64=0.3
    anisotropy::Float64=1.0
end

function transformation_strains(p)
    e,t = p.eps_eta,p.tetra
    strains = [zeros(3,3) for _ in 1:4]
    strains[1][1,1] = strains[1][2,2] = p.eps_c
    strains[2][1,1] = e
    strains[2][2,2] = e*t
    for (index,sign) in ((3,1),(4,-1))
        strains[index][1,1] = e*(1+3t)/4
        strains[index][2,2] = e*(3+t)/4
        strains[index][1,2] = strains[index][2,1] = sign*e*sqrt(3)*(1-t)/4
    end
    return strains
end

function elastic_tensor(p)
    a,mu,nu = p.anisotropy,p.mu,p.nu
    mu > 0 && a > 0 && -1 < nu < 0.5 || error("Invalid elastic parameters.")
    c11 = mu*(2*(2+a)/(1+a)-(1-4nu)/(1-2nu))
    c12 = mu*(2a/(1+a)-(1-4nu)/(1-2nu))
    c44 = mu*2a/(1+a)
    C = zeros(3,3,3,3)
    for i in 1:3,j in 1:3,k in 1:3,l in 1:3
        if i == j && k == l
            C[i,j,k,l] = i == k ? c11 : c12
        elseif i != j && ((i == k && j == l) || (i == l && j == k))
            C[i,j,k,l] = c44
        end
    end
    return C
end

function make_grid(nx,ny;dx=0.5,dy=0.5)
    kx = collect(2pi .* FFTW.fftfreq(nx,1/dx))
    ky = collect(2pi .* FFTW.fftfreq(ny,1/dy))
    k2 = reshape(kx.^2,nx,1) .+ reshape(ky.^2,1,ny)
    return (;kx,ky,k2,dx,dy)
end

function elastic_kernels(grid,p)
    C = elastic_tensor(p)
    strains = transformation_strains(p)
    stress = [[sum(C[i,j,k,l]*e[k,l] for k in 1:3,l in 1:3)
        for i in 1:3,j in 1:3] for e in strains]
    direct = [sum(strains[a].*stress[b]) for a in 1:4,b in 1:4]
    B = zeros(length(grid.kx),length(grid.ky),4,4)
    for i in eachindex(grid.kx),j in eachindex(grid.ky)
        i == 1 && j == 1 && continue # I allow the mean strain to relax, so I omit the zero-mode energy.
        n = [grid.kx[i],grid.ky[j],0.0]
        n ./= norm(n)
        acoustic = [sum(C[r,s,t,u]*n[s]*n[t] for s in 1:3,t in 1:3)
            for r in 1:3,u in 1:3]
        traction = [s*n for s in stress]
        for a in 1:4,b in 1:4
            B[i,j,a,b] = direct[a,b]-dot(traction[a],acoustic\traction[b])
        end
    end
    # I enforce conjugate symmetry on the Nyquist lines to obtain real fields.
    for a in 1:4,b in 1:4
        field = copy(B[:,:,a,b])
        for i in axes(field,1),j in axes(field,2)
            ii = mod(-(i-1),size(field,1))+1
            jj = mod(-(j-1),size(field,2))+1
            B[i,j,a,b] = (field[i,j]+field[ii,jj])/2
        end
    end
    return B
end

function elastic_derivatives(c,eta,B)
    fields = [fft(c),(fft(e.^2) for e in eta)...]
    potential = [real.(ifft(sum(B[:,:,a,b].*fields[b] for b in 1:4))) for a in 1:4]
    return potential[1],[2 .* eta[a] .* potential[a+1] for a in 1:3]
end

function bulk_derivatives(c,eta,p)
    q = [e.^2 for e in eta]
    dc = p.A1 .* c .- p.A2 .* sum(q)
    de = []
    for a in 1:3
        b,d = filter(index -> index != a,1:3)
        e = eta[a]
        derivative = @. 2p.A2*(1-c)*e + 4p.A41*e^3 +
            2p.A42*e*(q[b]+q[d]) + 6p.A61*e^5 + 2p.A62*e*q[b]*q[d]
        push!(de,derivative)
    end
    return dc,de
end

function total_energy(c,eta,grid,B,p)
    q = [e.^2 for e in eta]
    bulk = @. p.A1*c^2/2 + p.A2*(1-c)*(q[1]+q[2]+q[3]) +
        p.A41*(q[1]^2+q[2]^2+q[3]^2) +
        p.A42*(q[1]*q[2]+q[1]*q[3]+q[2]*q[3]) +
        p.A61*(q[1]^3+q[2]^3+q[3]^3) + p.A62*q[1]*q[2]*q[3]
    hats = [fft(c),(fft(e) for e in eta)...]
    gradient = sum(grid.k2 .* (p.kappa_c .* abs2.(hats[1]) .+
        sum(p.kappa_eta[a-1] .* abs2.(hats[a]) for a in 2:4)))/length(c)
    strain_hats = [hats[1],(fft(e) for e in q)...]
    elastic = 0.0
    for a in 1:4,b in 1:4
        elastic += real(sum(conj.(strain_hats[a]).*B[:,:,a,b].*strain_hats[b]))/(2length(c))
    end
    return (sum(bulk)+gradient+elastic)*grid.dx*grid.dy
end

function step(c,eta,grid,B,p,dt)
    dc,de = bulk_derivatives(c,eta,p)
    ec,ee = elastic_derivatives(c,eta,B)
    new_c = real.(ifft((fft(c).-dt*p.mobility.*grid.k2.*fft(dc.+ec))./
        (1 .+ 2dt*p.mobility*p.kappa_c.*grid.k2.^2)))
    new_eta = [real.(ifft((fft(eta[a]).-dt*p.relaxation[a].*fft(de[a].+ee[a]))./
        (1 .+ 2dt*p.relaxation[a]*p.kappa_eta[a].*grid.k2))) for a in 1:3]
    return new_c,new_eta
end

Base.@kwdef struct Noise
    amplitude_c::Float64=0.1
    amplitude_eta::NTuple{3,Float64}=(0.1,0.1,0.1)
    initial_steps::Int=0
    interval::Int=2500
end

function noise_due(count, noise)
    return count <= noise.initial_steps || (noise.interval > 0 && count % noise.interval == 0)
end

function add_noise!(c, eta, rng, noise)
    # I subtract the composition-noise mean to preserve total composition.
    dc = noise.amplitude_c .* randn(rng,size(c))
    dc .-= mean(dc)
    c .+= dc
    # I add independent nonconserved noise to each variant.
    for a in 1:3
        eta[a] .+= noise.amplitude_eta[a] .* randn(rng,size(c))
    end
end

function initial_fields(nx=32,ny=nx;seed=494,rng=MersenneTwister(seed),
                        composition=0.3,amplitude=0.001)
    # I retain the original uniform composition perturbation and its zero mean.
    dc = amplitude*composition .* (2 .* rand(rng,nx,ny) .- 1)
    c = composition .+ mean(dc) .- dc
    eta = [amplitude .* randn(rng,nx,ny) for _ in 1:3]
    return c,eta
end

function run_simulation(c,eta; p=Parameters(), noise=Noise(), rng=MersenneTwister(494),
        dx=0.5,dy=0.5,steps=100,dt1=0.001,dt2=0.01,time_to_change=1000,
        initial_count=0,initial_time=0.0,directory=joinpath(@__DIR__,"results"))
    steps >= 0 && dt1 > 0 && dt2 > 0 || error("Use nonnegative steps and positive timesteps.")
    size(c) == size(eta[1]) == size(eta[2]) == size(eta[3]) || error("Field sizes differ.")
    c = copy(c); eta = [copy(e) for e in eta]
    grid = make_grid(size(c)...;dx,dy)
    B = elastic_kernels(grid,p)
    # I record energy before and after each noise event separately.
    history = zeros(steps+1,5)
    time = initial_time
    history[1,:] = [time,mean(c),total_energy(c,eta,grid,B,p),0,NaN]
    for s in 1:steps
        count = initial_count+s-1
        dt = count > time_to_change ? dt2 : dt1
        c,eta = step(c,eta,grid,B,p,dt)
        energy_before_noise = total_energy(c,eta,grid,B,p)
        event = noise_due(count,noise)
        event && add_noise!(c,eta,rng,noise)
        time += dt
        all(isfinite,c) && all(e -> all(isfinite,e),eta) || error("Non-finite fields at step $s.")
        history[s+1,:] = [time,mean(c),total_energy(c,eta,grid,B,p),Int(event),energy_before_noise]
    end
    if directory !== nothing
        mkpath(directory)
        writedlm(joinpath(directory,"history.csv"),history,',')
        writedlm(joinpath(directory,"composition.csv"),c,',')
        for a in 1:3
            writedlm(joinpath(directory,"variant-$a.csv"),eta[a],',')
        end
    end
    return c,eta,history
end

function run_example(;n=32,steps=100,dt=0.001,noise=Noise(),
        directory=joinpath(@__DIR__,"results"))
    # I keep one seeded generator for initialization and subsequent noise.
    rng = MersenneTwister(494)
    c,eta = initial_fields(n;rng)
    return run_simulation(c,eta;noise,rng,steps,dt1=dt,dt2=dt,directory)
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_example()
end

# I read the original labelled parameter file, including repeated noise labels.
function run_from_input(path;directory=joinpath(@__DIR__,"results-input"))
    values = Dict{String,Float64}()
    variant_noise = Float64[]
    for line in eachline(path)
        words = split(line)
        length(words) < 2 && continue
        number = tryparse(Float64,words[2])
        number === nothing && continue
        values[words[1]] = number
        startswith(words[1],"sustained_noise_level_uncons") && push!(variant_noise,number)
    end
    length(variant_noise) == 3 || error("Specify all three variant noise amplitudes.")
    get(values,"flag",0) == 0 || error("For this initialization branch, load the original fields and call run_simulation directly.")
    get(values,"initcount",0) == 0 || error("For a restart, load the saved fields and supply initial_count and initial_time.")
    p = Parameters(A1=values["A1"],A2=values["A2"],A41=values["A41"],A42=values["A42"],
        A61=values["A61"],A62=values["A62"],kappa_c=values["kappa"],
        kappa_eta=Tuple(values["kappa$a"] for a in 1:3),mobility=values["mobility"],
        relaxation=Tuple(values["L$a"] for a in 1:3),eps_c=values["epsc"],
        eps_eta=values["eps_eta"],tetra=values["tetra"],mu=values["mubar"],
        nu=values["nubar"],anisotropy=values["Aniso"])
    noise = Noise(amplitude_c=values["sustained_noise_level_cons"],
        amplitude_eta=Tuple(variant_noise),initial_steps=Int(values["noise_steps"]))
    rng = MersenneTwister(abs(Int(values["SEED"])))
    c,eta = initial_fields(Int(values["nx"]),Int(values["ny"]);rng,
        composition=values["alloycomp"],amplitude=values["noise_level"])
    return run_simulation(c,eta;p,noise,rng,dx=values["del_x"],dy=values["del_y"],
        steps=Int(values["num_steps"]),dt1=values["del_t1"],dt2=values["del_t2"],
        time_to_change=Int(values["time_to_change"]),directory)
end
