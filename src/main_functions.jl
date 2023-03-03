function net_coevolution(;A::Array{Float64}, net_info::DataFrame, m::Float64, α::Float64, ϕ::Float64, ρ::Float64, θ_max::Float64, tmax::Int, ϵ::Float64, sim::Int, scenario::String)

    net_name=net_info.network[1] #Name of the network
    n_a=net_info.n_a[1] #Number of animals
    n_p=net_info.n_p[1] #Number of plants
    n_sp=net_info.n_sp[1] #Total number of species
    degree=net_info.degree #Number of direct partners (degree) for each species
    cp=net_info.cp #Species position (core/periphery) in the network

    dz=zeros(tmax, n_sp) #Creating array to store trait values over time
    D=zeros(n_sp, n_sp) #Initializing matrices to help compute the coevolutionary dynamics
    Q=zeros(n_sp, n_sp) #Initializing matrices to help compute the coevolutionary dynamics
    Q_sum=zeros(n_sp) #Initializing matrices to help compute the coevolutionary dynamics
    w0=zeros(n_sp) #Initializing vector to store species initial fitness
    weq=zeros(n_sp) #Initializing vector to store species fitness at equilibrium
    wmax=zeros(n_sp) #Initializing vector to store the theoretical maximum fitness for each species

    teq=0 #Initializing variable to store the time until equilibrium
    zdif=zeros(n_sp) #Initializing variable that will be used to stop the coevolutionary dynamics (difference in trait values between successive generations)
    
    z0=rand(0.0:0.001:10.0, n_sp) #Sampling initial trait values
    θ=rand(0.0:0.001:θ_max, n_sp) #Sampling theta from a uniform distribution within the range U[0.0, θ_max]

    @views dz[1,:].=z0
    fitness!(w=w0, dz=dz, θ=θ, D=D, A=A, degree=degree, Q=Q, Q_sum=Q_sum, α=α, m=m, ρ=ρ, time=1)
    max_fitness!(wmax=wmax, θ=θ, degree=degree, α=α, m=m, ρ=ρ)

    #Coevolutionary dynamics
   @inbounds for t in 1:(tmax-1)

        pairwise_distances!(dz=dz, A=A, D=D, time=t)
        Qij!(D=D, A=A, Q=Q, Q_sum=Q_sum, α=α, m=m)
        evolve!(dz=dz, zdif=zdif, θ=θ, Q_sum=Q_sum, m=m, ϕ=ϕ, time=t)        

        if mean(zdif) < ϵ
            teq=t
            break
        end

        teq=tmax

    end

    fitness!(w=weq, dz=dz, θ=θ, D=D, A=A, degree=degree, Q=Q, Q_sum=Q_sum, α=α, m=m, ρ=ρ, time=teq) #Computing species fitness at equilibrium

    @views zeq=dz[teq,:] #Getting equilibrium trait values

    T=T_matrix(dz=dz, A=A, D=D, Q=Q, Q_sum=Q_sum, m=m, α=α, time=teq) #Computing T-matrix
    T[diagind(T)].=0.0 #Setting diagonal equal to zero to compute the total amount of incoming evolutionary effects
    Tin=dropdims(sum(T, dims=2), dims=2) #Computing the total amount of incoming evolutionary effects

    T_ind=T_diag.*(1.0.-A) # Computing the matrix of indirect evolutionary effects
    Tin_ind=dropdims(sum(T_ind, dims=2), dims=2) #Computing the incoming amounf of indirect evolutionary effects
    indirect_effects=Tin_ind./Tin #Computing the contribution of indirect evolutionary effects relative to the total amounf of incoming evolutionary effects
    
    type=vcat(repeat(["plant"], n_p), repeat(["animal"], n_a)) #Setting species type (animal or plants)
   

    sp_id=["SP$i" for i in 1:n_sp] #Setting an ID for species
    
    df=DataFrame(network=net_name, scenario=scenario, sp_id=sp_id, type=type, position=cp, 
    degree=degree, theta=θ, z0=z0, zeq=zeq, wmax=wmax, w0=w0, weq=weq, indirect_effects=indirect_effects, 
    teq=teq, m=m, alpha=α, phi=ϕ, rho=ρ, theta_max=θ_max, sim=sim) #Creating data frame to store results

    return df #Returning data frame with results

end