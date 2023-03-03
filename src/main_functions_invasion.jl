function net_coevolution(;A::Array{Float64}, A_n::Array{Float64}, net_info::DataFrame, m::Float64, α::Float64, ϕ::Float64, ρ::Float64, tmax::Int, ϵ::Float64, sim::Int, scenario::String)

    net_name=net_info.network[1]
    n_a=net_info.n_a[1]
    n_p=net_info.n_p[1]
    n_sp=net_info.n_sp[1]
    cp=net_info.cp
    degree=net_info.degree
    degree_native=net_info.degree_native
    invader=net_info.invader[1]

    dz=zeros(tmax, n_sp)
    D=zeros(n_sp, n_sp)
    Q=zeros(n_sp, n_sp)
    Q_sum=zeros(n_sp)
    w0=zeros(n_sp)
    weq_preinvasion=zeros(n_sp)
    weq_postinvasion=zeros(n_sp)

    teq=0
    teq_inv=0
    inv=0.0

    zdif=zeros(n_sp)
    
    z0=rand(0.0:0.001:10.0, n_sp)
    θ=rand(0.0:0.001:10.0, n_sp)

    @views dz[1,:].=z0
    fitness!(w=w0, dz=dz, θ=θ, D=D, A=A_n, degree=degree_native, Q=Q, Q_sum=Q_sum, α=α, m=m, ρ=ρ, time=1)

   @inbounds for t in 1:(tmax-1)

    if inv==0.0
        dz[t,invader]=0.0
        pairwise_distances!(dz=dz, A=A_n, D=D, time=t)
        Qij!(D=D, A=A_n, Q=Q, Q_sum=Q_sum, α=α, m=m)
        evolve!(dz=dz, zdif=zdif, θ=θ, Q_sum=Q_sum, m=m, ϕ=ϕ, time=t)
    
    elseif inv==1.0
        pairwise_distances!(dz=dz, A=A, D=D, time=t)
        Qij!(D=D, A=A, Q=Q, Q_sum=Q_sum, α=α, m=m)
        evolve!(dz=dz, zdif=zdif, θ=θ, Q_sum=Q_sum, m=m, ϕ=ϕ, time=t)
    end

        if (NaNMath.mean(zdif) < ϵ) && (inv==0.0)
            
            teq_inv=t
            dz[t+1,invader]=dz[1,invader]
            fitness!(w=weq_preinvasion, dz=dz, θ=θ, D=D, A=A, degree=degree, Q=Q, Q_sum=Q_sum, α=α, m=m, ρ=ρ, time=(teq_inv+1))
            inv=inv+1.0

        elseif (NaNMath.mean(zdif) < ϵ) && (inv>0.0)
            teq=t
            break            
        end

    end

    fitness!(w=weq_postinvasion, dz=dz, θ=θ, D=D, A=A, degree=degree, Q=Q, Q_sum=Q_sum, α=α, m=m, ρ=ρ, time=teq)

    @views zeq=dz[teq,:]

    teq_afterinv=teq-teq_inv

    wdif=weq_postinvasion.-weq_preinvasion

    T_post=T_matrix(dz=dz, A=A, D=D, Q=Q, Q_sum=Q_sum, m=m, α=α, time=teq)
    T_pre=T_matrix_native(dz=dz, A=A_n, D=D, inv_pos=invader, Q=Q, Q_sum=Q_sum, m=m, α=α, time=teq_inv)

    T_post[diagind(T_post)].=0.0
    T_pre[diagind(T_post)].=0.0

    Tin_post=dropdims(sum(T_post, dims=2), dims=2)
    Tin_pre=dropdims(sum(T_pre, dims=2), dims=2)

    T_ind_post=T_post.*(1.0.-A)
    T_ind_pre=T_pre.*(1.0.-A_n)

    indirect_effects_post=dropdims(sum(T_ind_post, dims=2), dims=2)./Tin_post
    indirect_effects_pre=dropdims(sum(T_ind_pre, dims=2), dims=2)./Tin_pre

    type=vcat(repeat(["plant"], n_p), repeat(["animal"], n_a))
   
    sp_id=["SP$i" for i in 1:n_sp]
    sp_id[invader]="SP_Apis"

    df=DataFrame(network=net_name, scenario=scenario, sp_id=sp_id, type=type, cp=cp, degree=degree, degree_native=degree_native, 
    theta=θ, z0=z0, zeq=zeq, w0=w0, weq_preinvasion=weq_preinvasion, weq_postinvasion=weq_postinvasion, wdif=wdif, 
    indirect_effects_post=indirect_effects_post, indirect_effects_pre=indirect_effects_pre, teq=teq, teq_inv=teq_inv, 
    teq_afterinv=teq_afterinv, m=m, alpha=α, phi=ϕ, rho=ρ, sim=sim)

    return df

end



function net_multisim(p)
    
    @unpack α,m,ρ,tmax,ϵ,ϕ,nsim,net_info,scenario = p
    
    A=net_info.A[1]
    A_n=net_info.A_n[1]
    net_name=net_info.network[1]

    r=[DataFrame() for _ in 1:nsim]

    for n in 1:nsim
        r[n]=net_coevolution(A=A, A_n=A_n, net_info=net_info, m=m, α=α, ϕ=ϕ, ρ=ρ, tmax=tmax, ϵ=ϵ, sim=n, scenario=scenario)
    end

    CSV.write(datadir("sims", "$(scenario)_$(net_name)_m$(m).csv"), vcat(r...))

end

function coevo_simulations(;A_list, N_list, p_list)

    Threads.@threads for i in 1:length(p_list)
        
        net_multisim(A_list=A_list, N_list=N_list, p=p_list[i])
        
    end
end

