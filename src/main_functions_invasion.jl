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
    #dw=zeros(tmax, n_sp)
    D=zeros(n_sp, n_sp)
    Q=zeros(n_sp, n_sp)
    Q_sum=zeros(n_sp)
    w0=zeros(n_sp)
    #weq_1=zeros(n_sp)
    weq_2=zeros(n_sp)
    weq_3=zeros(n_sp)

    teq=0
    teq_inv=0
    inv=0.0

    zdif=zeros(n_sp)
    
    z0=rand(0.0:0.001:10.0, n_sp)
    θ=rand(0.0:0.001:10.0, n_sp)
    #θ[invader]=5.0

    @views dz[1,:].=z0
    fitness!(w=w0, dz=dz, θ=θ, D=D, A=A_n, degree=degree_native, Q=Q, Q_sum=Q_sum, α=α, m=m, ρ=ρ, time=1)

    #@views dw[1,:].=fitness!(z=z0, θ=θ, w=w, D=D, A=A, degree=degree, Q=Q, Q_sum=Q_sum, α=α, m=m, ρ=ρ)

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
            #fitness!(w=weq_1, dz=dz, θ=θ, D=D, A=A_n, degree=degree_native, Q=Q, Q_sum=Q_sum, α=α, m=m, ρ=ρ, time=teq_inv)
            fitness!(w=weq_2, dz=dz, θ=θ, D=D, A=A, degree=degree, Q=Q, Q_sum=Q_sum, α=α, m=m, ρ=ρ, time=(teq_inv+1))
            inv=inv+1.0

        elseif (NaNMath.mean(zdif) < ϵ) && (inv>0.0)
            teq=t
            break            
        end

    end

    fitness!(w=weq_3, dz=dz, θ=θ, D=D, A=A, degree=degree, Q=Q, Q_sum=Q_sum, α=α, m=m, ρ=ρ, time=teq)

    @views zeq=dz[teq,:]

    teq_afterinv=teq-teq_inv
    #wdif=weq_3.-w0
    #w_loss_freq=ifelse.(wdif.<0.0, 1.0, 0.0)
    #w_loss_amount=ifelse.(wdif.<0.0, wdif, 0.0)

    #wdif_a1=weq_2.-weq_1
    #w_loss_freq_a1=ifelse.(wdif_a1.<0.0, 1.0, 0.0)
    #w_loss_amount_a1=ifelse.(wdif_a1.<0.0, wdif_a1, 0.0)

    wdif=weq_3.-weq_2
    wl_freq=ifelse.(wdif.<0.0, 1.0, 0.0)
    wl_mag=ifelse.(wdif.<0.0, wdif, 0.0)
    wl_mag.=abs.(wl_mag)
    wg_mag=ifelse.(wdif.>0.0, wdif, 0.0)
    wg_mag.=abs.(wg_mag)

    #wdif_a3=weq_3.-weq_1
    #w_loss_freq_a3=ifelse.(wdif_a3.<0.0, 1.0, 0.0)
    #w_loss_amount_a3=ifelse.(wdif_a3.<0.0, wdif_a3, 0.0)

    T=T_matrix(dz=dz, A=A, D=D, Q=Q, Q_sum=Q_sum, m=m, α=α, time=teq)
    Tn=T_matrix_native(dz=dz, A=A_n, D=D, inv_pos=invader, Q=Q, Q_sum=Q_sum, m=m, α=α, time=teq_inv)

    T_diag=deepcopy(T)
    T_diag[diagind(T_diag)].=0.0
    Tin_diag=dropdims(sum(T_diag, dims=2), dims=2)

    T_ind=T_diag.*(1.0.-A)
    Tin_ind=dropdims(sum(T_ind, dims=2), dims=2)
    ind_prop=Tin_ind./Tin_diag

    Tn_diag=deepcopy(Tn)
    Tn_diag[diagind(Tn_diag)].=0.0
    Tnin_diag=dropdims(sum(Tn_diag, dims=2), dims=2)

    Tn_ind=Tn_diag.*(1.0.-A_n)
    Tnin_ind=dropdims(sum(Tn_ind, dims=2), dims=2)
    ind_prop_native=Tnin_ind./Tnin_diag

    Tout=dropdims(mean(T, dims=1), dims=1)
    Tin_apis=T[:,invader]
    Tout_apis=Tout[invader]

    if n_p>0 && n_a>0
        type=vcat(repeat(["plant"], n_p), repeat(["animal"], n_a))
    elseif n_p==0 && n_a==0
        type=repeat(["unknown"], n_sp)
    end

    sp_id=["SP$i" for i in 1:n_sp]
    sp_id[invader]="SP_Apis"

    #df=DataFrame(network=net_name, scenario=scenario, sp_id=sp_id, type=type, cp=cp, degree=degree, degree_native=degree_native, theta=θ, z0=z0, zeq=zeq, w0=w0, weq_1=weq_1, weq_2=weq_2, weq_3=weq_3, wdif=wdif, 
    #w_loss_freq=w_loss_freq, w_loss_amount=w_loss_amount, wdif_a1=wdif_a1, w_loss_freq_a1=w_loss_freq_a1, w_loss_amount_a1=w_loss_amount_a1, wdif_a2=wdif_a2, w_loss_freq_a2=w_loss_freq_a2,
    #w_loss_amount_a2=w_loss_amount_a2, wdif_a3=wdif_a3, w_loss_freq_a3=w_loss_freq_a3, w_loss_amount_a3=w_loss_amount_a3,
    #evo_broadcast=Tout, Tin_apis=Tin_apis, Tout_apis=Tout_apis, T_theta=T_θ, teq=teq, teq_inv=teq_inv, teq_afterinv=teq_afterinv, m=m, alpha=α, phi=ϕ, rho=ρ, sim=sim)

    df=DataFrame(network=net_name, scenario=scenario, sp_id=sp_id, type=type, cp=cp, degree=degree, degree_native=degree_native, theta=θ, z0=z0, zeq=zeq, w0=w0, weq_inv=weq_2, weq_final=weq_3, wdif=wdif, 
    wl_freq=wl_freq, wl_mag=wl_mag, wg_mag=wg_mag, indp=ind_prop, indp_native=ind_prop_native, Tin_apis=Tin_apis, Tout_apis=Tout_apis, teq=teq, teq_inv=teq_inv, teq_afterinv=teq_afterinv, m=m, alpha=α, phi=ϕ, rho=ρ, sim=sim)

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

