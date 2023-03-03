function binary!(;networks)
    for n in 1:length(networks)

       @views net=networks[n]

       for j in 1:size(net,2)
            for i in 1:size(net,1)
                if net[i,j]>0.0
                    net[i,j]=1.0
                end
            end
        end
    end
end

function adjacency(;incidence_matrices::AbstractArray)

    adj=deepcopy(incidence_matrices)
    
    for n in 1:length(incidence_matrices)

        incidence_matrix=incidence_matrices[n]

        n_p=size(incidence_matrix, 1)
        n_a=size(incidence_matrix, 2)

        A=vcat(hcat(zeros(n_p, n_p), incidence_matrix), hcat(incidence_matrix', zeros(n_a, n_a)))

        adj[n]=A

    end

    return adj

end

function native(incidence_matrices, invader_position)

    nat=deepcopy(incidence_matrices)

    for n in 1:length(incidence_matrices)

        ap=invader_position[n]
        @views @. nat[n][:, ap]=0.0

    end

    return nat

end

function networks_description(;networks, networks_native, net_name, matrix_type, invader_position)

    results=[DataFrame() for _ in 1:length(networks)]

    for i in 1:length(networks)

        network=networks[i]
        net_native=networks_native[i]

        adj=adjacency(incidence_matrices=[network, net_native])

        if matrix_type=="incidence"

            n_p=size(network, 1)
            n_a=size(network, 2)
            n_sp=n_p+n_a
                    
        elseif matrix_type=="adjacency"
        
            n_sp=size(network,1)
            n_p=0
            n_a=0            

        end      
                
        if matrix_type=="incidence"

            p_degree=dropdims(sum(network, dims=2), dims=2)
            a_degree=dropdims(sum(network, dims=1), dims=1)

            p_native_degree=dropdims(sum(net_native, dims=2), dims=2)
            a_native_degree=dropdims(sum(net_native, dims=1), dims=1)

            degree=vcat(p_degree, a_degree)
            degree_native=vcat(p_native_degree, a_native_degree)

            cp_plant=ifelse.(p_degree.>=mean(p_degree), "core", "periphery")
            cp_animal=ifelse.(a_degree.>=mean(a_degree), "core", "periphery")

            cp=vcat(cp_plant, cp_animal)
        
        elseif matrix_type=="adjacency"
                
            cp=["unknown" for _ in 1:n_sp]
            degree=dropdims(sum(network, dims=2), dims=2)

        end

        results[i]=DataFrame(A=repeat([adj[1]], n_sp), A_n=repeat([adj[2]], n_sp), network=net_name[i], n_p=n_p, n_a=n_a, n_sp=n_sp, degree=degree, degree_native=degree_native, cp=cp, invader=invader_position[i])

        
    end

    return results
        
end

    

function pairwise_distances!(;dz::Array{Float64}, A::Array{Float64}, D::Array{Float64}, time::Int64)

    fill!(D, 0.0)
    
    for j in 1:size(D,2)
        for i in 1:size(D,1)
            D[i,j]=A[i,j]*(dz[time,j]-dz[time,i])
        end
    end
end

function Qij!(;D::Array{Float64}, A::Array{Float64}, Q::Array{Float64}, Q_sum::Array{Float64}, α::Float64, m::Float64)
    
    fill!(Q, 0.0)
    fill!(Q_sum, 0.0)

    @inbounds for j in 1:size(Q,2)
        @inbounds for i in 1:size(Q,1)

            Q[i,j]=A[i,j]*exp(-α*(D[i,j]^2))

        end
    end

    sum!(Q_sum, Q)

    @inbounds for j in 1:size(Q,2)
        @inbounds for i in 1:size(Q,1)

            Q[i,j]=(m*(Q[i,j]/Q_sum[i]))*D[i,j]

        end
    end

    sum!(Q_sum, Q)

end

function fitness!(;w::Array{Float64}, dz::Array{Float64}, θ::Array{Float64}, D::Array{Float64}, A::Array{Float64}, degree::Array{Float64}, Q::Array{Float64}, Q_sum::Array{Float64}, α::Float64, m::Float64, ρ::Float64, time::Int64)
    
    fill!(D, 0.0)
    fill!(Q, 0.0)
    fill!(Q_sum, 0.0)

    pairwise_distances!(dz=dz, D=D, A=A, time=time)
    @. Q=A*(exp(-α*(D^2)))
    sum!(Q_sum, Q)

    for i in 1:length(Q_sum)

        w[i]=exp((ρ/2.0)*(((m/α)*log(Q_sum[i]/degree[i])) - (1.0-m)*(θ[i] - dz[time,i])^2))

    end   

end

function evolve!(;dz::Array{Float64}, zdif::Array{Float64}, θ::Array{Float64}, Q_sum::Array{Float64}, m::Float64, ϕ::Float64, time::Int64)

    for i in 1:size(dz,2)

        dz[time+1, i] = dz[time, i]  + (ϕ*(Q_sum[i] + ((1.0-m)*(θ[i]-dz[time, i]))))

        zdif[i]=abs(dz[time+1, i] - dz[time, i])

    end
end

function T_matrix(;dz::Array{Float64}, A::Array{Float64}, D::Array{Float64}, Q::Array{Float64}, Q_sum::Array{Float64}, m::Float64, α::Float64, time::Int64)

    fill!(Q, 0.0)
    fill!(Q_sum, 0.0)

    pairwise_distances!(dz=dz, A=A, D=D, time=time)

    @. Q=A*(exp(-α*(D^2)))
    sum!(Q_sum, Q)
    @. Q=(m*(Q/Q_sum))

    Ψ=Diagonal(repeat([(1.0-m)], length(Q_sum)))

    T=Ψ/(I-Q)

    return T

end

function T_matrix_native(;dz::Array{Float64}, A::Array{Float64}, inv_pos::Int64, D::Array{Float64}, Q::Array{Float64}, Q_sum::Array{Float64}, m::Float64, α::Float64, time::Int64)

    fill!(Q, 0.0)
    fill!(Q_sum, 0.0)

    dz[time,inv_pos]=0.0

    pairwise_distances!(dz=dz, A=A, D=D, time=time)

    @. Q=A*(exp(-α*(D^2)))
    sum!(Q_sum, Q)
    Q_sum[inv_pos]=1.0
    @. Q=(m*(Q/Q_sum))

    Ψ=Diagonal(repeat([(1.0-m)], length(Q_sum)))

    T=Ψ/(I-Q)

    return T

end

function Qij2!(;D::Array{Float64}, A::Array{Float64}, Q::Array{Float64}, Q_sum::Array{Float64}, α::Float64, m::Float64)
    
    fill!(Q, 0.0)
    fill!(Q_sum, 0.0)

    @. Q=A*(exp(-α*(D^2)))
    sum!(Q_sum, Q)
    @. Q=(m*(Q/Q_sum))*D
    sum!(Q_sum, Q)

end