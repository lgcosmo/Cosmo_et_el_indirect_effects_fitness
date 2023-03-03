using Distributed

addprocs(20)

@everywhere using DrWatson
@everywhere @quickactivate "coevo_fit"

@everywhere begin

    using DataFrames
    using Distributions
    using Statistics
    using LinearAlgebra
    using DelimitedFiles
    using CSV
    using Random
    using NaNMath

    BLAS.set_num_threads(1)

    include(srcdir("main_functions_invasion.jl"))
    include(srcdir("aux_functions_invasion.jl"))

    mi_settings=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

    networks_paths_original=readdir(srcdir("empirical_networks", "apis_invasion"); join=true)

    networks_names=readdir(srcdir("empirical_networks", "apis_invasion"))
    networks_names=replace.(networks_names, r".csv.*" => "")
    networks_names=networks_names[1:73]

    files_invasion=readdlm.(networks_paths_original, ',')
    networks_invasion=files_invasion[1:73]
    apis_adj=files_invasion[74]
    apis_adj=Int64.(apis_adj)
    apis_inc=files_invasion[75]
    apis_inc=Int64.(apis_inc)

    binary!(networks=networks_invasion)
    networks_native=native(networks_invasion, apis_inc)

    invasion_descriptors=networks_description(networks=networks_invasion, networks_native=networks_native, net_name=networks_names, matrix_type="incidence", invader_position=apis_adj)

    c=Dict(:net_info => invasion_descriptors, :m => mi_settings, :scenario => "invasion", :α => 0.2, :ϕ => 0.2, :ρ => 0.2, :tmax => 10000, :ϵ => 1e-5, :nsim => 1000)
    c_list=dict_list(c)

end

pmap(net_multisim, c_list)


