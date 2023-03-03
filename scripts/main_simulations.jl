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

    BLAS.set_num_threads(1)

    include(srcdir("main_functions.jl"))
    include(srcdir("aux_functions.jl"))

    mi_settings=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

    networks_paths_original=readdir(srcdir("empirical_networks", "original"); join=true)

    networks_names=readdir(srcdir("empirical_networks", "original"))
    networks_names=replace.(networks_names, r".csv.*" => "")

    networks_original=readdlm.(networks_paths_original, ',')

    binary!(networks=networks_original)

    adj_original=adjacency(incidence_matrices=networks_original)

    original_descriptors=networks_description(networks=networks_original, adj=adj_original, net_name=networks_names, matrix_type="incidence")

    c1=Dict(:net_info => original_descriptors, :m => mi_settings, :scenario => "original", :α => 0.2, :ϕ => 0.2, :ρ => 0.2, :θ_max =>5.0, :tmax => 10000, :ϵ => 1e-5, :nsim => 1000)
    c_list1=dict_list(c1)

    c2=Dict(:net_info => original_descriptors[[61,136]], :m => mi_settings, :scenario => "original", :α => 0.2, :ϕ => 0.2, :ρ => 0.2, :θ_max =>5.0, :tmax => 10000, :ϵ => 1e-5, :nsim => 1000)
    c_list2=dict_list(c2)

end

pmap(net_multisim, c_list1)
pmap(net_multisim, c_list2)

findall(x->x=="M_PL_057", networks_names)

original_descriptors[136]
original_descriptors[61]
using Random
using Distributions
d=Normal(5.0, 1.0)