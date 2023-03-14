# This script runs the numerical simulations for the invasion by A. mellifera of the model used in Cosmo et al. 2023 (Nature)
# It uses multiple threads of the computer processor (addprocs function)

# Usage instructions:
# Set the number of desired x threads to use for the numerical simulations (addprocs(x)).
# Set the parameter values of the model to be used in the numerical simulations. Their description are as follows:
    # m: Overall strength of mutualistic interactions
    # α: Sensitivity of the function qij to differences in the trait values between species.
    # ϕ: Compound parameter ϕ=σρ where σ is the additive genetic variance of the trait z
    # ρ: Sensitivity of the adaptive landscape of species to the traits that mediate mutualistic interactions.
    # θ_max: Upper bound of the uniform distribution, U[0.0, θ_max], to be used to sample the environmental optima of species in the beginning of each simulation.
    # tmax: Maximum amount of generations to run the coevolutionary model.
    # ϵ: Difference in trait values between successive generations so that species traits are considered to be at equilibrium.
    # nsim: Number of simulations to run for each empirical network
# Run the code below

# The function saves a dataframe containing the desired number of simulations for a particular combination of parameter values/empirical network
# This data frame contains the following variables:

#network: ID of the network
#scenario: scenario of the simulation (normal or invasion)
#sp_id: ID of each species within the network
#type: whether species are plants or animals
#position: core or peripheral position of species in the network
#degree: number of mutualistic partners of the species after the invasion by A. mellifera
#degree_native: number of mutualistic partners of the species before the invasion by A. mellifera
#theta: environmental optimum of the species
#z0: initial trait values
#zeq: trait values at the coevolutionary equilibrium after the invasion
#w0: initial fitness
#weq_preinvasion: fitness at coevolutionary equilibrium before the invasion by A. mellifera
#weq_postinvasion: fitness at coevolutionary equilibrium after the invasion by A. mellifera
#wdif: differences in fitness values before and after the invasion by A. mellifera
#indirect_effects_post: contribution of indirect effects to the species after the invasion by A. mellifera
#indirect_effects_pre: contribution of indirect effects to the species before the invasion by A. mellifera
#teq: time until equilibrium before the invasion by A. mellifera
#teq_inv: time of the invasion by A. mellifera
#teq_afterinv: time until equilibrium after the invasion by A. mellifera
#m: parameter controlling the proportional contribution of mutualisms as selective pressures
#alpha: parameter α of the model
#rho: parameter ρ of the model
#phi: compound parameter ϕ=σρ where σ is the additive genetic variance of the trait z
#sim: ID of the simulation

using Distributed

addprocs(20)

@everywhere using DrWatson
@everywhere @quickactivate "Cosmo_et_al_indirect_effects_fitness"

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

    include(srcdir("main_functions_invasion.jl")) #Loading required functions
    include(srcdir("aux_functions_invasion.jl")) #Loading required functions

    mi_settings=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9] # Setting values of m to be used in the numerical simulations

    networks_paths_original=readdir(srcdir("empirical_networks", "apis_invasion"); join=true) #Importing networks to simulate the invasion scenario

    networks_names=readdir(srcdir("empirical_networks", "apis_invasion")) #Importing networks to simulate the invasion scenario
    networks_names=replace.(networks_names, r".csv.*" => "") #Importing networks to simulate the invasion scenario
    networks_names=networks_names[1:73] #Importing networks to simulate the invasion scenario

    files_invasion=readdlm.(networks_paths_original, ',') #Importing networks to simulate the invasion scenario
    networks_invasion=files_invasion[1:73] #Importing networks to simulate the invasion scenario
    apis_adj=files_invasion[74] #Loading file with the position of A. mellifera in each network (adjacency matrix)
    apis_adj=Int64.(apis_adj) #Loading file with the position of A. mellifera in each network (adjacency matrix)
    apis_inc=files_invasion[75] #Loading file with the position of A. mellifera in each network (incidence matrix)
    apis_inc=Int64.(apis_inc) #Loading file with the position of A. mellifera in each network (incidence matrix)

    binary!(networks=networks_invasion) # Setting entries in the empirical networks to be binary (aij=1 if species interact and 0 otherwise)
    networks_native=native(networks_invasion, apis_inc) #Function to remove A. mellifera from the networks to simulate the native networks

    invasion_descriptors=networks_description(networks=networks_invasion, networks_native=networks_native, net_name=networks_names, matrix_type="incidence", invader_position=apis_adj)

    c=Dict(:net_info => invasion_descriptors, :m => mi_settings, :scenario => "invasion", :α => 0.2, :ϕ => 0.2, :ρ => 0.2, :tmax => 10000, :ϵ => 1e-5, :nsim => 1000)
    c_list=dict_list(c)

end

pmap(net_multisim, c_list)


