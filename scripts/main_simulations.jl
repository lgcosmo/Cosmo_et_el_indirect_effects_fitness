# This script runs the numerical simulations of the model used in Cosmo et al. 2023 (Nature)
# It uses multiple threads of the computer processor (addprocs function)

# Usage instructions:
# Set the number of desired x threads to use for the numerical simulations (addprocs(x)).
# Set the parameter values of the model to be used in the numerical simulations. Their description are as follows:
    # m: Overall strength of mutualistic interactions
    # α: Sensitivity of the function qij to differences in the trait values between species.
    # ϕ: Compound parameter ϕ=σρ
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
#degree: number of mutualistic partners of the species
#theta: environmental optimum of the species
#z0: initial trait values
#zeq: trait values at the coevolutionary equilibrium
#wmax: maximum theoretical fitness of the species
#w0: initial fitness
#weq: fitness at coevolutionary equilibrium
#indirect_effects: contribution of indirect effects to the species
#teq: time until equilibrium
#m: parameter controlling the proportional contribution of mutualisms as selective pressures
#alpha: parameter α of the model
#rho: parameter ρ of the model
#theta_max: upper bound of the uniform distribution used to sample environmental optima of species
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

    BLAS.set_num_threads(1)

    include(srcdir("main_functions.jl")) #Loading required functions
    include(srcdir("aux_functions.jl")) #Loading required functions

    mi_settings=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9] # Setting values of m to be used in the numerical simulations

    networks_paths=readdir(srcdir("empirical_networks", "original"); join=true) # Loading dataset of empirical networks

    networks_names=readdir(srcdir("empirical_networks", "original")) # Loading dataset of empirical networks
    networks_names=replace.(networks_names, r".csv.*" => "") # Loading dataset of empirical networks

    networks_original=readdlm.(networks_paths, ',') # Loading dataset of empirical networks


    binary!(networks=networks_original) # Setting entries in the empirical networks to be binary (aij=1 if species interact and 0 otherwise)

    adj_original=adjacency(incidence_matrices=networks_original) # Transforming the incidence matrices into adjacency matrices to be used in the numerical simulations

    original_descriptors=networks_description(networks=networks_original, adj=adj_original, net_name=networks_names, matrix_type="incidence") # Getting general descriptive measures of each empirical network

    c=Dict(:net_info => original_descriptors, :m => mi_settings, :scenario => "original", :α => 0.2, :ϕ => 0.2, :ρ => 0.2, :θ_max =>10.0, :tmax => 10000, :ϵ => 1e-5, :nsim => 1000) # Setting combinations of parameter values to run the simulations
    c_list=dict_list(c1) # Setting combinations of parameter values to run the simulations

end

pmap(net_multisim, c_list1) # Running the numerical simulations
