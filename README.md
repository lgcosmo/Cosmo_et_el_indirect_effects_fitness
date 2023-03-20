# Indirect effects shape species fitness in coevolved mutualistic networks

This README_Cosmo_et_al_indirect_effects_fitness.txt file was generated on 2023-03-14 by Leandro Giacobelli Cosmo

GENERAL INFORMATION

1. Title of Dataset: 

Cosmo_et_al_indirect_effects_fitness

2. Author information:

Leandro G. Cosmo: Programa de Pós-Graduação em Ecologia, Departamento de Ecologia, Instituto de Biociências, Universidade de São Paulo, São Paulo, Brazil.

Ana Paula A. Assis: Departamento de Genética e Biologia Evolutiva, Instituto de Biociências, Universidade de São Paulo, São Paulo, Brazil.

Marcus A. Aguiar: Instituto de Física ‘Gleb Wataghin’, Universidade Estadual de Campinas, Campinas, Brazil

Mathias M. Pires: Departamento de Biologia Animal, Instituto de Biologia, Universidade Estadual de Campinas, Campinas, Brazil.

Alfredo Valido: Island Ecology and Evolution Research Group, Institute of Natural Products and Agrobiology (IPNA-CSIC), C/Astrofísico Francisco Sánchez 3, 38206 La Laguna, Tenerife (Canary Islands, Spain).

Pedro Jordano: Estación Biológica de Doñana, CSIC, av. Americo Vespucio 26, 41092, Sevilla, Spain; Departamento de Biología Vegetal y Ecología, Universidad de Sevilla, Sevilla, Spain.

John N. Thompson: Department of Ecology and Evolutionary Biology, University of California, Santa Cruz, 95064, California, USA.

Jordi Bascompte: Department of Evolutionary Biology and Environmental Studies, University of Zurich, Winterthurerstrasse 190, Zurich, CH-8057, Switzerland.

Paulo R. Guimaraes Jr.: Departamento de Ecologia, Instituto de Biociências, Universidade de São Paulo, São Paulo, Brazil.

Corresponding author: Leandro G. Cosmo, E-Mail: lgcosmo@usp.br

2. Information about funding sources that supported the collection of the data:

This study was financed in part by the Coordenação de Aperfeiçoamento de Pessoal de Nível Superior – Brasil (CAPES) – Finance Code 001. LGC is funded by a São Paulo Research Foundation PhD scholarship (FAPESP; grants # 2019/22146-3 and #2022/07939-0). PRG is funded by CNPq (307134/2017-2), FAPESP (São Paulo Research Foundation; grant: # 2018/14809-0), and the Royal Society, London (CHL/R1/180156). APA was supported by Sao Paulo Research Foundation (FAPESP grant #2016/14277-2). MMP is funded by Sao Paulo Research Foundation (FAPESP grant #2019/25478-7). MAA is funded by São Paulo Research Foundation (FAPESP, grants #2016/01343‐7 and #2021/14335-0 - ICTP‐SAIFR), and by Brazil's Council for Scientific and Technological Development (CNPq, grant #301082/2019‐7). A.V. is supported by grant PGC2018-099772-B-100 from Spanish Ministry of Science, Innovation, and Universities (AEI). JB’s research is supported by the Swiss National Science Foundation (grant 310030_197201).

DATA & FILE OVERVIEW

1. File List: 

Data files:

Folder empirical_networs/original: files of the networks used to parameterize the model in the numerical simulations used to build figures 1-2.
Folder empirical_networs/invasion_apis: files of the networks used to parameterize the model in the numerical simulations used to build figures 3-4.

Scripts/Source functions:

main_functions.jl
aux_functions.jl
main_functions_invasion.jl
aux_functions_invasion.jl
main_simulations.jl
main_simulations_invasion.jl

DATA-SPECIFIC INFORMATION:

main_functions.jl/aux_functions.jl: Julia functions used to run the model numerical simulations.
main_functions_invasion.jl/aux_functions_invasion.jl: Julia functions used to run the model numerical simulations of the scenarios of invasion by A. mellifera

main_simulations.jl: script to reproduce the numerical simulations of the model used in the main text.
main_simulations_invasion.jl: script to reproduce the numerical simulations of scenario of the invasion by A. mellifera used in the main text.

USAGE INSTRUCTIONS:

This code base is using the Julia Language and [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> Cosmo_et_al_indirect_effects_fitness

It is authored by Cosmo et al.

To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and everything should work out of the box, including correctly finding local paths. 

After installing everything, run the script "main_simulations.jl"  or "main_simulations_invasion.jl" located at the "scripts" folder.
For specific details about the scripts and the data frames that the simulations output, please refer to the instructions at each of the scripts mentioned above.
