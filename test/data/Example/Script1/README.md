# Script 1
In this example we show how it is possible to use  functions of J-SPACE to performa a parameter-scan of a dynamical parameter (e.g., birth rate, death rate, etc.).  
In this example we use in particular 4 functions:
1. ***spatial_graph***(row, n_col, seed, dim = dim, n_cell=N_starting_cells)  
This is the function that generates the meta graph where the dynamics occurs.
- `row`. Integer number of the rows of the regular lattice.  
- `col`. Integer number of the columns of the regular lattice.   
- `seed`. Integer, the seed of the simulation.     
- `dim`. Integer number of the height of the regular lattice.   
- `N_starting_cells`. Integer  the number of starting cells.    
This function has the following output: **g_meta**, an undirected INT64 meta graph with Float64 weights. Object MetaGraph{Int64, Float64} in package *MetaGraphs*

2. ***simulate_evolution***(g_meta,Max_time,rate_birth, rate_death, rate_migration, drive_mut_rate, average_driver_mut_rate, std_driver_mut_rate,  Model, seed)  
This is the function that simulates the dynamics of the evolution cancer.    
- `g_meta`. A meta graph that represents the space of dynamics.  
- `Max_time`. Real number, it is the maximum time of the simulation.  
- `rate_birth`. Real number,  Birth rate per cell per unit of times.
- `rate_death`. Real number,  Death rate per cell per unit of times.  
- `rate_migration`. Real number,  Migration rate per cell per unit of time.  
- `drive_mut_rate`. Real number, probability of generate a new driver (i.e., subclones) after a division event.  
- `average_driver_mut_rate`. Real number,  average birth rate advantage a driver mutation.  
- `std_driver_mut_rate`. Real number, standard deviation of the birth rate advantage of a driver mutation.  
- `Model`. String, the possible values are `["contact", "voter", "hvoter"]` they are the possible different Models of interaction.  
- `seed`, previously described.  
This function has the following outputs: 
  - **df**, a data frame containing all events done during the dynamics. Object DataFrame in package *DataFrames*  
  - **G**, an undirected INT64 meta graph with Float64 weights, is the final configuration of the graph at the end of the dynamics  
  - **n_cell_alive**, a Vector{Any} containing the number of the cells alive at each recorded event. 
  - **set_mut**, a Vector{Any} containing all driver mutations occurred.
  - **Gs_conf**, a Vector{Any} containing the configuration of the graph at a specific time point(it will is emply).  
  - **CA_subpop**, a Vector of Vectors containig the number of cells divided by subpopulation.  
  - **α_subpop**, a Vector{Any} containing the birth rate of all subpopolation occurred.  

3. ***sampling_phylogentic_relation***(G,mode_sampling,df,num_cell,set_mut,seed,Driver_Tree)
This is the function that performs a sampling and it generates the relationals matrix of the sampled cells and tree driver mutations 
- `G`, previously described.  
- `mode_sampling`. String, if "Random" performs Random sampling, if "Neighbourhood" circular/spherical sampling.   
- `df`, previously described.  
- `num_cell`. Integer, number of sampled cells.  
- `set_mut`, previously described. 
- `seed`, previously described.
- `Driver_Tree`. Integer, if 1 returns the plot tree of the driver mutations.  
This function has the following outputs:
  - **matrix_R**, a data frame containing the relationals of the sampled cells "Father-Son", time of birth of the child, subpopolation and type of event.  
  - **tree_mutations**, a directed Int64 metagraph with Float64 weights, object MetaGraph{Int64, Float64}.  

4. ***create_tree***(matrix_R, Tree_Newick, Max_time)
This is the function that generates the phylogenetic tree.  
- `matrix_R`, previously described.  
- `Tree_Newick`. Boolean, if true the phylogenentic tree of the cells is saved as newick file.
- `Max_time`. It is the same in simulate_evolution.  
This function has the following outputs:
  - **tree_phyl**, a directed Int64 metagraph with Float64 weights, object MetaGraph{Int64, Float64}. 
  - **newick**, is the format Newick of tree_phyl, object *PhyloNetworks.HybridNetwork, Rooted Network* in package PhyloNetworks.    

The script  work as follow:
- Firstly, we use the command 'using' will load the module or package and make its exported names available for direct use.  
In particular, for this example, we load :   
  1) J_Space: for principal function  
  2) Random, PhyloNetworks, Graph, MetaGraph: to auxialiry function  
  3) CSV, Tables and DataFrames: to save data  

- We Initialize the seed of simulation with function MersenneTwister(Integer Number)
- In this case we want to run the same simulation using different value of the birth rate of the wild type. We set the range of the values into two variables.(i.e. birth_rate_min and birth_rate_max)  
- We choose the number of studied  values in the range  (n_rep)

- Create a range of values with function range by specifying the chosen values(i.e. values_of_birth_rate)

- We initialize all variables that we want to save (e.g., newick files, Newick_finals = [])
- We create the regular graph with  'g_meta = spatial_graph(200, 200, seed, dim = 3, n_cell=1)'  
- With a for loop we excute a simulation per point of the range (pay attention to indentation)
- Inside the loop, for each iteration we performs:
  1) the dynamics with 'df, G, n_cell_alive, set_mut, Gs_conf, CA_subpop, α_subpop = simulate_evolution(g_meta,200.0,0.3,0.01,0.01,0.00001,0.2,0.1,"contact",seed)'
  2) an sampling with 'matrix_R, tree_mut = sampling_phylogentic_relation(G,"Random",df,100,set_mut,seed,1)'
  3) a generation of a Newick file with 'tree_red, net = create_tree(matrix_R, true, 200.0)'   
- Finally, we save all outputfiles that we have obtained from called function.
