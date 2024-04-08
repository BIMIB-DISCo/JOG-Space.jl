# Script3
In this example we show how it is possible to use functions of J-SPACE to perform a parameter-scan of a parameter of a molecular evolution in the same phylogenetic tree and how the generated sequences of the cells of the tree changed(e.g., mutation rate average, indel rate, etc.).
In this example we use in particular 5 functions:
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
  - **tree_phylo**, a directed Int64 metagraph with Float64 weights, object MetaGraph{Int64, Float64}. 
  - **newick**, is the format Newick of tree_phyl, object *PhyloNetworks.HybridNetwork, Rooted Network* in package PhyloNetworks.  

5. ***experiment_noISA_sign***(tree_phylo, length_genome, sub_model, mut_rate_avg, indel_rate, indel_size, seed, set_mut, Lavelette_par, used_sign, vector_change_points, vector_activities, ratio_background_signature)  
This is the function that simulates the molecular evolution on the phylogenetic tree using signature.  
- `tree_phylo`, previously described.
- `length_genome`.  Length of the ancestral genome.
- `sub_model`. A string, variable that select the subistituion model, possible value for signature ["SBS-37","SBS-38"].  
- `mut_rate_avg`. A real number, the average mutational rate per trinucleotide and unit of time.
- `indel_rate`. Rate of indel per site and per unit of time. To esclude indel in the simulation put this parameter to 0.0 .  
- `indel_size`. Integer number, maximum size of indel.
- `seed`, previously described.  
- `set_mut`, previously described.  
- `Lavelette_par`. Real number, the parameter of  the Lavelette distribution for the size of indels.  
- `used_sign`. An array of strings. The list of the labels of the used signatures (e.g., `used_sign = ["SBS1","SBS4","SBS16"]`)
- `vector_change_points`. An array of real number . It contains the list of the time in the change points (e.g., one change-point at time  50 should be specified `vector_change_points = [0.0, 50.0]`). If no change point is desired (constant signature activity) one should use `vector_change_points = [0.0]`
- `vector_activities`. A vector of vectors. For each elements of `vector_change_points` it is necessary insert an array with the values of the activities for each signature (e.g., if `vector_change_points` has two elements `vector_activities = [[0.7,0.2,0.1], [0.0,0.3,0.7]]
`)
- `ratio_background_signature`. A real number between 1.0 and 0.0. If 1.0 all mutation will be due to mutational signatures, if 0 all ther mutations will be due to the background  process.  
This function has the following outputs:
  - **g_seq**, it is the generated ancestral genome or given. Object LongDNASeq in package *BioSequences*.  
  - **fasta_sample**, a Vector{Any} containing the Fasta file associated with the leafs of the phylogenetic tree. Each element is an object LongDNASeq.  
  - **tree_phylo**, is the same tree given as input but with new metadata: the Fasta file associated with the nodes.
  - **mutation_df**, a data frame containing driver and neutral mutations and in which leaves they are present.  

The script  work as follow:
- Firstly, we use the command 'using' will load the module or package and make its exported names available for direct use.  
In particular, for this example, we load :   
  1) J_Space: for principal function  
  2) Random, PhyloNetworks, Graph, MetaGraph: to auxialiry function  
  3) CSV, Tables and DataFrames, FASTX: to save data  
- We Initialize the seed of simulation with function MersenneTwister(Integer Number)
- We create the regular graph with  'g_meta = spatial_graph(200, 200, seed, dim = 3, n_cell=1)' 
- We simulate the dynamics with 'df, G, n_cell_alive, set_mut, Gs_conf, CA_subpop, α_subpop = simulate_evolution(g_meta,200.0,0.3,0.01,0.01,0.00001,0.2,0.1,"contact",seed)'
- We do a sampling with 'matrix_R, tree_mut = sampling_phylogentic_relation(G,"Random",df,100,set_mut,seed,1)'
- We generate the Newick file and phylogenetic tree with 'tree_red, net = create_tree(matrix_R, true, 200.0)'
- In this case we want to run the same molecular evolution using different value of the mutation rate average of the wild type. We set the range of the values into two variables.(i.e. mut_rate_avg_min and mut_rate_avg_max)  
- We choose the number of studied  values in the range  (n_rep)

- Create a range of values with function range by specifying the chosen values(i.e. values_of_mut_rate_avg)
- With a for loop we excute a simulation per point of the range *n_rep* (pay attention to indentation)  
- Inside the loop, for each iteration we performs:  
  1) a molecular evolution using singature with 'g_seq, fasta_samples, Tree, mutations_tot = experiment_noISA_sign(tree_red, 10000, "SBS-37", 0.001, 0.0, 100, seed, set_mut, 0.5, ["SBS6","SBS22"], [0.0], [0.5 0.5], 0.8)'  
  2) saving the Fasta file into variable *fasta_sample*
- Finally, we save all outputfiles that we have obtained from called function. 
