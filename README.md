# J-SPACE 

## POSSIBLE PROBLEMS

Since the library GLMakie uses the GPU, a possible error on virtual machines is the following:
```
 LoadError: InitError: Exception[GLFW.GLFWError(GLFW.PLATFORM_ERROR, "X11: Failed to open display localhost:20.0"), ErrorException("glfwInit failed")]
```

## INTRODUCTION 
J-SPACE is a Julia package to simulate the spatial growth and the genomic evolution  of a cell population and the experiment of sequencing the genome of the sampled cells.
Firstly, the software simulates the spatial dynamics of the cells as a continuous-time multi-type birth-death stochastic process on a graph employing different rules of interaction and an optimised Gillespie algorithm. 
After mimicking a spatial sampling of the tumour cells, J-SPACE  returns the phylogenetic tree of the sample and simulates molecular evolution of the genome under the infinite-site models or a set of different substitution models. Ther is also  the possibility of include indels. Finally, employing ART, J-SPACE generates the  synthetic single-end, paired-/mate-pair end reads of the next-generation sequencing platforms.


 ### Spatial clonal dynamics
 In J-SPACE the dynamics of the spatio-temporal evolution of a tumour is modelled by a stochastic multi-type Birth-Death process over an
arbitrary graph. J-SPACE could generate by itself a 2D or 3D regular lattice. In addition, it is possible to give as input any graph as an adjacency matrix ( an example of the format needed for this matrix is given in path "Example_adj_matrix", it must is symmetric and with only values(0,1) separeted from space between them).   
In this part, the user can tune the birth rate of wild type cells, the death rate of the cells, the rate of migration of cells (not tested), the rule of contact between cells (to simulate different mechanical interactions), the probability to develop a driver mutation per division, and the average birth rate advantage of a driver mutation.
Additionally, there is the possibility to performe an excision by specifying the timing and the ratio of cells that will die of the event associated. 
If the user want a specific the clonal dynamics (i.e., `Tree_Driver_Configure = 1`), it is possible to indicate the edge list representing the mutational tree of drivers and the path where this file is supplied in txt (the parameter of the config file `edgelist_treedriver` ). In this case the user should also specify the birth rate of each subpopulation and the path where this file is supplied in txt (the parameter of the config file  `driver_birth_rates`  ). For example, a linear tree with tree drivers is described by the following:

Driver_1 Driver_2

Driver_2 Driver_3

Driver_3 Driver_4

The file with the birth rate must have the following format:

Driver_1 0.2

Driver_2 0.4

Driver_3 0.5

Driver_4 0.6

In this case J-SPACE accept only the events that respect the mutational tree given.

Note that every rate inserted in J-SPACE must have the same unit of time both for the spatial dynamics and the molecular evolution.  

 ### Molecular evolution
 
J-SPACE simulate the evolution of the sequence of the sample after the simulation of the clonal dynamics. The user can sample the whole population or a subset of it, and the J-SPACE evaluate the phylogenetic tree of the samples. This GT tree is returned as a Newick file in the folder specified by the variable "/path_to_save_files/" .  

The molecular evolution of an ancestral genome  (which can be given by the user as FASTA file or generated randomly) is simulated along the sampled tree via the Doob-Gillespie algorithm.
The user can use an infinite-site model to have fast simulations of situations where the genome is long, the mutational rate is very low (e.g.,<10^-8 substitution for unit of time per site), and the total simulated time is long.

In the case of finite-site models, J-SPACE takes as input the matrix of instantaneous rates for different substitution models: JC69, F81, K80, HKY85, TN93, and K81.
We suppose that the indels have a size distributed as  a Lavalette distribution.

The user can also generate a custom time-dependent substitution model based on a linear combination of the Mutational signatures of the COSMIC database. In this case the user should provide the list of labels of the desired SBS signature in the COSMIC database (https://cancer.sanger.ac.uk/signatures/) (e.g.,  `used_sign = ["SBS1","SBS4","SBS16"]`
 ), the list of change points (e.g., one change-point at time  50 should be specified `vector_change_points = [0.0, 50.0]`), the values of the activations for each signature in each of the time span defined by the change-points (e.g., `vector_activities = [[0.7,0.2,0.1], [0.0,0.3,0.7]]`)
) and the ratio of mutation due to the background uniform process or due to the mutational signatures (e.g., `ratio_background_signature = 0.8`).
 Note that using finite-site for long genomes come at the cost of computational performance.
 After this computation, the sequences of the samples (i.e., the leafs of the phylogenetic tree) are returned as FASTA file in the folder  /"path_to_save_files"/Fasta output.
 
 

 ### Sequencing experiment

To simulate the reads of a sequencing experiment J-SPACE calls ART (https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm). The user can use a configuration file to specify the error model (for Illumina platforms), the number of reads, the length of the reads, and if the experiment uses single-end paired-end/mate-pair reads.
In addition, in the configuration file, there is the option to insert custom "calls" for ART  with the possibility to use it in any possible configuration.
If the user simulate the experiment, J-SPACE returns for each cell, the simulated reads as FASTQ file, the alignment map of the reads over the genome of the sampled cells in SAM and/or ALN format. 
If the infinite-site model is used,  it is possible to obtain the VCF file directly without simulating the reads with ART.
## REQUIRED  SOFTWARE AND PACKAGE
- Operating system: Linux, Windows.
- Julia (https://julialang.org/) Version 1.6.3 or higher.
- ART, install it from https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm.

## INSTALLATION OF J-SPACE

J-SPACE can be downloaded from Github. 
First, it is necessary to install the Julia from https://julialang.org/.   
Next, the user need to copy the project folder in the chosen working directory. To install J-SPACE follow the steps:
1. Using REPL or the COMMAND LINE move to the working directory.  
2. If you use the COMMAND LINE, to start a Julia session run the command:

```julia
> julia
```

3. To enter in the Pkg REPL  type `]`

```julia
(@v1.10) pkg>
```

4. Type the command

```julia
(@v1.10) pkg> activate .
```

5. To activate the J-SPACE project, type

```julia
(@v1.10) pkg> instantiate
```

## RUN J-SPACE

### RUN A SINGLE SIMULATION
The parameters and the configuration of the simulation are managment by the user by modifing the files "Parameters.toml" and "Config.toml" (the name of the file is not mandatory),  that are detailed in the next sections.  
To run a simulation of J-SPACE using the ".toml" file for the paramet follow the following step:

1. Load the J-SPACE package using
```julia
> using J_Space
```

2. Start the simulation
```julia
> Start_J_Space("Parameters.toml", "Config.toml")
```

NOTE: the simulation does not start if in the working folder are absent the two .toml files. 

### RUN THE EXAMPLES
To run the examples, in the main folder of J-SPACE, from command line digit

```julia
> julia --project=.  ./test/data/Experiments/Experiment_2D/experiment_2D.jl  
```

or

```julia
> julia --project=.  ./test/data/Experiments/Experiment_3D/experiment_3D.jl
```

#### Run the variant calling pipeline
##### Necessary package
- download the file `environment_j_space.yml`
- download the file `j_space_pipeline.sh`
- conda https://docs.conda.io/en/latest/
- gatk https://storage.googleapis.com/gatk-software/package-archive/gatk/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2

##### Run the variant calling
First move where is file `environment_j_space.yml`

Change field name and prefix into `environment_j_space.yml`

Open the conda environment
```julia
> conda env create -f environment_j_space.yml --prefix "path to the enviroment directory"
```

Activate the conda environment
```julia
> conda activate "path to the enviroment directory"  
```

Then register to gatk (Not necessary in the same working folder)
```julia
> gatk3-register "path to the gatk directory"
```

Move into working folder where you have `j_space_pipeline.sh` 
Run pipeline
```julia
> ./j_space_pipeline.sh "path/to/reference/" "path/to/FastaQ" "path/working/directory"
```

- the first, "path/to/reference/", indicate where is the reference fasta, must specificated file.fasta
- the second, "path/to/FastaQ", indicate where are the file FastaQ to analyzed, must specificated the folder
- the third, "path/working/directory", indicate your working directory, where will save output files

NOTE: the paths are absolute paths 


## OUTPUTS OF J-SPACE
 J-SPACE provides the following outputs. 

  - The state of the lattice at any time of the simulation (as plot).
  - The Ground Truth (GT) genotype of the sampled cells as FASTA files.
  - The GT phylogenetic  tree  of the sampled cells in Newick format.
  - The plot of the mutational tree of the driver mutations (if present).
  - The list of the driver mutation with their birth rate advantage as .csv (if present).
  - The simulated next-generation sequencing reads as FASTQ files in the folder "/"path_to_save_files"/Fasta output/sample_#", where sample_# is the #-th sample.
  - The alignment file, which maps the noisy reads on the sequences of the sampled cells both in formats SAM and ALN.
  - The alignment file, without noise in format SAM.
  - The list of the mutations for each sample as csv.


## THE CONFIGURATION FILE OF J-SPACE
In the file "Config.toml" the user can manage the configuration of J-SPACE. This file is useful to choose the path where save the files, the desired plots and output files.
We provide an example in the main folder of J-SPACE.

The following are the paramenters of such file:

- `seed`. Integer, the seed of the simulation.
- `generate_reference`. Integer, If 0 the user should inser the ancestral genome in fasta format in the folder "path_reference". If 1 J-SPACE generate a random sequence.
- `path_reference`. A string,  the path of the reference given by user.
- `path_to_save_files`.  A string, path of the folder where the output files should be saved.
- `path_to_save_plot`. A string,  path of the folder where the output plots should be saved.
- `Generate_graph`. Integer, if 0 the user should inser the graph of the dynamics as an adjacency matrix. If 1 J-SPACE generate a regular lattice, the paramenters of such lattice are specified in the file "Parameters.toml".
- `Path_to_Graph`.  Path to the adjacency matrix of the graph if imported.
- `Tree_Driver_Configure`. Integer number. if 1 J-SPACE takes as input a driver mutational tree. With 0 J-SPACE generates the tree randomly. 
- `Tree_Newick`. Integer, if 1  the phylogenentic tree of the cells is saved as newick file.
- `Final_configuration`. Integer,  if 1 return the plot of the final configuration of the lattice.
- `Driver_list`. Integer,  if 1 returns the list of the driver mutations.
- `Driver_Tree`. Integer, if 1 returns the plot tree of the driver mutations.
- `Single_cell_fasta`. Integer, if 1 save the fasta of the GT sequences of the sampled cells.
- `Single_cell_noise`. Integer, if 0  the sequencing experiment is not performed, if 1 it returns the FASTQ  files of the reads and SAM files with noise  , if 2 it returns the FASTQ  files of the reads, the SAM files with noise and without noise.
- `Alignment`. Integer, if 1 it returns the GT aligment file in ALN format.
- `VAF_GT`. Integer, if 1 it returns the VAF of the sampled cells  (working only if `type_isa = 1`).
- `Bulk_noise`. Integer, if 1 it returns an approximate bulk experiment not using ART (working only if `type_isa = 1`).
- `Time_of_sampling`. Array of times in an incresing order (e.g., [10.0 , 20.3 , 50.2]). J-SPACE perform the plot of the state of the lattice in that times.
- `Graph_configuration`. Integer,  if 1  returns the plot of the state of the lattice.

#### If `Tree_Driver_Configure = 1` 
- `edgelist_treedriver`. A string. The path to the txt file containg the edge list of the driver mutational tree.
- `driver_birth_rates`. A string. The path to the txt file containg the edge list of the drivers birth rates.

## THE PARAMETERS FILE OF J-SPACE
In the file "Parameters.toml" the user will find all the paramenters of the dynamics, molecular evolution and experiment. 

### Parameters of the generation of the lattice
- `row`. Integer number of the rows of the regular lattice, not used if a graph is imported
- `col`. Integer number of the columns of the regular lattice, not used if a graph is imported
- `dim`. Integer number of the height of the regular lattice, not used if a graph is imported. 
- `N_starting_cells`. Integer  the number of starting cells.

### Parameters of the clonal spatial dynamics
- `Model`. String, the possible values are `["contact", "voter", "hvoter"]` they are the possible different  Models of interaction.
- `Max_time`. Real number, it is the maximum time of the simulation 
- `rate_birth`. Real number,  Birth rate per cell per unit of times (not used if `Tree_Driver_Configure = 1`)
- `rate_death`. Real number,  Death rate per cell per unit of times
- `rate_migration`. Real number,  Migration rate per cell per unit of times
- `drive_mut_rate`. Real number, probability of generate a new driver (i.e., subclones) after a division event.
- `average_driver_mut_rate`. Real number,  average birth rate advantage a driver mutation 
- `std_driver_mut_rate`. Real number, standard deviation of the birth rate advantage of a driver mutation 
- `t_bottleneck`. An array of real number, it contains the list of times of the bottlenecks.
- `ratio_bottleneck`. An array of real number between 0 and 1, it contains the ratio of cells that will die at times `t_bottleneck`
### Parameters of the sampling
- `Random_sampling`. Integer, if 1  Random sampling, if 0 circular/spherical sampling 
- `num_cell`. Integer, number of sampled cells

#### If `Random_sampling = 0`
-  `pos_center`. Integer number, it is the center of the sampling
- `radius_sampling`. Integer number, it is size radius of the sampling

### Parameters of the molecular evolution
- `length_genome`.  Length of the ancestral genome. Used  if ancestral genome is not given. If the ancestral genome is given as fasta file this line is ignored.
- `prob_base`. An array of three real numbers (e.g., [0.3,0.2,0.2]).   the proportion of A, C, and G in the genome randomly generated. The proportion of T is calculated by normalization.
- `type_isa`. Integer number, if 1 J-SPACE use the ISA to simulate the molecular evolution, if 0 it is necessary to specify the substituion model below.

#### if `type_isa = 1`
- `neut_mut_rate`. Rate of mutation per site and per unit of time

#### if `type_isa = 0`
- `approx_snv_indel`. Integer number, if 0 SNV and INDEL are computed togheter, if 1 SNV and INDEL are computed separately. The first choice is very slow, and we recommend the use only for small and short trees.
- `sub_model`. A string, variable that select the subistituion model, possible value [ "JC69","F81","K80", "HKY","TrN93ef","TrN","K81","K81uf"].
- `indel_size`. Integer number, maximum size of indel .
- `Lavelette_par`. Real number, the parameter of  the Lavelette distribution for the size of indels
- `indel_rate`. Rate of indel per site and per unit of time. To esclude indel in the simulation put this parameter to 0.0 .
- `params`. Rates of the substitution models  in units of time of the simulation.  Below, some examples:
 if `sub_model= "JC69" `-> `params = [{"alpha" = 0.5}]`.   
 if `sub_model= "F81" ` -> `params = [{"alpha" = 0.5}]`.    
 if `sub_models= "K80"  `-> `params = [{"alpha" = 0.5, "beta" = 0.3}]`.    
 if `sub_models= "HKY" ` -> `params = [{"alpha" = 0.5, "beta" = 0.3}]`.    
 if `sub_models = "TrN93ef"`  -> `params = [{"alpha" = 0.5, "alpha2"=0.1,"beta" = 0.3}]`.     
 if `sub_models = "TrN" `-> `params = [{"alpha" = 0.5, "alpha2"=0.1,"beta" = 0.3}]`.     
 if `sub_models = "K81" ` -> `params = [{"alpha" = 0.5, "beta"=0.1,"beta2" = 0.3}]`.    
 if `sub_models = "K81uf"`  -> `params = [{"alpha" = 0.5, "beta"=0.1,"beta2" = 0.3}]`.

#### if sub_model = "SBS-37" or "SBS-38" J-Space uses a mutational signature based substitution model

- `mut_rate_avg`. A real number, the average mutational rate per trinucleotide and unit of time. 
- `used_sign`. An array of strings. The list of the labels of the used signatures (e.g., `used_sign = ["SBS1","SBS4","SBS16"]
`) 
- `vector_change_points`. An array of real number . It contains the list of the time in the change points (e.g., one change-point at time  50 should be specified 
- -`vector_change_points = [0.0, 50.0]`). If no change point is desired (constant signature activity) one should use `vector_change_points = [0.0]`
- `vector_activities`. A vector of vectors. For each elements of `vector_change_points` it is necessary insert an array with the values of the activities for each signature (e.g., if `vector_change_points` has two elements `vector_activities = [[0.7,0.2,0.1], [0.0,0.3,0.7]]
`)
- `ratio_background_signature`. A real number between 1.0 and 0.0. If 1.0 all mutation will be due to mutational signatures, if 0 all ther mutations will be due to the background  process.

### Parameters of the bulk experiment (approximate version, very fast does not need ART, but it generates the VAF, not the reads. Working only if `type_isa = 1`)
- `coverage`. Real number, average coverage of the simulate bulk experiment.
 - `FP`. Real number,  false positive rate.
- `FN`. Real number, false negative rate.

### Parameters of the sequencing experiment (ART)
- `command`. A string, if the user want to do custom calls of ART, e.g., `command = "art_illumina -ss HS25 -sam -i reference.fa -l 150 -f 10 -o single_dat"`.

#### If `command = ""`  for Illumina sequencing system is possible to compile the following parameters
- `profile`. String, the name of Illumina sequencing system of the built-in profile used for simulation, e.g.,`profile = "HS25"`.
- `len_read`. Integer, the length of reads to be simulated.
- `tot_num_reads`. Integer, number of reads/read pairs to be generated per sequence.
- `outfile_prefix`. String, the prefix of output filename.
- `paired_end`. Integer number,  0  indicate a single-end read simulation, if 1 a paried_end simulation is performed. If `paired_end = 1`, it is necessary to set  `mate_pair = 0`.
   NOTE: if paired end is equal to 1 you will find 2 FASTQ files for sample.             		    
- `mate_pair`. Integer, if  1 mate-pair read simulation. If `mate_pair = 1`, it is necessary to set  `paired_end = 0`.
   Art will automatically switch to a mate-pair simulation if the given mean fragment size >= 2000.
   
#### if `paired_end = 1`,  are required the following 
- `mean_fragsize`. Integer, the mean size of DNA/RNA fragments for paired-end simulations.
- `std_fragsize`. Integer, the standard deviation of DNA/RNA fragment size for paired-end simulations.

For all paramenters of ART  please see: https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm


## LICENSE
See the file `COPYING` for license information.


## CONTACTS
Please feel free to contact us if you have problems running our tool at fabrizio.angaroni@unimib.it and a.guidi@campus.unimib.it .
