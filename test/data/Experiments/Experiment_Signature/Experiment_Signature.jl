using J_Space
using Random
using CSV, PhyloNetworks, Graphs, MetaGraphs, Tables, FASTX, DataFrames
using Distributed
using DelimitedFiles, LinearAlgebra

println("load library")

if Sys.iswindows()
    edge_list_path = "utility\\Tree_driver_example.txt"
    driv_adv_path = "utility\\driver_advantage.txt"
else
    edge_list_path = "utility/Tree_driver_example.txt"
    driv_adv_path = "utility/driver_advantage.txt"
end

seed = MersenneTwister(1234)

g_meta = spatial_graph(200, 200, seed, dim = 3 , n_cell=1)

df, G, n_cell_alive, set_mut, Gs_conf, CA_subpop, α_subpop = simulate_evolution(
                                                        g_meta,
                                                        200.0,#time
                                                        0.01,#rate_death
                                                        0.0,#rate_migration
                                                        0.01,#rate_mut
                                                        "hvoter", #model
                                                        edge_list_path,
                                                        driv_adv_path,
                                                        seed)


matrix_R, tree_mut = sampling_phylogentic_relation(G,
                                                   "Random",
                                                   df,
                                                   100,
                                                   set_mut,
                                                   seed,
                                                   0)

tree_red, net = create_tree(matrix_R, true, 200.0)
leafs = J_Space.get_leafs(tree_red)
#save
if Sys.iswindows()
    mkpath(".\\Experiments\\Experiment_Signature\\Fileoutput\\")
    writeTopology(net,
                  ".\\Experiments\\Experiment_Signature\\Fileoutput\\formatNewick")
    CSV.write(".\\Experiments\\Experiment_Signature\\Fileoutput\\Dinamica.csv",
                              df,
                              delim = ",")
    CSV.write(".\\Experiments\\Experiment_Signature\\Fileoutput\\CA_subpop.csv",
                Tables.table(CA_subpop),
                header=false)
else
    mkpath("./Experiments/Experiment_Signature/Fileoutput/")
    writeTopology(net,
                  "./Experiments/Experiment_Signature/Fileoutput/formatNewick")
    CSV.write("./Experiments/Experiment_Signature/Fileoutput/Dinamica.csv",
                              df,
                              delim = ",")
    CSV.write("./Experiments/Experiment_Signature/Fileoutput/CA_subpop.csv",
                Tables.table(CA_subpop),
                header=false)
end



#####   EXP1
seed = MersenneTwister(1234)
g_sign = copy(tree_red)
g_seq, fasta_samples, position_used, mutations_tot = J_Space.experiment_noISA_sign(g_sign,
                              10000, #len genome
                              "96-SBS", #model
                              0.001, #mut_rate_avg
                              0.0, #rate indel
                              100, #size_indel
                              seed,
                              set_mut,
                              0.5,#lavalette_par
                              ["SBS6","SBS22"], #used_sign
                              [0.0], #vector_change_points
                              [0.5 0.5], #vector_activities
                              0.8) #ratio_background_signature

if Sys.iswindows()
    mkpath(".\\Experiments\\Experiment_Signature\\Fasta1\\") # Create folder
    path_for_fasta = ".\\Experiments\\Experiment_Signature\\Fasta1\\"
else
    mkpath("./Experiments/Experiment_Signature/Fasta1/")
    path_for_fasta = "./Experiments/Experiment_Signature/Fasta1/"
end
for le in 1:length(leafs)
    w = FASTA.Writer(open(path_for_fasta
                              * "exp1"
                              * "_Sign"
                              * "_sample"
                              * string(leafs[le])
                              * ".fasta",
                              "w"))
    rec = FASTA.Record("Sample" * string(leafs[le]), fasta_samples[le])
    write(w, rec)
    close(w)
end

call_ART("HS25",
        path_for_fasta,
        150,
        100,
        "exp_1",
        true,
        seed,
        sam = false,
        ef = true,
        mate_pair = false,
        mean_fragsize = 200,
        std_fragsize = 10,
        no_ALN = false)
cd("..")
#####   EXP2
seed = MersenneTwister(1234)
g_sign = copy(tree_red)
g_seq, fasta_samples, position_used, mutations_tot2 = J_Space.experiment_noISA_sign(g_sign,
                                  10000, #len genome
                                  "96-SBS", #model
                                  0.001, #mut_rate_avg
                                  0.0, #rate indel
                                  100, #size_indel
                                  seed,
                                  set_mut,
                                  0.5,#lavalette_par
                                  ["SBS6","SBS22"], #used_sign
                                  [0.0, 100.0], #vector_change_points
                                  [0.0 1.0; 1.0 0.0], #vector_activities
                                  0.8, #ratio_background_signature
                                  frequency_dna=[0.3,0.2,0.2])
mutations_tot
println("save_fasta")
if Sys.iswindows()
    mkpath(".\\Experiments\\Experiment_Signature\\Fasta2\\") # Create folder
    path_for_fasta = ".\\Experiments\\Experiment_Signature\\Fasta2\\"
else
    mkpath("./Experiments/Experiment_Signature/Fasta2/")
    path_for_fasta = "./Experiments/Experiment_Signature/Fasta2/"
end
for le in 1:length(leafs)
    w = FASTA.Writer(open(path_for_fasta
                                  * "exp2_"
                                  * "_Sign"
                                  * "_sample"
                                  * string(leafs[le])
                                  * ".fasta",
                                  "w"))
    rec = FASTA.Record("Sample" * string(leafs[le]), fasta_samples[le])
    write(w, rec)
    close(w)
end
call_ART("HS25",
        path_for_fasta,
        150,
        100,
        "exp_1",
        true,
        seed,
        sam = false,
        ef = true,
        mate_pair = false,
        mean_fragsize = 200,
        std_fragsize = 10,
        no_ALN = false)
cd("..")
########### EXP 3
seed = MersenneTwister(1234)
g_sign = copy(tree_red)
g_seq, fasta_samples, position_used, mutations_tot3 = J_Space.experiment_noISA_sign(g_sign,
                                  10000, #len genome
                                  "96-SBS", #model
                                  0.001, #mut_rate_avg
                                  0.0, #rate indel
                                  100, #size_indel
                                  seed,
                                  set_mut,
                                  0.5,#lavalette_par
                                  ["SBS6","SBS22"], #used_sign
                                  [0.0, 100.0], #vector_change_points
                                  [1.0 0.0; 0.0 1.0], #vector_activities
                                  0.8, #ratio_background_signature
                                  frequency_dna=[0.3,0.2,0.2])
mutations_tot

println("save_fasta")
if Sys.iswindows()
    mkpath(".\\Experiments\\Experiment_Signature\\Fasta3\\") # Create folder
    path_for_fasta = ".\\Experiments\\Experiment_Signature\\Fasta3\\"
else
    mkpath("./Experiments/Experiment_Signature/Fasta3/")
    path_for_fasta = "./Experiments/Experiment_Signature/Fasta3/"
end
for le in 1:length(leafs)
    w = FASTA.Writer(open(path_for_fasta
                          * "exp3_"
                          * "_Sign"
                          * "_sample"
                          * string(leafs[le])
                          * ".fasta",
                          "w"))
    rec = FASTA.Record("Sample" * string(leafs[le]), fasta_samples[le])
    write(w, rec)
    close(w)
end
call_ART("HS25",
        path_for_fasta,
        150,
        100,
        "exp_1",
        true,
        seed,
        sam = false,
        ef = true,
        mate_pair = false,
        mean_fragsize = 200,
        std_fragsize = 10,
        no_ALN = false)

cd("..")
#Save
if Sys.iswindows()
    CSV.write(".\\Experiments\\Experiment_Signature\\Fileoutput\\Mutation_tot_exp1.csv",
                              mutations_tot,
                              delim = ",")

    CSV.write(".\\Experiments\\Experiment_Signature\\Fileoutput\\Mutation_tot_exp2.csv",
                              mutations_tot2,
                              delim = ",")

    CSV.write(".\\Experiments\\Experiment_Signature\\Fileoutput\\Mutation_tot_exp3.csv",
                              mutations_tot3,
                              delim = ",")
else
    CSV.write("./Experiments/Experiment_Signature/Fileoutput/Mutation_tot_exp1.csv",
                              mutations_tot,
                              delim = ",")

    CSV.write("./Experiments/Experiment_Signature/Fileoutput/Mutation_tot_exp2.csv",
                              mutations_tot2,
                              delim = ",")

    CSV.write("./Experiments/Experiment_Signature/Fileoutput/Mutation_tot_exp3.csv",
                              mutations_tot3,
                              delim = ",")
end
