using J_Space, Random, CSV, PhyloNetworks, Graphs, MetaGraphs, Tables, DataFrames
using FASTX
println("loading libraries")

if Sys.iswindows()
    mkpath(".\\Example\\script2\\Fileoutput\\") # Create folder
    mkpath(".\\Example\\script2\\Plot\\") # Create folder
else
    mkpath("./Example/script2/Fileoutput/") # Create folder
    mkpath("./Example/script2/Plot/") # Create folder
end

seed = MersenneTwister(1234)

println("generate graph...")
g_meta = spatial_graph(200, 200, seed, dim = 3, n_cell=1)

println("simulation...")
df, G, n_cell_alive, set_mut, Gs_conf, CA_subpop, α_subpop =
                                    simulate_evolution(g_meta,
                                                               200.0,
                                                               0.3,
                                                               0.01,
                                                               0.01,
                                                               0.00001,
                                                               0.2,
                                                               0.1,
                                                               "contact",
                                                               seed)
matrix_R, tree_mut = sampling_phylogentic_relation(G,
                                                   "Random",
                                                    df,
                                                    100,
                                                    set_mut,
                                                    seed,
                                                    1)

tree_red, net = create_tree(matrix_R, true, 200.0)

println("save data dynamics")
if Sys.iswindows()
    CSV.write(".\\Example\\script3\\Fileoutput" *
              "\\Dynamics.csv",
              df,
              delim = ",")
    CSV.write(".\\Example\\script3\\Fileoutput" *
              "\\n_cell_alive_final.csv",
              Tables.table(n_cell_alive),
              header=false)
    savegraph(".\\Example\\script3\\Plot" *
              "\\Final_conf.mg",G)
    List_driver = DataFrame(Driver = set_mut,
                            Fitness = α_subpop)
    CSV.write(".\\Example\\script3" *
              "\\Fileoutput\\DriverList.csv",
              List_driver,delim = ",")
    CSV.write(".\\Example\\script3\\Fileoutput" *
              "\\CA_subpop_final.csv",
              Tables.table(CA_subpop),
              header=false)
    writeTopology(net,".\\Example\\script3\\Fileoutput\\formatNewick")
else
    CSV.write("./Example/script3/Fileoutput" *
              "/Dynamics.csv",
              df,
              delim = ",")
    CSV.write("./Example/script3/Fileoutput" *
              "/n_cell_alive_final.csv",
              Tables.table(n_cell_alive),
              header=false)
    savegraph("./Example/script3/Plot" *
              "/Final_conf.mg",
              G)
    List_driver = DataFrame(Driver = set_mut,
                            Fitness = α_subpop)
    CSV.write("./Example/script3/Fileoutput" *
              "/DriverList.csv",
              List_driver,delim = ",")
    CSV.write("./Example/script3/Fileoutput" *
              "/CA_subpop_final.csv",
              Tables.table(CA_subpop),
              header=false)
    writeTopology(net,"./Example/script3/Fileoutput/formatNewick")
end

path_for_fasta = ""
if Sys.iswindows()
    mkpath(".\\Example\\script3\\Fasta\\") # Create folder
    path_for_fasta = ".\\Example\\script3\\Fasta\\"
else
    mkpath("./Example/script3/Fasta/")
    path_for_fasta = "./Example/script3/Fasta/"
end

leafs = J_Space.get_leafs(tree_red)
mut_rate_avg_max = 0.001
mut_rate_avg_min = 0.00001
n_rep = 10
values_of_mut_rate_avg = range(mut_rate_avg_min,mut_rate_avg_max,length=n_rep)

for mut_rate_avg in values_of_mut_rate_avg
    println("runned configuration with mut_rate_avg -> ", mut_rate_avg)
    seed = MersenneTwister(1234)
    g_sign = copy(tree_red)
    g_seq, fasta_samples, position_used, mutations_tot = experiment_noISA_sign(
                                  g_sign,
                                  10000, #len genome
                                  "SBS-37", #model
                                  mut_rate_avg, #mut_rate_avg
                                  0.0, #rate indel
                                  100, #size_indel
                                  seed,
                                  set_mut,
                                  0.5,#lavalette_par
                                  ["SBS6","SBS22"], #used_sign
                                  [0.0], #vector_change_points
                                  [0.5 0.5], #vector_activities
                                  0.8) #ratio_background_signature
    println("save fasta")
    for le in 1:length(leafs)
        w = FASTA.Writer(open(path_for_fasta
                                  * "exp_$mut_rate_avg"
                                  * "_sample"
                                  * string(leafs[le])
                                  * ".fasta",
                                  "w"))
        rec = FASTA.Record("Sample" * string(leafs[le]), fasta_samples[le])
        write(w, rec)
        close(w)
    end
end
