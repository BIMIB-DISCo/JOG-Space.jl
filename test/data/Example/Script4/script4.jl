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
    CSV.write(".\\Example\\script4\\Fileoutput" *
              "\\Dynamics.csv",
              df,
              delim = ",")
    CSV.write(".\\Example\\script4\\Fileoutput" *
              "\\n_cell_alive_final.csv",
              Tables.table(n_cell_alive),
              header=false)
    savegraph(".\\Example\\script4\\Plot" *
              "\\Final_conf.mg",G)
    List_driver = DataFrame(Driver = set_mut,
                            Fitness = α_subpop)
    CSV.write(".\\Example\\script4" *
              "\\Fileoutput\\DriverList.csv",
              List_driver,delim = ",")
    CSV.write(".\\Example\\script4\\Fileoutput" *
              "\\CA_subpop_final.csv",
              Tables.table(CA_subpop),
              header=false)
    writeTopology(net,".\\Example\\script4\\Fileoutput\\formatNewick")
else
    CSV.write("./Example/script4/Fileoutput" *
              "/Dynamics.csv",
              df,
              delim = ",")
    CSV.write("./Example/script4/Fileoutput" *
              "/n_cell_alive_final.csv",
              Tables.table(n_cell_alive),
              header=false)
    savegraph("./Example/script4/Plot" *
              "/Final_conf.mg",
              G)
    List_driver = DataFrame(Driver = set_mut,
                            Fitness = α_subpop)
    CSV.write("./Example/script4/Fileoutput" *
              "/DriverList.csv",
              List_driver,delim = ",")
    CSV.write("./Example/script4/Fileoutput" *
              "/CA_subpop_final.csv",
              Tables.table(CA_subpop),
              header=false)
    writeTopology(net,"./Example/script4/Fileoutput/formatNewick")
end

path_for_fasta = ""
if Sys.iswindows()
    mkpath(".\\Example\\script4\\Fasta\\") # Create folder
    path_for_fasta = ".\\Example\\script4\\Fasta\\"
else
    mkpath("./Example/script4/Fasta/")
    path_for_fasta = "./Example/script4/Fasta/"
end
leafs = J_Space.get_leafs(tree_red)

g_sign = copy(tree_red)
g_seq, fasta_samples, Tree, mutations_tot = experiment_noISA_sign(
                              g_sign,
                              10000, #len genome
                              "SBS-37", #model
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
println("save fasta")
for le in 1:length(leafs)
    w = FASTA.Writer(open(path_for_fasta
                              * "_sample"
                              * string(leafs[le])
                              * ".fasta",
                              "w"))
    rec = FASTA.Record("Sample" * string(leafs[le]), fasta_samples[le])
    write(w, rec)
    close(w)
end

if Sys.iswindows()
    path_for_fasta = ".\\Example\\script4\\Fasta\\"
    mkpath(".\\Example\\script4\\Fasta\\$coverage\\") # Create folder
    path_for_output = ".\\Example\\script4\\Fasta\\$coverage\\"
else
    path_for_fasta = "./Example/script4/Fasta//"
    mkpath("./Example/script4/Fasta/$coverage/") # Create folder
    path_for_output = "./Example/script4/Fasta/$coverage/"
end


coverage_max = 100
coverage_min = 10
n_rep = 10
#ATTENCTION, Coverage is INTEGER
values_of_covarage = range(coverage_min,coverage_max,length=n_rep)
for coverage in values_of_covarage
    println("runned configuration with coverage -> ", coverage)
    coverage = convert(Int64, coverage)

    if Sys.iswindows()
        path_for_fasta = pwd()*"\\Example\\script4\\Fasta\\"
        mkpath(".\\Example\\script4\\Fasta\\coverage_$coverage\\") # Create folder
        path_for_output = ".\\Example\\script4\\Fasta\\coverage_$coverage\\"
    else
        path_for_fasta = pwd()*"/Example/script4/Fasta/"
        mkpath("./Example/script4/Fasta/coverage_$coverage/") # Create folder
        path_for_output = "./Example/script4/Fasta/coverage_$coverage/"
    end

    call_ART("HS25",
            path_for_fasta,
            path_for_output,
            150,
            coverage,
            "exp_1",
            true,
            seed,
            sam = false,
            ef = true,
            mate_pair = false,
            mean_fragsize = 200,
            std_fragsize = 10,
            no_ALN = false)
end
