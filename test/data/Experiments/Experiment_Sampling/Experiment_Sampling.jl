using J_Space
using Random
using CSV, PhyloNetworks, Graphs, MetaGraphs, Tables, FASTX, DataFrames
using Distributed
function get_leafs(tree::AbstractGraph)
    leafs = []
    for v in vertices(tree)
        if indegree(tree, v) == 1 && outdegree(tree, v) == 0
            push!(leafs, v)
        end
    end
    return leafs
end

println("load library")

seed = MersenneTwister(1234)

if Sys.iswindows()
    ref = ".\\reference.fasta"
else
    ref = "./reference.fasta"
end

# only for contact mode
params = IdDict( "alpha" => 0.5)

g_meta = spatial_graph(1000, 1000, seed, dim = 3, n_cell=1)
println("simulation...")
#df, G, n_cell_alive, set_mut, Gs_conf, CA_subpop, α_subpop, times =
                                        #J_Space.simulate_evolution_2(g_meta,
df, G, n_cell_alive, set_mut, Gs_conf, CA_subpop, α_subpop=
                                        J_Space.simulate_evolution(g_meta,
                                                           200.0,
                                                           0.4,
                                                           0.01,
                                                           0.0,
                                                           0.0,
                                                           0.3,
                                                           0.1,
                                                           "contact",
                                                           seed)
sampling_cell = [10, 100, 1000]#,5000]


#@distributed for sampling in sampling_cell
@distributed for sampling in sampling_cell
    Times_tree_phy = []
    Memory_tree_phy = []
    Tree_fil_final =[]
    Newick_final = []
    Times_ISA_tot = []
    Times_NoIsa_tot = []
    memory_ISA = []
    memory_NoISA = []

    for i in 1:50
        println("runned sampling: n_cell -> ", sampling ," iteration: ", i)
        println("sampling...")
        value, time, allocated, meta = @timed sampling_phylogentic_relation(G,
                                                          "Random",
                                                          df,
                                                          sampling,
                                                          set_mut,
                                                          seed,
                                                          0)
        matrix_R = value[1]
        #value, times, allocated2, meta = @timed J_Space.create_tree_2(matrix_R, true)
        value, time2, allocated2, meta = @timed create_tree(matrix_R, true, 200.0)
        #println("value[1]", value[1])
        push!(Tree_fil_final, copy(value[1]))
        push!(Newick_final, value[2])
        #time_f = time+value[3]
        #push!(Times_tree_phy, time_f)
        #mem = value[4]
        #mem_phy = allocated+(allocated2 - mem)
        #push!(Memory_tree_phy, mem_phy)
        tree_red = value[1]
        g_isa = copy(tree_red)
        println("ISA...")
        mutation_driver = Dict{}()
        value, time, allocated, meta = @timed experiment_ISA(g_isa,
                                                            0.000001,
                                                            seed,
                                                            ref,
                                                            set_mut)
        fastaX = value[2]
        push!(Times_ISA_tot, time)
        push!(memory_ISA, allocated)

        println("save_fasta")
        leafs = get_leafs(tree_red)
        if Sys.iswindows()
            mkpath(".\\Experiments\\Experiment_Sampling\\Fasta\\") # Create folder
            path_for_fasta = ".\\Experiments\\Experiment_Sampling\\Fasta\\"
        else
            mkpath("./Experiments/Experiment_Sampling/Fasta/")
            path_for_fasta = "./Experiments/Experiment_Sampling/Fasta/"
        end
        for le in 1:length(leafs)
            w = FASTA.Writer(open(path_for_fasta
                                      * "$sampling"
                                      * "_"
                                      * string(i)
                                      * "_ISA"
                                      * "_sample"
                                      * string(leafs[le])
                                      * ".fasta",
                                      "w"))
            rec = FASTA.Record("Sample" * string(leafs[le]), fastaX[le])
            write(w, rec)
            close(w)
        end

        println("NoISA")
        mutation_driver = Dict{}()
        value, time, allocated, meta =
                                    @timed experiment_noISA(tree_red,
                                                            ref,
                                                            "JC69",
                                                            params,
                                                            0.000001,
                                                            100,
                                                            seed,
                                                            set_mut,
                                                            0.5,
                                                            1)

        fastaX = value[2]
        push!(Times_NoIsa_tot, time)
        push!(memory_NoISA, allocated)
        println("Save fasta noISA")

        for le in 1:length(leafs)
            w = FASTA.Writer(open(path_for_fasta
                                      * "$sampling"
                                      * "_"
                                      * string(i)
                                      * "_noISA_indel"
                                      * "_sample"
                                      * string(leafs[le])
                                      * ".fasta",
                                      "w"))
            rec = FASTA.Record("Sample"
                                   * string(leafs[le]),
                                   fastaX[le])
            write(w, rec)
            close(w)
        end
    end

    println("SAVE final...")
    for l in 1:length(Newick_final)
        if Sys.iswindows()
            mkpath(".\\Experiments\\Experiment_Sampling\\Plot")
            mkpath(".\\Experiments\\Experiment_Sampling\\Fileoutput")
            savegraph(".\\Experiments\\Experiment_Sampling\\Plot" *
                      "\\Tree_Fil_$sampling-$l.mg",
                      Tree_fil_final[l])
            writeTopology(Newick_final[l],
                      ".\\Experiments\\Experiment_Sampling\\Fileoutput" *
                      "\\formatNewick_$sampling-$l")
            CSV.write(".\\Experiments\\Experiment_Sampling\\Fileoutput" *
                      "\\time_tree_phy_$sampling.csv",
                      Tables.table(Times_tree_phy),
                      header=false)
            CSV.write(".\\Experiments\\Experiment_Sampling\\Fileoutput" *
                  "\\Memory_tree_phy_$sampling.csv",
                  Tables.table(Memory_tree_phy),
                  header=false)
            CSV.write(".\\Experiments\\Experiment_Sampling\\Fileoutput" *
                  "\\Times_ISA_$sampling.csv",
                  Tables.table(Times_ISA_tot),
                  header=false)
            CSV.write(".\\Experiments\\Experiment_Sampling\\Fileoutput" *
                  "\\Times_NoISA_Indel_$sampling.csv",
                  Tables.table(Times_NoIsa_tot),
                  header=false)
            CSV.write(".\\Experiments\\Experiment_Sampling\\Fileoutput" *
                  "\\memory_ISA_$sampling.csv",
                  Tables.table(memory_ISA),
                  header=false)
            CSV.write(".\\Experiments\\Experiment_Sampling\\Fileoutput" *
                  "\\memory_NoISA_Indel_$sampling.csv",
                  Tables.table(memory_NoISA),
                  header=false)
        else
            mkpath("./Experiments/Experiment_Sampling/Plot")
            mkpath("./Experiments/Experiment_Sampling/Fileoutput")
            savegraph("./Experiments/Experiment_Sampling/Plot" *
                  "/Tree_Fil_$sampling-$l-contact.mg",
                  Tree_fil_final[l])
            writeTopology(Newick_final[l],
                      "./Experiments/Experiment_Sampling/Fileoutput" *
                      "/formatNewick_$sampling-$l-contact")
            CSV.write("./Experiments/Experiment_Sampling/Fileoutput" *
                      "/time_tree_phy_$sampling.csv",
                      Tables.table(Times_tree_phy),
                      header=false)
            CSV.write("./Experiments/Experiment_Sampling/Fileoutput" *
                  "/Memory_tree_phy_$sampling.csv",
                  Tables.table(Memory_tree_phy),
                  header=false)
            CSV.write("./Experiments/Experiment_Sampling/Fileoutput" *
                      "/Times_ISA_$sampling.csv",
                      Tables.table(Times_ISA_tot),
                      header=false)
            CSV.write("./Experiments/Experiment_Sampling/Fileoutput" *
                      "/Times_NoISA_Indel_$sampling.csv",
                      Tables.table(Times_NoIsa_tot),
                      header=false)
            CSV.write("./Experiments/Experiment_Sampling/Fileoutput" *
                      "/memory_ISA_$sampling.csv",
                      Tables.table(memory_ISA),
                      header=false)
            CSV.write("./Experiments/Experiment_Sampling/Fileoutput" *
                      "/memory_NoISA_Indel_$sampling.csv",
                      Tables.table(memory_NoISA),
                      header=false)
        end
    end
end
