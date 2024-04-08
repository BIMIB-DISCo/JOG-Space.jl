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

if Sys.iswindows()
    ref = ".\\reference.fasta"
else
    ref = "./reference.fasta"
end

# only for contact mode
params = IdDict( "alpha" => 0.5)

seed = MersenneTwister(1234)

g_meta = spatial_graph(1000, 1000, seed, dim = 3, n_cell=1)
println("simulation...")
df, G, n_cell_alive, set_mut, Gs_conf, CA_subpop, α_subpop =
                                        simulate_evolution(g_meta,
                                                           200.0,
                                                           0.4,
                                                           0.01,
                                                           0.0,
                                                           0.0,
                                                           0.3,
                                                           0.1,
                                                           "contact",
                                                           seed)

matrix_R, tree_mut = sampling_phylogentic_relation(G,
                                                  "Random",
                                                  df,
                                                  100,
                                                  set_mut,
                                                  seed,
                                                  0)

tree_red, newick = J_Space.create_tree(matrix_R, true, 200.0)
length_genome = [1000, 10000, 100000, 1000000]

leafs = get_leafs(tree_red)

@distributed for len in length_genome
    Times_ISA_tot = []
    Times_NoIsa_tot = []
    Times_NoIsa_indel_tot = []
    memory_ISA = []
    memory_NoISA = []
    memory_NoISA_indel = []
    for i in 1:50
        println("runned: len_genome -> ", len ," iteration: ", i)
        seed_trake = copy(seed)
        g_isa = copy(tree_red)
        println("ISA...")
        value, time, allocated, meta =@timed J_Space.experiment_ISA(g_isa,
                               0.000001, #neutral mutational rate
                               seed,
                               len, #length genome
                               set_mut,
                               frequency_dna = [0.25,0.25,0.25])
        fastaX = value[2]
        push!(Times_ISA_tot, time)
        push!(memory_ISA, allocated)

        println("save_fasta")

        if Sys.iswindows()
            mkpath(".\\Experiments\\Experiment_Genome\\Fasta\\") # Create folder
            path_for_fasta = ".\\Experiments\\Experiment_Genome\\Fasta\\"
        else
            mkpath("./Experiments/Experiment_Genome/Fasta/")
            path_for_fasta = "./Experiments/Experiment_Genome/Fasta/"
        end
        for le in 1:length(leafs)
            w = FASTA.Writer(open(path_for_fasta
                                      * "$len"
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

        seed_noISa = copy(seed_trake)
        g_noISA = copy(tree_red)
        value, time, allocated, meta = @timed J_Space.experiment_noISA(g_noISA,
                                 len, #len genome
                                 "JC69", #model
                                 params,
                                 0.0, #rate indel
                                 100, #max size indel
                                 seed_noISa,
                                 set_mut,
                                 0.5, #lavalette_par
                                 1, #approx_snv_indel
                                 frequency_dna = [0.25,0.25,0.25])
        fastaX = value[2]
        push!(Times_NoIsa_tot, time)
        push!(memory_NoISA, allocated)

        println("Save fasta noISA")

        for le in 1:length(leafs)
            w = FASTA.Writer(open(path_for_fasta
                                      * "$len"
                                      * "_"
                                      * string(i)
                                      * "_noISA"
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

        println("NoISA Indel")
        seed_noISa_indel = copy(seed_trake)
        g_noISA = copy(tree_red)
        value, time, allocated, meta = @timed J_Space.experiment_noISA(g_noISA,
                                 len, #len genome
                                 "JC69", #model
                                 params,
                                 0.000001, #rate indel
                                 100, #max size indel
                                 seed_noISa_indel,
                                 set_mut,
                                 0.5, #lavalette_par
                                 1, #approx_snv_indel
                                 frequency_dna = [0.25,0.25,0.25])

        fastaX = value[2]
        push!(Times_NoIsa_indel_tot, time)
        push!(memory_NoISA_indel, allocated)

        println("Save fasta noISA Indel")

        for le in 1:length(leafs)
            w = FASTA.Writer(open(path_for_fasta
                                      * "$len"
                                      * "_"
                                      * string(i)
                                      * "_noISA"
                                      * "_Indel"
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
    if Sys.iswindows()
        mkpath(".\\Experiments\\Experiment_Genome\\Fileoutput\\")
        CSV.write(".\\Experiments\\Experiment_Genome\\Fileoutput" *
              "\\Times_ISA_$len.csv",
              Tables.table(Times_ISA_tot),
              header=false)
        CSV.write(".\\Experiments\\Experiment_Genome\\Fileoutput" *
              "\\Times_NoISA_$len.csv",
              Tables.table(Times_NoIsa_tot),
              header=false)
        CSV.write(".\\Experiments\\Experiment_Genome\\Fileoutput" *
              "\\Times_NoISA_indel_$len.csv",
              Tables.table(Times_NoIsa_indel_tot),
              header=false)
        CSV.write(".\\Experiments\\Experiment_Genome\\Fileoutput" *
              "\\memory_ISA_$len.csv",
              Tables.table(memory_ISA),
              header=false)
        CSV.write(".\\Experiments\\Experiment_Genome\\Fileoutput" *
              "\\memory_NoISA_$len.csv",
              Tables.table(memory_NoISA),
              header=false)
        CSV.write(".\\Experiments\\Experiment_Genome\\Fileoutput" *
              "\\memory_NoISA_indel_$len.csv",
              Tables.table(memory_NoISA_indel),
              header=false)
    else
        mkpath("./Experiments/Experiment_Genome/Fileoutput")
        CSV.write("./Experiments/Experiment_Genome/Fileoutput" *
                  "/Times_ISA_$len.csv",
                  Tables.table(Times_ISA_tot),
                  header=false)
        CSV.write("./Experiments/Experiment_Genome/Fileoutput" *
                  "/Times_NoISA_$len.csv",
                  Tables.table(Times_NoIsa_tot),
                  header=false)
        CSV.write("./Experiments/Experiment_Genome/Fileoutput" *
                  "/Times_NoISA_indel_$len.csv",
                  Tables.table(Times_NoIsa_indel_tot),
                  header=false)
        CSV.write("./Experiments/Experiment_Genome/Fileoutput" *
                  "/memory_ISA_$len.csv",
                  Tables.table(memory_ISA),
                  header=false)
        CSV.write("./Experiments/Experiment_Genome/Fileoutput" *
                  "/memory_NoISA_$len.csv",
                  Tables.table(memory_NoISA),
                  header=false)
        CSV.write("./Experiments/Experiment_Genome/Fileoutput" *
                  "/memory_NoISA_indel_$len.csv",
                  Tables.table(memory_NoISA_indel),
                  header=false)
    end
end
