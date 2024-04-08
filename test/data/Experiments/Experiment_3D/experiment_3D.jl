using J_Space
using Random
using CSV, PhyloNetworks, Graphs, MetaGraphs, Tables, FASTX, DataFrames

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

mut_driver_rate = [0.0, 0.00000001, 0.000001, 0.0001]

# only for contact mode
params = IdDict( "alpha" => 0.5)
T_isa = []
T_noIsa = []

for mut_rate in mut_driver_rate
    println("runned configuration: mut_driver_rate -> ", mut_rate ,
                    " in 3D mode -> contact")
    dinamica_final = []
    n_cell_alive_final = []
    Set_mut_final = []
    CA_subpop_final = []
    α_subpop_final = []
    G_state_final = []
    Tree_mut_final = []
    Tree_fil_final =[]
    Newick_final = []

    for i in 1:10
        println(i)
        g_meta = spatial_graph(71, 71, seed, dim = 3, n_cell=1)
        println("simulation...")
        df, G, n_cell_alive, set_mut, Gs_conf, CA_subpop, α_subpop =
                                                simulate_evolution(g_meta,
                                                                   200.0,
                                                                   0.4,
                                                                   0.01,
                                                                   0.0,
                                                                   mut_rate,
                                                                   0.2,
                                                                   0.1,
                                                                   "contact",
                                                                   seed)
        if n_cell_alive[end] < 1000
            println("Simulation with less than 1000 cells: simulation ",
                    "discarded")
            continue
        end
        push!(G_state_final, copy(G))
        push!(dinamica_final, df)
        push!(n_cell_alive_final, n_cell_alive)
        push!(Set_mut_final, set_mut)
        push!(CA_subpop_final, CA_subpop)
        push!(α_subpop_final, α_subpop)

        println("sampling...")
        matrix_R, tree_mut = sampling_phylogentic_relation(G,
                                                           "Random",
                                                           df,
                                                           100,
                                                           set_mut,
                                                           seed,
                                                           1)
        push!(Tree_mut_final, copy(tree_mut))
        tree_red, net = create_tree(matrix_R, true, 200.0)
        push!(Tree_fil_final, copy(tree_red))
        push!(Newick_final, net)

        g_isa = copy(tree_red)
        println("ISA...")
        mutation_driver = Dict{}()
        g_seq, fastaX, position_used, mutation_driver =
                                            experiment_ISA(g_isa,
                                                           0.00000001,
                                                           seed,
                                                           ref,
                                                           set_mut)
        if isempty(mutation_driver) == false
            if Sys.iswindows()
                mkpath(".\\Experiments\\Experiment_3D\\Fileoutput\\")
                CSV.write(".\\Experiments\\Experiment_3D\\Fileoutput" *
                          "\\Mutation_driver_ISA_$mut_rate-$i-contact.csv",
                          mutation_driver)
            else
                mkpath("./Experiments/Experiment_3D/Fileoutput/")
                CSV.write("./Experiments/Experiment_3D/Fileoutput" *
                          "/Mutation_driver_ISA_$mut_rate-$i-contact.csv",
                          mutation_driver)
            end
        end
        println("save_fasta")
        leafs = get_leafs(tree_red)
        if Sys.iswindows()
            mkpath(".\\Experiments\\Experiment_3D\\Fasta\\") # Create folder
            path_for_fasta = ".\\Experiments\\Experiment_3D\\Fasta\\"
        else
            mkpath("./Experiments/Experiment_3D/Fasta/")
            path_for_fasta = "./Experiments/Experiment_3D/Fasta/"
        end
        for le in 1:length(leafs)
            w = FASTA.Writer(open(path_for_fasta
                                          * "$mut_rate"
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
        g_seq, fastaX, Tree_SC, mutation_driver =
                                        experiment_noISA(tree_red,
                                                        ref,
                                                        "JC69",
                                                        params,
                                                        0.00000001,
                                                        100,
                                                        seed,
                                                        set_mut,
                                                        0.5,
                                                        1)
        if isempty(mutation_driver) == false
            if Sys.iswindows()
                CSV.write(".\\Experiments\\Experiment_3D\\Fileoutput" *
                          "\\Mutation_driver_noISA_$mut_rate-$i-contact.csv",
                          mutation_driver)
            else
                CSV.write("./Experiments/Experiment_3D/Fileoutput" *
                          "/Mutation_driver_noISA_$mut_rate-$i-contact.csv",
                          mutation_driver)
            end
        end

        println("Save fasta noISA")

        for le in 1:length(leafs)
            w = FASTA.Writer(open(path_for_fasta
                                          * "$mut_rate"
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

    end
    println("SAVE final...")
    for l in 1:length(Newick_final)
        if Sys.iswindows()
            mkpath(".\\Experiments\\Experiment_3D\\Plot\\")
            CSV.write(".\\Experiments\\Experiment_3D\\Fileoutput" *
                      "\\Dinamica_$mut_rate-$l-contact.csv",
                      dinamica_final[l],
                      delim = ",")
            CSV.write(".\\Experiments\\Experiment_3D\\Fileoutput" *
                      "\\n_cell_alive_final_$mut_rate-$l-contact.csv",
                      Tables.table(n_cell_alive_final[l]),
                      header=false)
            savegraph(".\\Experiments\\Experiment_3D\\Plot" *
                      "\\Final_conf_$mut_rate-$l-contact.lgz",
                       G_state_final[l])
            List_driver = DataFrame(Driver = Set_mut_final[l],
                                    Fitness = α_subpop_final[l])
            CSV.write(".\\Experiments\\Experiment_3D" *
                      "\\Fileoutput\\DriverList_$mut_rate-$l-contact.csv",
                      List_driver,delim = ",")
            CSV.write(".\\Experiments\\Experiment_3D\\Fileoutput" *
                      "\\CA_subpop_final_$mut_rate-$l-contact.csv",
                      Tables.table(CA_subpop_final[l]),
                      header=false)
            savegraph(".\\Experiments\\Experiment_3D\\Plot" *
                      "\\Tree_Muts_$mut_rate-$l-contact.lgz",
                      Tree_mut_final[l])
            savegraph(".\\Experiments\\Experiment_3D\\Plot" *
                      "\\Tree_Fil_$mut_rate-$l-contact.lgz",
                      Tree_fil_final[l])
            writeTopology(Newick_final[l],
                          ".\\Experiments\\Experiment_3D\\Fileoutput" *
                          "\\formatNewick_$mut_rate-$l-contact")
        else
            mkpath("./Experiments/Experiment_3D/Plot/")
            CSV.write("./Experiments/Experiment_3D/Fileoutput" *
                      "/Dinamica_$mut_rate-$l-contact.csv",
                      dinamica_final[l],
                      delim = ",")
            CSV.write("./Experiments/Experiment_3D/Fileoutput" *
                      "/n_cell_alive_final_$mut_rate-$l-contact.csv",
                      Tables.table(n_cell_alive_final[l]),
                      header=false)
            savegraph("./Experiments/Experiment_3D/Plot" *
                      "/Final_conf_$mut_rate-$l-contact.lgz",
                      G_state_final[l])
            List_driver = DataFrame(Driver = Set_mut_final[l],
                                    Fitness = α_subpop_final[l])
            CSV.write("./Experiments/Experiment_3D/Fileoutput" *
                      "/DriverList_$mut_rate-$l-contact.csv",
                      List_driver,delim = ",")
            CSV.write("./Experiments/Experiment_3D/Fileoutput" *
                      "/CA_subpop_final_$mut_rate-$l-contact.csv",
                      Tables.table(CA_subpop_final[l]),
                      header=false)
            savegraph("./Experiments/Experiment_3D/Plot" *
                      "/Tree_Muts_$mut_rate-$l-contact.lgz",
                      Tree_mut_final[l])
            savegraph("./Experiments/Experiment_3D/Plot" *
                      "/Tree_Fil_$mut_rate-$l-contact.lgz",
                      Tree_fil_final[l])
            writeTopology(Newick_final[l],
                          "./Experiments/Experiment_3D/Fileoutput" *
                          "/formatNewick_$mut_rate-$l-contact")
        end
    end
end




for mut_rate in mut_driver_rate
    println("runned configuration: mut_driver_rate -> ", mut_rate ,
                    " in 3D mode -> voter")

    dinamica_final = []
    n_cell_alive_final = []
    Set_mut_final = []
    CA_subpop_final = []
    α_subpop_final = []
    G_state_final = []
    Tree_mut_final = []
    Tree_fil_final =[]
    Newick_final = []

    for i in 1:10
        println(i)
        g_meta = spatial_graph(71, 71, seed, dim = 3, n_cell=1)
        println("simulation...")
        df, G, n_cell_alive, set_mut, Gs_conf, CA_subpop, α_subpop =
                                                    simulate_evolution(g_meta,
                                                                       200.0,
                                                                       0.4,
                                                                       0.01,
                                                                       0.0,
                                                                       mut_rate,
                                                                       0.2,
                                                                       0.1,
                                                                       "voter",
                                                                       seed)
        if n_cell_alive[end] < 1000
            println("Simulation with less than 1000 cells: simulation ",
                    "discarded")
            continue
        end
        push!(G_state_final, copy(G))
        push!(dinamica_final, df)
        push!(n_cell_alive_final, n_cell_alive)
        push!(Set_mut_final, set_mut)
        push!(CA_subpop_final, CA_subpop)
        push!(α_subpop_final, α_subpop)

        println("sampling...")
        matrix_R, tree_mut = sampling_phylogentic_relation(G,
                                                           "Random",
                                                           df,
                                                           100,
                                                           set_mut,
                                                           seed,
                                                           1)

        push!(Tree_mut_final, copy(tree_mut))
        tree_red, net = create_tree(matrix_R, true, 200.0)
        push!(Tree_fil_final, copy(tree_red))
        push!(Newick_final, net)

    end

    println("SAVE final...")
    for l in 1:length(Newick_final)
        if Sys.iswindows()
            CSV.write(".\\Experiments\\Experiment_3D\\Fileoutput" *
                      "\\Dinamica_$mut_rate-$l-voter.csv",
                      dinamica_final[l],
                      delim = ",")
            CSV.write(".\\Experiments\\Experiment_3D\\Fileoutput" *
                      "\\n_cell_alive_final_$mut_rate-$l-voter.csv",
                      Tables.table(n_cell_alive_final[l]),
                      header=false)
            savegraph(".\\Experiments\\Experiment_3D\\Plot" *
                      "\\Final_conf_$mut_rate-$l-voter.lgz",
                        G_state_final[l])
            List_driver = DataFrame(Driver = Set_mut_final[l],
                                    Fitness = α_subpop_final[l])
            CSV.write(".\\Experiments\\Experiment_3D" *
                      "\\Fileoutput\\DriverList_$mut_rate-$l-voter.csv",
                      List_driver,delim = ",")
            CSV.write(".\\Experiments\\Experiment_3D\\Fileoutput" *
                      "\\CA_subpop_final_$mut_rate-$l-voter.csv",
                      Tables.table(CA_subpop_final[l]),
                      header=false)
            savegraph(".\\Experiments\\Experiment_3D\\Plot" *
                      "\\Tree_Muts_$mut_rate-$l-voter.lgz",
                      Tree_mut_final[l])
            savegraph(".\\Experiments\\Experiment_3D\\Plot" *
                      "\\Tree_Fil_$mut_rate-$l-voter.lgz",
                      Tree_fil_final[l])
            writeTopology(Newick_final[l],
                          ".\\Experiments\\Experiment_3D\\Fileoutput" *
                          "\\formatNewick_$mut_rate-$l-voter")
        else
            CSV.write("./Experiments/Experiment_3D/Fileoutput" *
                      "/Dinamica_$mut_rate-$l-voter.csv",
                      dinamica_final[l],
                      delim = ",")
            CSV.write("./Experiments/Experiment_3D/Fileoutput" *
                      "/n_cell_alive_final_$mut_rate-$l-voter.csv",
                      Tables.table(n_cell_alive_final[l]),
                      header=false)
            savegraph("./Experiments/Experiment_3D/Plot" *
                      "/Final_conf_$mut_rate-$l-voter.lgz",
                      G_state_final[l])
            List_driver = DataFrame(Driver = Set_mut_final[l],
                                    Fitness = α_subpop_final[l])
            CSV.write("./Experiments/Experiment_3D/Fileoutput" *
                      "/DriverList_$mut_rate-$l-voter.csv",
                      List_driver,delim = ",")
            CSV.write("./Experiments/Experiment_3D/Fileoutput" *
                      "/CA_subpop_final_$mut_rate-$l-voter.csv",
                      Tables.table(CA_subpop_final[l]),
                      header=false)
            savegraph("./Experiments/Experiment_3D/Plot" *
                      "/Tree_Muts_$mut_rate-$l-voter.lgz",
                      Tree_mut_final[l])
            savegraph("./Experiments/Experiment_3D/Plot" *
                      "/Tree_Fil_$mut_rate-$l-voter.lgz",
                      Tree_fil_final[l])
            writeTopology(Newick_final[l],
                          "./Experiments/Experiment_3D/Fileoutput" *
                          "/formatNewick_$mut_rate-$l-voter")
        end
    end
end

for mut_rate in mut_driver_rate
    println("runned configuration: mut_driver_rate -> ", mut_rate ,
                    " in 3D mode -> h_voter")
    dinamica_final = []
    n_cell_alive_final = []
    Set_mut_final = []
    CA_subpop_final = []
    α_subpop_final = []
    G_state_final = []
    Tree_mut_final = []
    Tree_fil_final =[]
    Newick_final = []

    for i in 1:10
        println(i)
        g_meta = spatial_graph(71, 71, seed, dim = 3, n_cell=1)
        println("simulation...")
        df, G, n_cell_alive, set_mut, Gs_conf, CA_subpop, α_subpop =
                                                   simulate_evolution(g_meta,
                                                                      200.0,
                                                                      0.4,
                                                                      0.01,
                                                                      0.0,
                                                                      mut_rate,
                                                                      0.2,
                                                                      0.1,
                                                                      "h_voter",
                                                                      seed)
        if n_cell_alive[end] < 1000
            println("Simulation with less than 1000 cells: simulation ",
                    "discarded")
            continue
        end
        push!(G_state_final, copy(G))
        push!(dinamica_final, df)
        push!(n_cell_alive_final, n_cell_alive)
        push!(Set_mut_final, set_mut)
        push!(CA_subpop_final, CA_subpop)
        push!(α_subpop_final, α_subpop)

        println("sampling...")
        matrix_R, tree_mut = sampling_phylogentic_relation(G,
                                                              "Random",
                                                              df,
                                                              100,
                                                              set_mut,
                                                              seed,
                                                              1)
        push!(Tree_mut_final, copy(tree_mut))
        tree_red, net = create_tree(matrix_R, true, 200.0)
        push!(Tree_fil_final, copy(tree_red))
        push!(Newick_final, net)

    end
    println("SAVE final...")
    for l in 1:length(Newick_final)
        if Sys.iswindows()
            CSV.write(".\\Experiments\\Experiment_3D\\Fileoutput" *
                      "\\Dinamica_$mut_rate-$l-hvoter.csv",
                      dinamica_final[l],
                      delim = ",")
            CSV.write(".\\Experiments\\Experiment_3D\\Fileoutput" *
                      "\\n_cell_alive_final_$mut_rate-$l-hvoter.csv",
                      Tables.table(n_cell_alive_final[l]),
                      header=false)
            savegraph(".\\Experiments\\Experiment_3D\\Plot" *
                      "\\Final_conf_$mut_rate-$l-hvoter.lgz",
                        G_state_final[l])
            List_driver = DataFrame(Driver = Set_mut_final[l],
                                    Fitness = α_subpop_final[l])
            CSV.write(".\\Experiments\\Experiment_3D" *
                      "\\Fileoutput\\DriverList_$mut_rate-$l-hvoter.csv",
                      List_driver,delim = ",")
            CSV.write(".\\Experiments\\Experiment_3D\\Fileoutput" *
                      "\\CA_subpop_final_$mut_rate-$l-hvoter.csv",
                      Tables.table(CA_subpop_final[l]),
                      header=false)
            savegraph(".\\Experiments\\Experiment_3D\\Plot" *
                      "\\Tree_Muts_$mut_rate-$l-hvoter.lgz",
                      Tree_mut_final[l])
            savegraph(".\\Experiments\\Experiment_3D\\Plot" *
                      "\\Tree_Fil_$mut_rate-$l-hvoter.lgz",
                      Tree_fil_final[l])
            writeTopology(Newick_final[l],
                          ".\\Experiments\\Experiment_3D\\Fileoutput" *
                          "\\formatNewick_$mut_rate-$l-hvoter")
        else
            CSV.write("./Experiments/Experiment_3D/Fileoutput" *
                      "/Dinamica_$mut_rate-$l-hvoter.csv",
                      dinamica_final[l],
                      delim = ",")
            CSV.write("./Experiments/Experiment_3D/Fileoutput" *
                      "/n_cell_alive_final_$mut_rate-$l-hvoter.csv",
                      Tables.table(n_cell_alive_final[l]),
                      header=false)
            savegraph("./Experiments/Experiment_3D/Plot" *
                      "/Final_conf_$mut_rate-$l-hvoter.lgz",
                      G_state_final[l])
            List_driver = DataFrame(Driver = Set_mut_final[l],
                                    Fitness = α_subpop_final[l])
            CSV.write("./Experiments/Experiment_3D/Fileoutput" *
                      "/DriverList_$mut_rate-$l-hvoter.csv",
                      List_driver,delim = ",")
            CSV.write("./Experiments/Experiment_3D/Fileoutput" *
                      "/CA_subpop_final_$mut_rate-$l-hvoter.csv",
                      Tables.table(CA_subpop_final[l]),
                      header=false)
            savegraph("./Experiments/Experiment_3D/Plot" *
                      "/Tree_Muts_$mut_rate-$l-hvoter.lgz",
                      Tree_mut_final[l])
            savegraph("./Experiments/Experiment_3D/Plot" *
                      "/Tree_Fil_$mut_rate-$l-hvoter.lgz",
                      Tree_fil_final[l])
            writeTopology(Newick_final[l],
                          "./Experiment_3D/Fileoutput" *
                          "/formatNewick_$mut_rate-$l-hvoter")
        end
    end
end
