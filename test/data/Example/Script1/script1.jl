using J_Space, Random, CSV, PhyloNetworks, Graphs, MetaGraphs, Tables, DataFrames
println("loading libraries")

if Sys.iswindows()
    mkpath(".\\Example\\script2\\Fileoutput\\") # Create folder
    mkpath(".\\Example\\script2\\Plot\\") # Create folder
else
    mkpath("./Example/script2/Fileoutput/") # Create folder
    mkpath("./Example/script2/Plot/") # Create folder
end

seed = MersenneTwister(1234)

birth_rate_min = 0.1
birth_rate_max = 0.6
n_rep = 10


values_of_birth_rate = range(birth_rate_min,birth_rate_max,length=n_rep)
G_state_final = []
dinamica_final = []
n_cell_alive_final = []
Set_mut_final = []
CA_subpop_final = []
α_subpop_final = []
Tree_mut_final = []
Tree_fil_final =[]
Newick_final = []

println("generate graph...")
g_meta = spatial_graph(200, 200, seed, dim = 3, n_cell=1)

#@distributed for birth_rate in values_of_birth_rate
for birth_rate in values_of_birth_rate
    println("runned configuration with birth_rate -> ", birth_rate ,
                                                       " in 3D mode -> contact")
    println("simulation...")
    df, G, n_cell_alive, set_mut, Gs_conf, CA_subpop, α_subpop =
                                        simulate_evolution(copy(g_meta),
                                                                   200.0,
                                                                   birth_rate,
                                                                   0.01,
                                                                   0.01,
                                                                   0.00001,
                                                                   0.2,
                                                                   0.1,
                                                                   "contact",
                                                                   seed)

    if n_cell_alive[end] < 100
         println("Simulation with less than 1000 cells: simulation discarded")
         continue
    end
    push!(G_state_final, copy(G))
    push!(dinamica_final, copy(df))
    push!(n_cell_alive_final, copy(n_cell_alive))
    push!(Set_mut_final, copy(set_mut))
    push!(CA_subpop_final, copy(CA_subpop))
    push!(α_subpop_final, copy(α_subpop))

    println("Sampling...")
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
println("SAVE DATA")
for l in 1:length(G_state_final)
    if Sys.iswindows()
        CSV.write(".\\Example\\script1\\Fileoutput" *
                  "\\Dynamics-$l.csv",
                  dinamica_final[l],
                  delim = ",")
        CSV.write(".\\Example\\script1\\Fileoutput" *
                  "\\n_cell_alive_final-$l.csv",
                  Tables.table(n_cell_alive_final[l]),
                  header=false)
        savegraph(".\\Example\\script1\\Plot" *
                  "\\Final_conf-$l.mg",G_state_final[l])
        List_driver = DataFrame(Driver = Set_mut_final[l],
                                Fitness = α_subpop_final[l])
        CSV.write(".\\Example\\script1" *
                  "\\Fileoutput\\DriverList-$l.csv",
                  List_driver,delim = ",")
        CSV.write(".\\Example\\script1\\Fileoutput" *
                  "\\CA_subpop_final-$l.csv",
                  Tables.table(CA_subpop_final[l]),
                  header=false)
    else
        CSV.write("./Example/script1/Fileoutput" *
                  "/Dynamics-$l.csv",
                  dinamica_final[l],
                  delim = ",")
        CSV.write("./Example/script1/Fileoutput" *
                  "/n_cell_alive_final-$l.csv",
                  Tables.table(n_cell_alive_final[l]),
                  header=false)
        savegraph("./Example/script1/Plot" *
                  "/Final_conf-$l.mg",
                  G_state_final[l])
        List_driver = DataFrame(Driver = Set_mut_final[l],
                                Fitness = α_subpop_final[l])
        CSV.write("./Example/script1/Fileoutput" *
                  "/DriverList-$l.csv",
                  List_driver,delim = ",")
        CSV.write("./Example/script1/Fileoutput" *
                  "/CA_subpop_final-$l.csv",
                  Tables.table(CA_subpop_final[l]),
                  header=false)
    end
end
