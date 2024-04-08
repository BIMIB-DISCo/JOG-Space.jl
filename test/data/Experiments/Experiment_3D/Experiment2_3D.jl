using J_Space
using Random
using CSV, Graphs, MetaGraphs, Tables, DataFrames
using Distributed


println("load library")

seed = MersenneTwister(1234)
n_cells = [50, 500, 1000, 5000, 10000, 20000, 30000, 40000, 50000, 70000,
          100000, 150000, 200000, 300000, 400000, 500000]

times_cell = []
dinamica_final = []
n_cell_alive_final = []
Set_mut_final = []
CA_subpop_final = []
α_subpop_final = []
G_state_final = []
Tree_mut_final = []
Tree_fil_final =[]
Newick_final = []
allocated_tot = []
time_dynamics = []

@distributed for i in 1:50
    println("runned configuration -> ", i ," in 3D mode -> contact")
    g_meta = spatial_graph(1000, 1000, seed, dim = 3, n_cell=1)
    println("simulation...")
    #df, G, n_cell_alive, set_mut, Gs_conf, CA_subpop, α_subpop, times =
    value, times, bytes, meta = @timed J_Space.simulate_evolution_2(g_meta,
                                                                   300.0,
                                                                   0.4,
                                                                   0.01,
                                                                   0.0,
                                                                   0.0,
                                                                   0.2,
                                                                   0.1,
                                                                   "contact",
                                                                   seed)
    if value[3][end] < 1
        println("Simulation with less than 1000 cells: simulation discarded")
        continue
    end
    println("singola simulazione salvata")
    push!(G_state_final, copy(value[2]))
    push!(dinamica_final, value[1])
    push!(n_cell_alive_final, value[3])
    push!(Set_mut_final, value[4])
    push!(CA_subpop_final, value[6])
    push!(α_subpop_final, value[7])
    push!(times_cell, value[8])
    push!(allocated_tot, bytes)
    push!(time_dynamics, times)
end

println("SAVE final...")
for l in 1:length(G_state_final)
    if Sys.iswindows()
        mkpath(".\\Experiments\\Experiment_3D\\Fileoutput\\")
        mkpath(".\\Experiments\\Experiment_3D\\Plot\\")
        CSV.write(".\\Experiments\\Experiment_3D\\Fileoutput" *
                  "\\Dynamics-$l.csv",
                  dinamica_final[l],
                  delim = ",")
        CSV.write(".\\Experiments\\Experiment_3D\\Fileoutput" *
                  "\\n_cell_alive_final-$l.csv",
                  Tables.table(n_cell_alive_final[l]),
                  header=false)
        savegraph(".\\Experiments\\Experiment_3D\\Plot" *
                  "\\Final_conf-$l.mg",G_state_final[l])
        List_driver = DataFrame(Driver = Set_mut_final[l],
                                Fitness = α_subpop_final[l])
        CSV.write(".\\Experiments\\Experiment_3D" *
                  "\\Fileoutput\\DriverList-$l.csv",
                  List_driver,delim = ",")
        CSV.write(".\\Experiments\\Experiment_3D\\Fileoutput" *
                  "\\CA_subpop_final-$l.csv",
                  Tables.table(CA_subpop_final[l]),
                  header=false)
        Time_Dynamic = DataFrame(N_cell = n_cells[1:length(times_cell[l])],
                                 T_real = times_cell[l])
        CSV.write(".\\Experiments\\Experiment_3D" *
                  "\\Fileoutput\\Time_Dynamic-$l.csv",
                  Time_Dynamic, delim = ",")
        CSV.write(".\\Experiments\\Experiment_3D\\Fileoutput" *
                  "\\allocated-$l.csv",
                  Tables.table(allocated_tot[l]),
                  header=false)
    else
        mkpath("./Experiments/Experiment_3D/Fileoutput/")
        mkpath("./Experiments/Experiment_3D/Plot/")
        CSV.write("./Experiments/Experiment_3D/Fileoutput" *
                  "/Dynamics-$l.csv",
                  dinamica_final[l],
                  delim = ",")
        CSV.write("./Experiments/Experiment_3D/Fileoutput" *
                  "/n_cell_alive_final-$l.csv",
                  Tables.table(n_cell_alive_final[l]),
                  header=false)
        savegraph("./Experiments/Experiment_3D/Plot" *
                  "/Final_conf-$l.mg",
                  G_state_final[l])
        List_driver = DataFrame(Driver = Set_mut_final[l],
                                Fitness = α_subpop_final[l])
        CSV.write("./Experiments/Experiment_3D/Fileoutput" *
                  "/DriverList-$l.csv",
                  List_driver,delim = ",")
        CSV.write("./Experiments/Experiment_3D/Fileoutput" *
                  "/CA_subpop_final-$l.csv",
                  Tables.table(CA_subpop_final[l]),
                  header=false)
        Time_Dynamic = DataFrame(N_cell = n_cells[1:length(times_cell[l])],
                                 T_real = times_cell[l])
        CSV.write("./Experiments/Experiment_3D" *
                  "/Fileoutput/Time_Dynamic-$l.csv",
                  Time_Dynamic, delim = ",")
        CSV.write("./Experiments/Experiment_3D/Fileoutput" *
                  "/allocated-$l.csv",
                  Tables.table(allocated_tot[l]),
                  header=false)
    end
end
