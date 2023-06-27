import gomea
f = gomea.fitness.DeceptiveTrapFunction(60,trap_size=5)
lm = gomea.linkage.LinkageTree()
gd = gomea.DiscreteGOMEA(fitness=f,linkage_model=lm,max_number_of_populations=1,base_population_size=100)
out = gd.run()

out.printAllStatistics()
out.printFinalStatistics()