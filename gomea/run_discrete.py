import gomea
import sys

n_vars = 60
k = 5
n_pop = 100
if(len(sys.argv) == 2):
    n_vars = int(sys.argv[1])
elif(len(sys.argv) == 3):
    n_vars = int(sys.argv[1])
    k = int(sys.argv[2])
elif(len(sys.argv) == 4):
    n_vars = int(sys.argv[1])
    k = int(sys.argv[2])
    n_pop = int(sys.argv[3])

print("Running discrete GOMEA with n_vars = {}, k = {}, n_pop = {}".format(n_vars,k,n_pop))
f = gomea.fitness.DeceptiveTrapFunction(n_vars,trap_size=k)
lm = gomea.linkage.LinkageTree()
gd = gomea.DiscreteGOMEA(fitness=f,linkage_model=lm,max_number_of_populations=1,base_population_size=n_pop)
out = gd.run()

out.printAllStatistics()
out.printFinalStatistics()