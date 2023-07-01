#_________________ Función 1
import numpy
def build_population(N, p):
    population = []
    for i in range(N):
        allele1 = "A"
        if numpy.random.rand() > p:
            allele1 = "a"
        allele2 = "A"
        if numpy.random.rand() > p:
            allele2 = "a"
        population.append((allele1, allele2))
    return population

#_________________ Función 2
def compute_frequencies(population):
    AA = population.count(("A", "A"))
    Aa = population.count(("A", "a"))
    aA = population.count(("a", "A"))
    aa = population.count(("a", "a"))
    return({"AA": AA, "aa": aa, "Aa": Aa, "aA": aA})

#_________________ Función 3
import numpy as np
def reproduce_population(population):
    new_generation = []
    N = len(population)
    for i in range(N):
        dad = np.random.randint(N)
        mom = np.random.randint(N)
        chr_mom = np.random.randint(2)
        offspring = (population[mom][chr_mom], population[dad][1 - chr_mom])
        new_generation.append(offspring)
    return new_generation

#_________________ Función 4
def simulate_drift(N, p):
    my_pop = build_population(N, p)
    fixation = False
    num_generations = 0
    while fixation == False:
        genotype_counts = compute_frequencies(my_pop)
        if (genotype_counts["AA"] == N or genotype_counts["aa"] == N):
            print("An allele reached fixation at generation", num_generations)
            print("The genotype counts are")
            print(genotype_counts)
            fixation == True
            break
        my_pop = reproduce_population(my_pop)
        num_generations += 1
    return num_generations, genotype_counts
