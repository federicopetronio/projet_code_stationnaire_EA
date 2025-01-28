import numpy as np
import random
import sys
import time

# Define the bounds for each gene
gene_bounds = {
    "B_max": (0, 300e-4),   # max magnetic field
    "beta_mag": (0, 4),     # magnetic field gradient
    "Q_mgs": (0, 5),        # charge-to-mass ratio
    "L_ch": (0, 0.024),     # channel length
    "d_ave": (0, 0.056),    # average diameter
    "Delta_r": (0, 0.015)   # radial difference
}

# Number of individuals in the population
population_size = 100
# Number of genes
num_genes = len(gene_bounds)

# Parameters for the genetic algorithm
num_generations = 50
mutation_rate = 0.1  # Probability of mutating a gene
crossover_rate = 0.5  # Probability of crossover between two parents
elitism_rate = 0.1    # Fraction of the best individuals to carry over

# Randomly generate an individual
def random_individual():
    return {gene: np.random.uniform(low, high) for gene, (low, high) in gene_bounds.items()}

# Create the initial population
def generate_initial_population():
    return [random_individual() for _ in range(population_size)]

# The test function now takes individual gene values as parameters

# Display a loading bar for the current generation
def show_progress(current, total):
    bar_length = 40
    progress = float(current) / total
    block = int(round(bar_length * progress))
    text = f"\rProgress: [{'#' * block + '-' * (bar_length - block)}] {current}/{total} individuals"
    sys.stdout.write(text)
    sys.stdout.flush()

# Selection: Roulette-wheel selection based on fitness
def select(population, fitnesses):
    total_fitness = sum(fitnesses)
    selection_probs = [fitness / total_fitness for fitness in fitnesses]
    return random.choices(population, weights=selection_probs, k=2)

# Crossover: Single-point crossover
def crossover(parent1, parent2):
    if random.random() < crossover_rate:
        crossover_point = random.randint(1, num_genes-1)
        child1 = {**parent1, **dict(list(parent2.items())[crossover_point:])}
        child2 = {**parent2, **dict(list(parent1.items())[crossover_point:])}
        return child1, child2
    else:
        return parent1, parent2

# Mutation: Randomly change a gene
def mutate(individual):
    for gene, (low, high) in gene_bounds.items():
        if random.random() < mutation_rate:
            individual[gene] = np.random.uniform(low, high)
    return individual

# Main genetic algorithm function
def genetic_algorithm():
    population = generate_initial_population()

    for generation in range(num_generations):
        # Evaluate each individual in the population
        fitnesses = []
        for idx, individual in enumerate(population):
            # Pass individual genes to the test function
            fitness = test(
                individual["B_max"],
                individual["beta_mag"],
                individual["Q_mgs"],
                individual["L_ch"],
                individual["d_ave"],
                individual["Delta_r"]
            )
            fitnesses.append(fitness)
            show_progress(idx + 1, population_size)

        # Sort population based on fitness
        sorted_population = [x for _, x in sorted(zip(fitnesses, population), key=lambda pair: pair[0], reverse=True)]

        # Elitism: Carry over the best individuals
        num_elite = int(elitism_rate * population_size)
        new_population = sorted_population[:num_elite]

        # Generate new population through selection, crossover, and mutation
        while len(new_population) < population_size:
            parent1, parent2 = select(sorted_population, fitnesses)
            child1, child2 = crossover(parent1, parent2)
            new_population.append(mutate(child1))
            if len(new_population) < population_size:
                new_population.append(mutate(child2))

        population = new_population

        # Output best result in each generation
        best_individual = sorted_population[0]
        best_fitness = max(fitnesses)
        print(f"\nGeneration {generation+1}: Best fitness = {best_fitness}")

    # Return the best individual from the final generation
    return sorted_population[0]

# Run the genetic algorithm
best_solution = genetic_algorithm()
print("Best solution found:")
print(best_solution)
