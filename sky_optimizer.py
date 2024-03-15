import random
from tqdm import trange

class SkyOptimizer:
    def __init__(self, regions_data, init_order, **kwargs):
        self.regions = regions_data
        self.order = init_order
        self.optimizer = None
    
    def set_optimizer(self, optimizer, **kwargs):
        match optimizer:
            case 'genetic':
                self.pop_size = kwargs['pop_size']
                self.offsp_size = kwargs['offsp_size']
                self.gen_num = kwargs['gen_num']
                self.mut_prob = kwargs['mut_prob']
                self.mut_rate = kwargs['mut_rate']
                self.optimizer = self.opt_genetic
                self.set_pop_generator(kwargs['pop_gen'])
                self.set_par_selector(kwargs['selector'])
                self.set_crossover(kwargs['crossover_op'])
                self.set_mutop(kwargs['mut_op'])
                self.set_survop(kwargs['surv_op'])
    
    def set_pop_generator(self, generator):
        match generator:
            case 'random':
                self.new_pop = self.popgen_random
    
    def set_par_selector(self, selector):
        match selector:
            case 'natural':
                self.selection = self.parsel_natural
    
    def set_crossover(self, crossover_op):
        match crossover_op:
            case 'single_point':
                self.crossover = self.crossover_single
    
    def set_mutop(self, mut_op):
        match mut_op:
            case 'permutate':
                self.mutation = self.mutop_permutate

    def set_survop(self, surv_op):
        match surv_op:
            case 'replace':
                self.survival = self.survop_replace

    def run_optimizer(self, print_output=True):
        self.optimizer()

    # Population generators
    def popgen_random(self):
        population = []
        for _ in range(self.pop_size):
            new_pop = self.order.copy()
            random.shuffle(new_pop)
            population.append(new_pop)
        return population

    # Parent selectors
    def parsel_natural(self, population, num):
        parents = []
        for _ in range(num):
            par1 = random.choice(population)
            par2 = random.choice(population)
            parents.append((par1, par2))
        return parents

    # Crossover operators
    def crossover_single(self, parents):
        children = []
        for p1, p2 in parents:
            r = random.randint(1, len(p1) - 1)
            new_child = p1[:r]
            for s in p2:
                if s not in new_child:
                    new_child.append(s)
            children.append(new_child)
        return children

    # Survival selection operators
    def survop_replace(self, pop, children):
        return children

    # Mutation operators
    def mutop_permutate(self, children):
        for c in children:
            if self.mut_prob > random.random():
                n = len(c)
                mut_num = int(self.mut_rate * n / 2)
                for _ in range(mut_num):
                    i = random.randint(0, n - 1)
                    j = random.randint(0, n - 1)
                    c[i], c[j] = c[j], c[i]

    def opt_genetic(self):
        pop = self.new_pop()
        pbar = trange(self.gen_num)
        for gen in pbar:
            parents = self.selection(pop, self.offsp_size)
            children = self.crossover(parents)
            children = self.mutation(children)
            pop = self.survival(pop, children)
        return pop