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
                # self.parent_num = kwargs['parent_num']
                self.offsp_size = kwargs['offsp_size']
                self.mut_prob = kwargs['mut_prob']
                self.gen_num = kwargs['gen_num']
                self.mut_prob = kwargs['mut_prob']
                self.optimizer = self.opt_genetic
                self.set_pop_generator(kwargs['pop_gen'])
                self.set_par_selector(kwargs['selector'])
                self.set_crossover(kwargs['crossover'])
                self.set_survop(kwargs['surv_op'])
    
    def set_pop_generator(self, generator):
        match generator:
            case 'random':
                self.new_pop = self.popgen_random
    
    def set_par_selector(self, selector):
        match selector:
            case 'natural':
                self.select_par = self.parsel_natural
    
    def set_crossover(self, crossover_op):
        match crossover_op:
            case 'single_point':
                self.crossover = self.crossover_single

    def set_survop(self, surv_op):
        match surv_op:
            case 'replace':
                self.survop = self.survop_replace

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

    def opt_genetic(self):
        pop = self.new_pop()
        pbar = trange(self.gen_num)
        for gen in pbar:
            parents = self.select_par(pop, self.offsp_size)
            children = self.crossover(parents)
            pop = self.survop(pop, children)
        return pop