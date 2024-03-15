import random
import copy
from tqdm import trange

class SkyOptimizer:
    def __init__(self, init_order, obj_func, **kwargs):
        self.order = init_order
        self.obj_func = obj_func
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
                self.set_mutop(kwargs['mutation_op'])
                self.set_survop(kwargs['survival_op'])
    
    def set_pop_generator(self, generator):
        match generator:
            case 'random':
                self.new_pop = self.popgen_random
            case 'identical':
                self.new_pop = self.popgen_identical
    
    def set_par_selector(self, selector):
        match selector:
            case 'natural':
                self.selection = self.parsel_natural
            case 'rank':
                self.selection = self.parsel_rank
    
    def set_crossover(self, crossover_op):
        match crossover_op:
            case 'single_point':
                self.crossover = self.crossover_single
            case 'clone':
                self.crossover = self.crossover_clone
    
    def set_mutop(self, mut_op):
        match mut_op:
            case 'permutate':
                self.mutation = self.mutop_permutate
            case 'swap':
                self.mutation = self.mutop_swap

    def set_survop(self, surv_op):
        match surv_op:
            case 'replace':
                self.survival = self.survop_replace
            case 'evolution':
                self.survival = self.survop_evol

    def run_optimizer(self, print_output=True):
        pop = self.optimizer()
        ranked_pop = sorted(pop, key=self.obj_func)
        self.order = ranked_pop[0]
        return self.order

    # Population generators
    def popgen_random(self):
        population = []
        for _ in range(self.pop_size):
            new_pop = self.order.copy()
            random.shuffle(new_pop)
            population.append(new_pop)
        return population
    
    def popgen_identical(self):
        population = []
        for _ in range(self.pop_size):
            population.append(self.order.copy())
        return population

    # Parent selectors
    def parsel_natural(self, population, num):
        parents = []
        for _ in range(num):
            par1 = random.choice(population)
            par2 = random.choice(population)
            parents.append((par1, par2))
        return parents

    def parsel_rank(self, population, num):
        pop_rank = sorted(population, key=self.obj_func)
        parents = []
        for i in range(1, 2*num, 2):
            par1 = pop_rank[i-1]
            par2 = pop_rank[i]
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

    def crossover_clone(self, parents):
        children = []
        for p1, p2 in parents:
            children.append(p1.copy())
        return children

    # Survival selection operators
    def survop_replace(self, pop, children):
        return children

    def survop_evol(self, pop, children):
        pop.extend(children)
        # print(len(pop), len(children))
        pop_rank = sorted(pop, key=self.obj_func)
        # for i in children:
        #     print(self.obj_func(i), i)
        return pop_rank[:self.pop_size]

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
        return children

    def mutop_swap(self, children):
        for c in children:
            if self.mut_prob > random.random():
                n = len(c)
                i = random.randint(0, n - 2)
                # print('mutation', i, i+1)
                c[i], c[i+1] = c[i+1], c[i]
        return children

    def opt_genetic(self):
        pop = self.new_pop()
        pbar = trange(self.gen_num)
        for gen in pbar:
            parents = self.selection(pop, self.offsp_size)
            children = self.crossover(parents)
            children = self.mutation(children)
            pop = self.survival(pop, children)
            min_sol = min(pop, key=self.obj_func)
            total_time = self.obj_func(min_sol)
            pbar.set_postfix({'total_time':f'{total_time:.2f}'})
        return pop