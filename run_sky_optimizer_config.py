# Sky optimizer config file

# Seed for the RNG
random_seed = 42

# List of calibrator regions
cal_list = ['136', '137', '138', '139', '140', '141', '142']

# Options for the genetic algorithm optimizer
optimizer_kwargs = {'pop_size':100,
                    'offsp_size':25,
                    'gen_num':50,
                    'mut_prob':0.9,
                    'mut_rate':0.05,
                    'pop_gen':'identical',
                    'selector':'rank',
                    'crossover_op':'single_point',
                    'mutation_op':'local_perm',
                    'survival_op':'evolution'
                    }

# Run the GA optimizer
run_optimizer = True
# Place calibrators
place_calibrators = True
# Perform permutations around calibrators
permutate_calibrators = True
# Amount of regions to permutate around calibrators
cal_perm_size = 8
# Number of calibrator permutations to try
cal_perm_iter = 10_000
# Perform permutations around peaks of time
permutate_peaks = True
# Amount of regions to permutate around peaks
peak_perm_size = 8
# Number of peak permutations to try
peak_perm_iter = 20_000
# Number of peaks to permutate
peak_max_n = 5
# Attempt filling wait times with nearby regions
fill_wait_times = True
# Plot region times
make_plot = False