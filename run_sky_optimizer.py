#!/usr/bin/env python3

import optimizer_utils as opt
from sky_optimizer import SkyOptimizer
import sys
import pickle
import random
import run_sky_optimizer_config as cfg
from plot_utils import plot_order
import os

random.seed(cfg.random_seed)

date_tag = sys.argv[1]
lst_start = float(sys.argv[2])
sun_jd = float(sys.argv[3])

if sys.argv[1] == 'None':
    date_tag = ''

time_cal = 575. / 3600. # Get from telescope_data

cal_list = cfg.cal_list

base_path = os.path.dirname(os.getcwd()) + '/'
intermediate_results_path = base_path + 'intermediate_results/'
# intermediate_results_path = '../intermediate_results/'
region_optimization_results_file = f"{intermediate_results_path}region_optimizer_{date_tag}.dat"
sky_optimization_results_file = f"{intermediate_results_path}sky_optimizer_{date_tag}.dat"

with open(region_optimization_results_file, 'rb') as file:
    print(f'Data loaded from {region_optimization_results_file}')
    sources = pickle.load(file, encoding='latin-1')
    regions_data = pickle.load(file, encoding='latin-1')

regions = {d['number']: d for d in regions_data}

for k in cal_list:
    regions[k]['order'] = [0]
    regions[k]['sources'] = [regions[k]['sources'][0]]
    regions[k]['obstime'] = time_cal

calibrators = {key: regions[key] for key in cal_list}

region_numbers = list(regions.keys())

start_order = [x for x in region_numbers if x not in cal_list]
random.shuffle(start_order)

start_order = ['85', '22', '70', '87', '42', '90', '117', '29', '6', '40', '76', '116', '62', '64', '96', '60', '61', '113', '79', '20', '110', '126', '63', '124', '81', '125', '49', '8', '39', '24', '107', '21', '51', '10', '37', '0', '5', '133', '1', '15', '132', '91', '114', '93', '92', '130', '129', '94', '115', '48', '65', '97', '111', '80', '34', '50', '35', '52', '36', '109', '3', '86', '102', '68', '9', '55', '53', '25', '71', '59', '14', '74', '28', '118', '57', '75', '131', '16', '77', '2', '32', '47', '18', '33', '66', '19', '82', '67', '98', '99', '69', '11', '83', '108', '123', '38', '54', '101', '84', '106', '121', '56', '13', '41', '4', '26', '73', '58', '27']

def obj_func(order):
    return opt.compute_total_time(regions, order, sources, lst_start)

def val_func(order):
    return opt.check_calibrators(order, regions, sources, calibrators, lst_start)

if cfg.run_optimizer:
    print('Running GA Optimizer')
    optimizer_kwargs = cfg.optimizer_kwargs
    optimizer = SkyOptimizer(start_order, obj_func)
    optimizer.set_optimizer('genetic', **optimizer_kwargs)
    order_opt = optimizer.run_optimizer()
else:
    order_opt = start_order.copy()

if cfg.place_calibrators:
    print('placing calibrator regions')
    while opt.check_calibrators(order_opt, regions, sources, calibrators, lst_start, prnt=False) is False:
        order_opt = opt.place_calibrator2(order_opt, regions, sources, calibrators, lst_start)
        print()

if cfg.permutate_calibrators:
    print('Permutating calibrator regions')
    perm_size = cfg.cal_perm_size
    perm_iter = cfg.cal_perm_iter
    for idx, reg in enumerate(order_opt):
        if reg in calibrators.keys():
            order_opt = opt.local_perm(order_opt, idx, perm_size, perm_iter, obj_func, val_f=val_func)

if cfg.permutate_peaks:
    print('Permutating peaks')
    perm_size = cfg.peak_perm_size
    perm_iter = cfg.peak_perm_iter
    max_n = cfg.peak_max_n
    for _ in range(max_n):
        Time = opt.get_time_detail(order_opt, regions, sources)
        idx = Time.index(max(Time, key=lambda x: x[1]))
        order_opt = opt.local_perm(order_opt, idx, perm_size, perm_iter, obj_func, val_f=val_func)

if cfg.fill_wait_times:
    print('Filling wait times')
    order_opt, added_regs = opt.fill_wait(order_opt, regions, sources, exclude_list=list(calibrators.keys()), output=True)

final_order = order_opt.copy()
print(f"{final_order=}")

final_order_lst = opt.get_lst_obs(final_order, regions, sources, lst_i=lst_start)

with open(sky_optimization_results_file, 'wb') as file:
        print('Saving sky optimization results file')
        pickle.dump(final_order, file, protocol=2)
        pickle.dump(final_order_lst, file, protocol=2)

if cfg.make_plot:
    print('Plotting results')
    fig, ax = plot_order(final_order, regions, sources, calibrators)
    fig.savefig(f'../plots_regions/plot_{date_tag}.png', dpi=200)