import pickle
import json
from telescope_data import time_cal

prefix = 'data/'
suffix = '_full'
cal_list = ['136', '137', '138', '139', '140', '141', '142']

dat_filename = f'{prefix}regopt_results{suffix}.dat'
src_filename = f'{prefix}sources{suffix}.json'
reg_filename = f'{prefix}regions{suffix}.json'
cal_filename = f'{prefix}calibrators{suffix}.json'

with open(dat_filename, 'rb') as file:
    print(f'Data loaded from {dat_filename}')
    sources = pickle.load(file)
    regions = pickle.load(file)

reg_dict = {d['number']: d for d in regions}

print(f'Calibrator regions: {cal_list}')
for k in cal_list:
    reg_dict[k]['order'] = [0]
    reg_dict[k]['sources'] = [reg_dict[k]['sources'][0]]
    reg_dict[k]['obstime'] = time_cal

cal_dict = {key: reg_dict[key] for key in cal_list}

with open(src_filename, 'w') as file:
        print(f'{len(sources)} sources saved to {src_filename}')
        json.dump(sources, file, indent=4, default=int)

with open(reg_filename, 'w') as file:
        print(f'{len(reg_dict)} regions saved to {src_filename}')
        json.dump(reg_dict, file, indent=4)

with open(cal_filename, 'w') as file:
        print(f'{len(cal_dict)} calibrators saved to {src_filename}')
        json.dump(cal_dict, file, indent=4)