import optimizer_utils as opt
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.ticker as ticker

def plot_order(order, reg_dict, src_dict, cal_dict, mark_idx=[], \
    reg_labels=False, mark_cals=True, plot_lst=True):
    Time = opt.get_time_detail(order, reg_dict, src_dict)
    fig, ax = plt.subplots()
    t_wait = 0
    t_slew = 0
    t_obs = 0
    reg_label = []
    for i, t in enumerate(Time):
        # Wait time
        ax.bar(i, t[1], color='k')
        # Slew time
        ax.bar(i, t[2], bottom=t[1], color='y')
        # Obs time
        ax.bar(i, t[3], bottom=t[1]+t[2], color='c')
        if mark_cals and (t[0] in cal_dict.keys()):
            # Calibrators
            ax.plot(i, t[1]+t[2]+t[3]+0.1, marker='v', color='r')
        if i in mark_idx:
            ax.plot(i, t[1]+t[2]+t[3]+0.1, marker='P', color='r')
        t_wait += t[1]
        t_slew += t[2]
        t_obs += t[3]
        reg_label.append(t[0])
    ax.set_xlim(-1, len(Time))
    ax.set_ylabel('Time [h]')
    ax.set_xlabel('Regions')
    ax.set_title(f'Total time: {t_wait + t_slew + t_obs:.2f} h\n' + \
        f'Wait time: {t_wait:.2f} h, Slew time: {t_slew:.2f} h, Obs. time: {t_obs:.2f} h')
    legends = [ Line2D([], [], color='k', lw=5, label='Wait time'),
                Line2D([], [], color='y', lw=5, label='Slew time'),
                Line2D([], [], color='c', lw=5, label='Obs. time')]
    if reg_labels:
        ax.set_xticks(list(range(len(Time))), labels=reg_label)
        ax.tick_params(axis='x', labelsize=6, rotation=45, labelrotation=90)
    if plot_lst:
        ax2 = ax.twinx()
        lst = [t[4] for t in Time]
        ax2.plot([i for i in range(len(Time))], lst, color='r', lw=1.5)
        ax2.yaxis.set_major_locator(ticker.MultipleLocator(24))
        ax2.yaxis.set_minor_locator(ticker.MultipleLocator(12))
        ax2.set_ylabel('LST [h]')
        ax2.set_ylim(bottom=lst[0])
        legends.append(Line2D([], [], color='r', lw=1.5, label='Current LST'))
    if mark_cals:
        legends.append(Line2D([], [], linestyle='None', color='r', marker='v', label='Calibrator'))
    ax.legend(handles=legends, loc='upper left')
    return fig, ax