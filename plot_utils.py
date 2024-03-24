import optimizer_utils as opt
import matplotlib.pyplot as plt

def plot_order(order, reg_dict, src_dict, cal_dict, path=None, **kwargs):
    Time = opt.get_time_detail(order, reg_dict, src_dict)
    fig, ax = plt.subplots()
    for i, t in enumerate(Time):
        ax.bar(i, t[1], color='k', align='edge')
        ax.bar(i, t[2], bottom=t[1], color='y', align='edge')
        ax.bar(i, t[3], bottom=t[1]+t[2], color='c', align='edge')
        if t[0] in cal_dict.keys():
            ax.plot(i+0.25, t[1]+t[2]+t[3]+0.1, marker='*', color='r')
    ax.set_xlim(-0.5, len(Time)+0.5)
    return fig, ax