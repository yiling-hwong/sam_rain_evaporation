"""
@author: Yi-Ling HWONG
"""
import matplotlib as mpl
#mpl.rcParams['figure.dpi'] = 300
import matplotlib.pyplot as plt
from applications.utils import *

def get_timeseries(config,agg_index):

    file_path = "../data/agg_index_timeseries/"+config+"_"+agg_index+".csv"
    timeseries = read_csv_to_list(file_path)

    return timeseries

def plot_timeseries(agg_index,agg_index_timeseries,plot_labels,line_colors):

    title_fontsize = 17
    label_fontsize = 15
    tick_fontsize = 14
    legend_fontsize = 12

    xlim_min = 0
    xlim_max = 39

    # --------- GET TIME AXIS
    times = [int(n) for n in range(len(agg_index_timeseries[0]))]

    fig = plt.figure(figsize=(8, 6))

    """
    L128
    """
    ax = fig.add_subplot(2, 1, 1)

    for n in range(0,8):

        ax.plot(times,agg_index_timeseries[n], label=plot_labels[n], color=line_colors[n])

        ax.set_xlim([xlim_min, xlim_max])
        ax.set_title("(a) L128", fontsize=title_fontsize)
        ax.tick_params(axis='x', labelsize=tick_fontsize)
        ax.tick_params(axis='y', labelsize=tick_fontsize)

        plt.xlabel("Time [day]", fontsize=label_fontsize)

        if agg_index == "org_pw":
            plt.ylabel(r"$\hat{\sigma}^{2}_{pw}$ [mm]", fontsize=label_fontsize)
        if agg_index == "crh_var":
            plt.ylabel(r"${\sigma}_{crh}$", fontsize=label_fontsize)

        plt.grid(True)

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles=handles, labels=labels, bbox_to_anchor=(1.31, 0.45), loc="center right", ncol=1,fontsize=legend_fontsize)


    """
    L256
    """
    ax = fig.add_subplot(2, 1, 2)

    for n in range(8,16):
        ax.plot(times, agg_index_timeseries[n], label=plot_labels[n], color=line_colors[n-8])

        ax.set_xlim([xlim_min, xlim_max])
        ax.set_title("(b) L256", fontsize=title_fontsize)
        ax.tick_params(axis='x', labelsize=tick_fontsize)
        ax.tick_params(axis='y', labelsize=tick_fontsize)

        plt.xlabel("Time [day]", fontsize=label_fontsize)

        if agg_index == "org_pw":
            plt.ylabel(r"$\hat{\sigma}^{2}_{pw}$ [mm]", fontsize=label_fontsize)
        if agg_index == "crh_var":
            plt.ylabel(r"${\sigma}_{crh}$", fontsize=label_fontsize)

        plt.grid(True)

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles=handles, labels=labels, bbox_to_anchor=(1.31, 0.45), loc="center right", ncol=1,fontsize=legend_fontsize)


    #-----------------------------------------
    fig.subplots_adjust(left=0.1,
                        bottom=0.1,
                        right=0.78,
                        top=0.92,
                        wspace=0.08,
                        hspace=0.55)

    """
    wspace and hspace specify the space reserved between Matplotlib subplots. They are the fractions of axis width and height, respectively.
    left, right, top and bottom parameters specify four sides of the subplots’ positions. They are the fractions of the width and height of the figure.
    top and bottom should add up to 1.0
    """

    plt.show()

def main():

    agg_index = "org_pw" # org_pw or crh_var

    configs = ["alpha0_L128", "alpha0.001_L128", "alpha0.002_L128","alpha0.005_L128",
               "alpha0.01_L128", "alpha0.02_L128", "alpha0.05_L128","alpha1_L128",
               "alpha0_L256", "alpha0.01_L256", "alpha0.02_L256", "alpha0.05_L256",
               "alpha0.1_L256", "alpha0.2_L256", "alpha0.5_L256", "alpha1_L256"]

    plot_labels = [r"${\alpha}$=0", r"${\alpha}$=0.001", r"${\alpha}$=0.002",r"${\alpha}$=0.005",
                   r"${\alpha}$=0.01", r"${\alpha}$=0.02", r"${\alpha}$=0.05",r"${\alpha}$=1",
                   r"${\alpha}$=0", r"${\alpha}$=0.01", r"${\alpha}$=0.02", r"${\alpha}$=0.05",
                   r"${\alpha}$=0.1", r"${\alpha}$=0.2", r"${\alpha}$=0.5", r"${\alpha}$=1"]

    line_colors = get_line_colors("gnuplot", int(len(configs)/2)+1)  # YlOrRd_r,Dark2,rainbow,jet,turbo,gnuplot,gnuplot2,tab10,Accent
    line_colors = line_colors[1:]


    agg_index_timeseries = []

    for config in configs:

        timeseries = get_timeseries(config,agg_index)
        agg_index_timeseries.append(timeseries)

    print (len(agg_index_timeseries))

    plot_timeseries(agg_index,agg_index_timeseries,plot_labels,line_colors)

if __name__ == "__main__":

    main()