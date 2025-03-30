"""
@author: Yi-Ling HWONG
"""
import matplotlib as mpl
#mpl.rcParams['figure.dpi'] = 300
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.text as mtext
from applications.utils import *

class LegendTitle(object):
    def __init__(self, text_props=None):
        self.text_props = text_props or {}
        super(LegendTitle, self).__init__()

    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        x0, y0 = handlebox.xdescent, handlebox.ydescent
        title = mtext.Text(x0, y0, orig_handle, usetex=True, **self.text_props)
        handlebox.add_artist(title)
        return title

def get_profile(config,var):

    file_path = "../data/qrad_profile_start/" + config + "_" +var+".csv"
    profile = read_csv_to_list(file_path)

    return profile

def get_timeseries_w(config,end_level_timeseries,plot_moving_average,movingavg_window):

    file_path = "../data/w_timeseries/" + config + "_" + end_level_timeseries + ".csv"
    timeseries = read_csv_to_list(file_path)

    if plot_moving_average == True:
        timeseries = movingaverage(timeseries, movingavg_window, "same")

    return timeseries

def get_timeseries_pw(config):

    file_path = "../data/pw_range_timeseries/" + config + ".csv"
    timeseries = read_csv_to_list(file_path)

    return timeseries

def plot_figure(qrad_profile_all, p_all, z_all, plot_labels, plot_pressure_flag, end_level_profile, line_colors, timeseries_pw, timeseries_w,plot_labels_ts,line_colors_w):

    title_fontsize = 18
    label_fontsize = 16
    tick_fontsize = 16
    legend_fontsize = 14.5

    xlim_min = -2.5
    xlim_max = -0.6
    time_lim_min = 0
    time_lim_max = 39

    # --------- GET TIME AXIS
    # PW RANGE
    len_times = len(timeseries_pw[0])
    times_pw = [n/8 for n in range(len_times)]

    # W
    len_times = len(timeseries_w[0])
    times_w = [n for n in range(len_times)]

    #---------------

    fig = plt.figure(figsize=(12, 7))
    gs = gridspec.GridSpec(2,3)

    linestyles = ["dashed","dashed","dashed","dashed","solid","solid","solid","solid"]

    """
    QRAD EARLY PROFILES FOR AGGREGATED CASES
    """
    ax = fig.add_subplot(gs[:2,0])

    for n in range(0,8):
        if plot_pressure_flag == True:
            ax.plot(qrad_profile_all[n][:end_level_profile],p_all[n][:end_level_profile], label=plot_labels[n], linestyle=linestyles[n],color=line_colors[n])

            y_descrip = "P"
            y_unit = "hPa"
            ax.set_ylim([p_all[n][:end_level_profile][-1], p_all[n][:end_level_profile][0]])
            ax.set_ylim(ax.get_ylim()[::-1])  # invert y axis when y axis is pressure

        elif plot_pressure_flag == False:
            ax.plot(qrad_profile_all[n][:end_level_profile],z_all[n][:end_level_profile], label=plot_labels[n], linestyle=linestyles[n],color=line_colors[n])

            y_descrip = "z"
            y_unit = "m"
            ax.set_ylim([z_all[n][:end_level_profile][0], z_all[n][:end_level_profile][-1]])

        ax.set_xlim([xlim_min, xlim_max])
        ax.legend(loc="upper left", fontsize=legend_fontsize)
        ax.set_title(r"(a) $\overline{Q}_{rad}$, aggregated, day 0-10", fontsize=title_fontsize)
        ax.tick_params(axis='x', labelsize=tick_fontsize)
        ax.tick_params(axis='y', labelsize=tick_fontsize)

        plt.xlabel(r"$\overline{Q}_{rad}$ [K day$^{-1}$]", fontsize=label_fontsize)
        plt.ylabel(y_descrip + ' [' + y_unit + ']', fontsize=label_fontsize)
        plt.grid(True)

    # get handles and labels
    handles, labels = ax.get_legend_handles_labels()
    handles.insert(0,"L128")
    labels.insert(0, "")
    handles.insert(5,"L256")
    labels.insert(5,"")

    # add legend to plot
    ax.legend(handles, labels, fontsize=legend_fontsize,handler_map={str: LegendTitle({'fontsize': legend_fontsize+4})})


    """
    PW RANGE TIMESERIES
    """
    ax = fig.add_subplot(gs[0,1:])

    for n in range(0, 8):
        ax.plot(times_pw, timeseries_pw[n], label=plot_labels[n], linestyle=linestyles[n], color=line_colors[n])
        ax.set_xlim([time_lim_min, time_lim_max])

        ax.tick_params(axis='x', labelsize=tick_fontsize)
        ax.tick_params(axis='y', labelsize=tick_fontsize)
        plt.ylabel(r"PW$_{moist-dry}$ [mm]", fontsize=label_fontsize)
        ax.set_title("(b) " + r"PW range time-series", fontsize=title_fontsize)
        ax.set_ylim([0.0,45])
        ax.set_xlim(left=0.0)
        ax.xaxis.set_tick_params(labelbottom=False)

        plt.grid(True)


    """
    WSUB TIMESERIES
    """
    ax = fig.add_subplot(gs[1,1:])

    for n in range(0,8):
        ax.plot(times_w, timeseries_w[n],label=plot_labels_ts[n], linestyle=linestyles[n],color=line_colors_w[n])
        ax.set_xlim([time_lim_min, time_lim_max])
        ax.tick_params(axis='x', labelsize=tick_fontsize)
        ax.tick_params(axis='y', labelsize=tick_fontsize)
        ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
        plt.xlabel("Time [day]", fontsize=label_fontsize)
        plt.ylabel(r"$\overline{w}_{dry,PBL}$ [m s$^{-1}$]", fontsize=label_fontsize)
        ax.set_title("(c) " + r"$\overline{w}_{dry,PBL}$ time-series", fontsize=title_fontsize)
        plt.axhline(y=0, color="darkred", linestyle="dotted",linewidth=1.5)
        ax.set_ylim([-0.01, 0.013])

        plt.grid(True)

    # get handles and labels
    handles, labels = plt.gca().get_legend_handles_labels()

    # specify order of items in legend
    order = [0, 4, 1, 5, 2, 6, 3, 7]

    # add legend to plot
    ax.legend([handles[idx] for idx in order], [labels[idx] for idx in order], fontsize=legend_fontsize-2.5, ncol=4, loc="upper center")

    #-----------------------------------------

    fig.subplots_adjust(left=0.085,
                        bottom=0.12,
                        right=0.99,
                        top=0.92,
                        wspace=0.43,
                        hspace=0.25)

    """
    wspace and hspace specify the space reserved between Matplotlib subplots. They are the fractions of axis width and height, respectively.
    left, right, top and bottom parameters specify four sides of the subplots’ positions. They are the fractions of the width and height of the figure.
    top and bottom should add up to 1.0
    """

    plt.show()


def main():

    # QRAD EARLY STAGE PROFILES & PW TIMESERIES
    configs = ["alpha0_L128", "alpha0.001_L128", "alpha0.002_L128","alpha0.005_L128",
               "alpha0_L256", "alpha0.01_L256", "alpha0.02_L256", "alpha0.05_L256"]

    plot_labels = [r"${\alpha}$=0", r"${\alpha}$=0.001", r"${\alpha}$=0.002",r"${\alpha}$=0.005",
                   r"${\alpha}$=0", r"${\alpha}$=0.01", r"${\alpha}$=0.02", r"${\alpha}$=0.05"]

    # W TIMESERIES
    configs_ts = ["alpha0_L128","alpha0.01_L128","alpha0.02_L128","alpha0.05_L128",
               "alpha0_L256",  "alpha0.01_L256", "alpha0.02_L256", "alpha0.05_L256"]
    plot_labels_ts = ["L128_"+r"${\alpha}$0", "L128_"+r"${\alpha}$0.01", "L128_"+r"${\alpha}$0.02","L128_"+r"${\alpha}$0.05",
                   "L256_"+r"${\alpha}$0", "L256_"+r"${\alpha}$0.01", "L256_"+r"${\alpha}$0.02","L256_"+r"${\alpha}$0.05"]

    # PLOT OPTIONS
    plot_pressure_flag = True # set to True for pressure on Y-axis, False for height on Y-axis
    end_level_profile = 25 # 9 for 850hPa, 11 for 800 hPa (2km), 21 for 500 hPa, 25 for 400 hPa, 37 for 200 hPa, 41 for 150hPa, 64 for total (for 1d file), 53 for total
    plot_moving_average = True
    movingavg_window = 7
    end_level_timeseries = "1km" # 1, 2, 3, 5km

    # LINE COLOURS
    line_colors = get_line_colors("gnuplot", int(len(configs))+1)
    line_colors = line_colors[1:]
    line_colors_w_L128 = get_line_colors("YlGnBu_r", int(len(configs_ts) / 2))
    line_colors_w_L256 = get_line_colors("YlGnBu_r", int(len(configs_ts) / 2))
    line_colors_w = line_colors_w_L128 + line_colors_w_L256


    # GET DATA
    qrad_profile_all = []
    p_all = []
    z_all = []

    for config in configs:

        qrad = get_profile(config,"qrad")
        p = get_profile(config,"p")
        z = get_profile(config,"z")

        qrad_profile_all.append(qrad)
        p_all.append(p)
        z_all.append(z)

    timeseries_pw = []
    for config in configs:
        timeseries = get_timeseries_pw(config)
        timeseries_pw.append(timeseries)

    timeseries_w = []
    for config in configs_ts:

        timeseries = get_timeseries_w(config,end_level_timeseries,plot_moving_average,movingavg_window)
        timeseries_w.append(timeseries)


    plot_figure(qrad_profile_all, p_all, z_all, plot_labels, plot_pressure_flag, end_level_profile, line_colors,
                timeseries_pw, timeseries_w,plot_labels_ts,line_colors_w)

if __name__ == "__main__":

    main()