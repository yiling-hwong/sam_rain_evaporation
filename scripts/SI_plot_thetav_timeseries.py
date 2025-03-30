"""
@author: Yi-Ling HWONG
"""
import matplotlib as mpl
#mpl.rcParams['figure.dpi'] = 300
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.legend import Legend
from applications.utils import *

def get_timeseries(configs,target):

    timeseries_list = []
    for config in configs:
        file_path = "../data/thetav_timeseries/alpha"+config+"_"+target+".csv"
        timeseries = read_csv_to_list(file_path)

        timeseries_list.append(timeseries)

    return timeseries_list


def plot_timeseries(to_plot,thetav_wet_timeseries_L128,thetav_dry_timeseries_L128,theta_wet_timeseries_L128,theta_dry_timeseries_L128, qv_wet_timeseries_L128,qv_dry_timeseries_L128,
                             thetav_wet_timeseries_L256, thetav_dry_timeseries_L256,theta_wet_timeseries_L256, theta_dry_timeseries_L256,
                             qv_wet_timeseries_L256, qv_dry_timeseries_L256):

    f, axarr = plt.subplots(nrows=3, ncols=2, sharex=True, figsize=(10, 8))
    ((ax1, ax2), (ax3, ax4), (ax5,ax6)) = axarr

    xlim_min = 0  # in days
    xlim_max = 20
    suptitle_fontsize = 19
    title_fontsize = 16
    label_fontsize = 15
    tick_fontsize = 12
    legend_fontsize = 13
    sublegend_fontsize = 12

    # --------- Get ylim min and max
    ylim_min_L128 = -1e-3
    ylim_max_L128 = -ylim_min_L128
    ylim_min_L256 = -1e-3
    ylim_max_L256 = -ylim_min_L256

    # --------- GET TIME AXIS AND ZERO ARRAY
    times = [int(n) for n in range(len(thetav_wet_timeseries_L128[0]))]
    zero_arr = np.zeros(len(times))

    """
    1. L128 THETA_V
    """
    ax1.plot(times, thetav_wet_timeseries_L128[0], color="red", linestyle="solid", label="Moist")
    ax1.plot(times, thetav_wet_timeseries_L128[1], color="red", linestyle="dashed")
    ax1.plot(times, thetav_dry_timeseries_L128[0], color="blue", linestyle="solid", label="Dry")
    ax1.plot(times, thetav_dry_timeseries_L128[1], color="blue", linestyle="dashed")

    ax1.set_xlim(xlim_min, xlim_max)
    ax1.set_ylim(ylim_min_L128, ylim_max_L128)
    ax1.set_yticks(np.linspace(ax1.get_ybound()[0], ax1.get_ybound()[1], 5))
    ax1.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    ax1.tick_params(axis='x', labelsize=tick_fontsize)
    ax1.tick_params(axis='y', labelsize=tick_fontsize)
    ax1.yaxis.get_offset_text().set_fontsize(tick_fontsize)
    ax1.grid()
    ax1.legend(loc="upper left", fontsize=legend_fontsize)
    ax1.plot(times, zero_arr, color="black", linestyle="solid", linewidth=1.5)
    ax1.axhline(y=0, color='black', linestyle='-', linewidth=1.5)

    # ----------
    # Additional legend

    custom_lines = [Line2D([0], [0], color="black", lw=1.5, linestyle="solid"),
                    Line2D([0], [0], color="black", lw=1.5, linestyle="dashed"),
                    Line2D([0], [0], color="black", lw=1.5, linestyle="dotted")]

    if to_plot == "agg":
        labels_leg2 = [r"$\alpha=0.001$", r"$\alpha=0.002$"]
    if to_plot == "nonagg":
        labels_leg2 = [r"$\alpha=0.01$", r"$\alpha=0.02$"]
    leg = Legend(ax1, custom_lines, labels_leg2, loc=3, fontsize=sublegend_fontsize, ncol=3)
    ax1.add_artist(leg)

    # -----------

    """
    2. L256 THETA_V
    """
    ax2.plot(times, thetav_wet_timeseries_L256[0], color="red", linestyle="solid", label="Moist")
    ax2.plot(times, thetav_wet_timeseries_L256[1], color="red", linestyle="dashed")
    ax2.plot(times, thetav_dry_timeseries_L256[0], color="blue", linestyle="solid", label="Dry")
    ax2.plot(times, thetav_dry_timeseries_L256[1], color="blue", linestyle="dashed")

    ax2.set_xlim(xlim_min, xlim_max)
    ax2.set_ylim(ylim_min_L256, ylim_max_L256)
    ax2.set_yticks(np.linspace(ax2.get_ybound()[0], ax2.get_ybound()[1], 5))
    ax2.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    ax2.tick_params(axis='x', labelsize=tick_fontsize)
    ax2.tick_params(axis='y', labelsize=tick_fontsize)
    ax2.yaxis.get_offset_text().set_fontsize(tick_fontsize)
    ax2.grid()
    ax2.legend(loc="upper left", fontsize=legend_fontsize)
    ax2.plot(times, zero_arr, color="black", linestyle="solid", linewidth=1.5)
    ax2.axhline(y=0, color='black', linestyle='-', linewidth=1.5)

    # ----------
    # Additional legend

    custom_lines = [Line2D([0], [0], color="black", lw=1.5, linestyle="solid"),
                    Line2D([0], [0], color="black", lw=1.5, linestyle="dashed"),
                    Line2D([0], [0], color="black", lw=1.5, linestyle="dotted")]

    if to_plot == "agg":
        labels_leg2 = [r"$\alpha=0.01$", r"$\alpha=0.02$"]
    if to_plot == "nonagg":
        labels_leg2 = [r"$\alpha=0.1$", r"$\alpha=0.2$"]
    leg = Legend(ax2, custom_lines, labels_leg2, loc=3, fontsize=sublegend_fontsize, ncol=3)
    ax2.add_artist(leg)

    # -----------

    """
    3. L128 THETA
    """
    ax3.plot(times, theta_wet_timeseries_L128[0], color="red", linestyle="solid", label="Moist")
    ax3.plot(times, theta_wet_timeseries_L128[1], color="red", linestyle="dashed")
    ax3.plot(times, theta_dry_timeseries_L128[0], color="blue", linestyle="solid", label="Dry")
    ax3.plot(times, theta_dry_timeseries_L128[1], color="blue", linestyle="dashed")

    ax3.set_xlim(xlim_min, xlim_max)
    ax3.set_ylim(ylim_min_L128, ylim_max_L128)
    ax3.set_yticks(np.linspace(ax3.get_ybound()[0], ax3.get_ybound()[1], 5))
    ax3.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    ax3.tick_params(axis='x', labelsize=tick_fontsize)
    ax3.tick_params(axis='y', labelsize=tick_fontsize)
    ax3.yaxis.get_offset_text().set_fontsize(tick_fontsize)
    ax3.grid()
    ax3.axhline(y=0, color='black', linestyle='-', linewidth=1.5)

    """
    4. L256 THETA
    """
    ax4.plot(times, theta_wet_timeseries_L256[0], color="red", linestyle="solid", label="Moist")
    ax4.plot(times, theta_wet_timeseries_L256[1], color="red", linestyle="dashed")
    ax4.plot(times, theta_dry_timeseries_L256[0], color="blue", linestyle="solid", label="Dry")
    ax4.plot(times, theta_dry_timeseries_L256[1], color="blue", linestyle="dashed")

    ax4.set_xlim(xlim_min, xlim_max)
    ax4.set_ylim(ylim_min_L256, ylim_max_L256)
    ax4.set_yticks(np.linspace(ax4.get_ybound()[0], ax4.get_ybound()[1], 5))
    ax4.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    ax4.tick_params(axis='x', labelsize=tick_fontsize)
    ax4.tick_params(axis='y', labelsize=tick_fontsize)
    ax4.yaxis.get_offset_text().set_fontsize(tick_fontsize)
    ax4.grid()
    ax4.axhline(y=0, color='black', linestyle='-', linewidth=1.5)

    """
    5. L128 QV
    """

    ax5.plot(times, qv_wet_timeseries_L128[0], color="red", linestyle="solid", label="Moist")
    ax5.plot(times, qv_wet_timeseries_L128[1], color="red", linestyle="dashed")
    ax5.plot(times, qv_dry_timeseries_L128[0], color="blue", linestyle="solid", label="Dry")
    ax5.plot(times, qv_dry_timeseries_L128[1], color="blue", linestyle="dashed")

    ax5.set_xlim(xlim_min, xlim_max)
    ax5.set_ylim(ylim_min_L128, ylim_max_L128)
    ax5.set_yticks(np.linspace(ax5.get_ybound()[0], ax5.get_ybound()[1], 5))
    ax5.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    ax5.tick_params(axis='x', labelsize=tick_fontsize)
    ax5.tick_params(axis='y', labelsize=tick_fontsize)
    ax5.yaxis.get_offset_text().set_fontsize(tick_fontsize)
    ax5.grid()
    ax5.axhline(y=0, color='black', linestyle='-', linewidth=1.5)

    """
    6. L256 QV
    """

    ax6.plot(times, qv_wet_timeseries_L256[0], color="red", linestyle="solid", label="Moist")
    ax6.plot(times, qv_wet_timeseries_L256[1], color="red", linestyle="dashed")
    ax6.plot(times, qv_dry_timeseries_L256[0], color="blue", linestyle="solid", label="Dry")
    ax6.plot(times, qv_dry_timeseries_L256[1], color="blue", linestyle="dashed")

    ax6.set_xlim(xlim_min, xlim_max)
    ax6.set_ylim(ylim_min_L256, ylim_max_L256)
    ax6.set_yticks(np.linspace(ax6.get_ybound()[0], ax6.get_ybound()[1], 5))
    ax6.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    ax6.tick_params(axis='x', labelsize=tick_fontsize)
    ax6.tick_params(axis='y', labelsize=tick_fontsize)
    ax6.yaxis.get_offset_text().set_fontsize(tick_fontsize)
    ax6.grid()
    ax6.axhline(y=0, color='black', linestyle='-', linewidth=1.5)

    # -------------

    ax1.set_title(r"(a) L128, ${\theta}^{'}_{v}$/${\overline{\theta}_{v}}$", fontsize=title_fontsize)
    ax3.set_title(r"(b) L128, ${\theta}^{'}$/${\overline{\theta}}$", fontsize=title_fontsize)
    ax5.set_title(r"(c) L128, ${\epsilon}q_{v}^{'}$ / (${1 + {\epsilon}\overline{q_{v}}}$)", fontsize=title_fontsize)

    ax2.set_title(r"(d) L256, ${\theta}^{'}_{v}$/${\overline{\theta}_{v}}$", fontsize=title_fontsize)
    ax4.set_title(r"(e) L256, ${\theta}^{'}$/${\overline{\theta}}$", fontsize=title_fontsize)
    ax6.set_title(r"(f) L256, ${\epsilon}q_{v}^{'}$ / (${1 + {\epsilon}\overline{q_{v}}}$)", fontsize=title_fontsize)

    ax5.set_xlabel("Time [day]", fontsize=label_fontsize)
    ax6.set_xlabel("Time [day]", fontsize=label_fontsize)

    if to_plot == "agg":
        plt.suptitle("Aggregated",fontsize=suptitle_fontsize)
    if to_plot == "nonagg":
        plt.suptitle("Non-aggregated", fontsize=suptitle_fontsize)

    plt.tight_layout()
    plt.show()


def main():

    # SELECT AGGREGATION STATE TO PLOT (AGG or NON-AGG)
    to_plot = "agg" #agg or nonagg

    if to_plot == "agg":
        configs_L128 = ["0.001_L128", "0.002_L128"]
        configs_L256 = ["0.01_L256", "0.02_L256"]
    if to_plot == "nonagg":
        configs_L128 = ["0.01_L128", "0.02_L128"]
        configs_L256 = ["0.1_L256", "0.2_L256"]


    thetav_wet_timeseries_L128 = get_timeseries(configs_L128,"thetav_wet")
    thetav_dry_timeseries_L128 = get_timeseries(configs_L128,"thetav_dry")
    theta_wet_timeseries_L128 = get_timeseries(configs_L128,"theta_wet")
    theta_dry_timeseries_L128 = get_timeseries(configs_L128,"theta_dry")
    qv_wet_timeseries_L128 = get_timeseries(configs_L128,"qv_wet")
    qv_dry_timeseries_L128 = get_timeseries(configs_L128,"qv_dry")

    thetav_wet_timeseries_L256 = get_timeseries(configs_L256,"thetav_wet")
    thetav_dry_timeseries_L256 = get_timeseries(configs_L256,"thetav_dry")
    theta_wet_timeseries_L256 = get_timeseries(configs_L256,"theta_wet")
    theta_dry_timeseries_L256 = get_timeseries(configs_L256,"theta_dry")
    qv_wet_timeseries_L256 = get_timeseries(configs_L256,"qv_wet")
    qv_dry_timeseries_L256 = get_timeseries(configs_L256,"qv_dry")

    plot_timeseries(to_plot,thetav_wet_timeseries_L128,thetav_dry_timeseries_L128,theta_wet_timeseries_L128,theta_dry_timeseries_L128, qv_wet_timeseries_L128,qv_dry_timeseries_L128,
                             thetav_wet_timeseries_L256, thetav_dry_timeseries_L256,theta_wet_timeseries_L256, theta_dry_timeseries_L256,
                             qv_wet_timeseries_L256, qv_dry_timeseries_L256)


if __name__ == "__main__":

    main()