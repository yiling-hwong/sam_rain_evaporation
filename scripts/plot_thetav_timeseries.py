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

def plot_timeseries(thetav_wet_timeseries_agg,thetav_dry_timeseries_agg,theta_wet_timeseries_agg,theta_dry_timeseries_agg, qv_wet_timeseries_agg,qv_dry_timeseries_agg,
                             thetav_wet_timeseries_nonagg, thetav_dry_timeseries_nonagg,theta_wet_timeseries_nonagg, theta_dry_timeseries_nonagg,
                             qv_wet_timeseries_nonagg, qv_dry_timeseries_nonagg):

    num_row = 3
    num_col = 2
    bigframe_pad = 35

    xlim_min = 0  # in days
    xlim_max = 39
    title_fontsize = 16
    label_fontsize = 15
    tick_fontsize = 12
    legend_fontsize = 13
    sublegend_fontsize = 13

    # --------- Get ylim min and max
    ylim_min_agg = -2e-3
    ylim_max_agg = -ylim_min_agg
    ylim_min_nonagg = -2e-3 #-5e04
    ylim_max_nonagg = -ylim_min_nonagg

    # --------- GET TIME AXIS AND ZERO ARRAY
    times = [int(n) for n in range(len(thetav_wet_timeseries_agg[0]))]
    zero_arr = np.zeros(len(times))

    fig, big_axes = plt.subplots(figsize=(10, 8), nrows=1, ncols=2, sharex=True)

    for col, big_ax in enumerate(big_axes, start=0):

        if col == 0:
            big_ax.set_title("Aggregated", fontsize=title_fontsize+1, weight="bold",pad=bigframe_pad)
        if col == 1:
            big_ax.set_title("Non-aggregated", fontsize=title_fontsize+1, weight="bold", pad=bigframe_pad)

        # Turn off axis lines and ticks of the big subplot
        # obs alpha is 0 in RGBA string!
        big_ax.set_xticks([])
        big_ax.set_yticks([])
        # removes the white frame
        big_ax._frameon = False

    """
    1. AGG THETA_V
    """
    ax1 = fig.add_subplot(num_row, num_col, 1)
    ax1.plot(times, thetav_wet_timeseries_agg[0], color="red", linestyle="solid", label="Moist")
    ax1.plot(times, thetav_wet_timeseries_agg[1], color="red", linestyle="dashed")
    #ax1.plot(times, thetav_wet_timeseries_agg[2], color="red", linestyle="dotted")
    ax1.plot(times, thetav_dry_timeseries_agg[0], color="blue", linestyle="solid", label="Dry")
    ax1.plot(times, thetav_dry_timeseries_agg[1], color="blue", linestyle="dashed")
    #ax1.plot(times, thetav_dry_timeseries_agg[2], color="blue", linestyle="dotted")

    ax1.set_xlim(xlim_min, xlim_max)
    ax1.set_ylim(ylim_min_agg, ylim_max_agg)
    ax1.set_yticks(np.linspace(ax1.get_ybound()[0], ax1.get_ybound()[1], 5))
    ax1.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    ax1.tick_params(axis='x', labelsize=tick_fontsize)
    ax1.tick_params(axis='y', labelsize=tick_fontsize)
    ax1.yaxis.get_offset_text().set_fontsize(tick_fontsize)
    ax1.xaxis.set_tick_params(labelbottom=False)
    ax1.grid()
    ax1.legend(loc="upper left", fontsize=legend_fontsize)
    ax1.plot(times, zero_arr, color="black", linestyle="solid", linewidth=1.5)
    ax1.axhline(y=0, color='black', linestyle='-', linewidth=1.5)

    # ----------
    # Additional legend

    # custom_lines = [Line2D([0], [0], color="black", lw=1.5, linestyle="solid"),
    #                 Line2D([0], [0], color="black", lw=1.5, linestyle="dashed"),
    #                 Line2D([0], [0], color="black", lw=1.5, linestyle="dotted")]
    #
    # labels_leg2 = [r"L256_$\alpha$0", r"L128_$\alpha$0", r"L256_$\alpha$1_IntRad"]
    custom_lines = [Line2D([0], [0], color="black", lw=1.5, linestyle="dashed"),
                    Line2D([0], [0], color="black", lw=1.5, linestyle="solid")]

    labels_leg2 = ["L128", "L256"]
    leg = Legend(ax1, custom_lines, labels_leg2, loc=3, fontsize=sublegend_fontsize, ncol=2)
    ax1.add_artist(leg)

    # -----------

    """
    2. NON-AGG THETA_V
    """
    ax2 = fig.add_subplot(num_row, num_col, 2)
    ax2.plot(times, thetav_wet_timeseries_nonagg[0], color="red", linestyle="solid", label="Moist")
    ax2.plot(times, thetav_wet_timeseries_nonagg[1], color="red", linestyle="dashed")
    #ax2.plot(times, thetav_wet_timeseries_nonagg[2], color="red", linestyle="dotted")
    ax2.plot(times, thetav_dry_timeseries_nonagg[0], color="blue", linestyle="solid", label="Dry")
    ax2.plot(times, thetav_dry_timeseries_nonagg[1], color="blue", linestyle="dashed")
    #ax2.plot(times, thetav_dry_timeseries_nonagg[2], color="blue", linestyle="dotted")

    ax2.set_xlim(xlim_min, xlim_max)
    ax2.set_ylim(ylim_min_nonagg, ylim_max_nonagg)
    ax2.set_yticks(np.linspace(ax2.get_ybound()[0], ax2.get_ybound()[1], 5))
    ax2.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    ax2.tick_params(axis='x', labelsize=tick_fontsize)
    ax2.tick_params(axis='y', labelsize=tick_fontsize)
    ax2.yaxis.get_offset_text().set_fontsize(tick_fontsize)
    ax2.xaxis.set_tick_params(labelbottom=False)
    ax2.grid()
    ax2.legend(loc="upper left", fontsize=legend_fontsize)
    ax2.plot(times, zero_arr, color="black", linestyle="solid", linewidth=1.5)
    ax2.axhline(y=0, color='black', linestyle='-', linewidth=1.5)

    # ----------
    # Additional legend

    # custom_lines = [Line2D([0], [0], color="black", lw=1.5, linestyle="solid"),
    #                 Line2D([0], [0], color="black", lw=1.5, linestyle="dashed"),
    #                 Line2D([0], [0], color="black", lw=1.5, linestyle="dotted")]
    #
    # #labels_leg2 = [r"L256_$\alpha$0.1", r"L128_$\alpha$0.01",r"L128_$\alpha$1_IntRad"]
    # labels_leg2 = [r"L256_$\alpha$1", r"L128_$\alpha$1", r"L128_$\alpha$1_IntRad"]
    # leg = Legend(ax2, custom_lines, labels_leg2, loc=3, fontsize=sublegend_fontsize, ncol=3)
    # ax2.add_artist(leg)

    custom_lines = [Line2D([0], [0], color="black", lw=1.5, linestyle="dashed"),
                    Line2D([0], [0], color="black", lw=1.5, linestyle="solid")]

    labels_leg2 = ["L128", "L256"]
    leg = Legend(ax2, custom_lines, labels_leg2, loc="lower center", fontsize=sublegend_fontsize, ncol=2)
    ax2.add_artist(leg)

    # -----------

    """
    3. AGG THETA
    """
    ax3 = fig.add_subplot(num_row, num_col, 3)
    ax3.plot(times, theta_wet_timeseries_agg[0], color="red", linestyle="solid", label="Moist")
    ax3.plot(times, theta_wet_timeseries_agg[1], color="red", linestyle="dashed")
    #ax3.plot(times, theta_wet_timeseries_agg[2], color="red", linestyle="dotted")
    ax3.plot(times, theta_dry_timeseries_agg[0], color="blue", linestyle="solid", label="Dry")
    ax3.plot(times, theta_dry_timeseries_agg[1], color="blue", linestyle="dashed")
    #ax3.plot(times, theta_dry_timeseries_agg[2], color="blue", linestyle="dotted")

    ax3.set_xlim(xlim_min, xlim_max)
    ax3.set_ylim(ylim_min_agg, ylim_max_agg)
    ax3.set_yticks(np.linspace(ax3.get_ybound()[0], ax3.get_ybound()[1], 5))
    ax3.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    ax3.tick_params(axis='x', labelsize=tick_fontsize)
    ax3.tick_params(axis='y', labelsize=tick_fontsize)
    ax3.yaxis.get_offset_text().set_fontsize(tick_fontsize)
    ax3.xaxis.set_tick_params(labelbottom=False)
    ax3.grid()
    ax3.axhline(y=0, color='black', linestyle='-', linewidth=1.5)

    """
    4. NON-AGG THETA
    """
    ax4 = fig.add_subplot(num_row, num_col, 4)
    ax4.plot(times, theta_wet_timeseries_nonagg[0], color="red", linestyle="solid", label="Moist")
    ax4.plot(times, theta_wet_timeseries_nonagg[1], color="red", linestyle="dashed")
    #ax4.plot(times, theta_wet_timeseries_nonagg[2], color="red", linestyle="dotted")
    ax4.plot(times, theta_dry_timeseries_nonagg[0], color="blue", linestyle="solid", label="Dry")
    ax4.plot(times, theta_dry_timeseries_nonagg[1], color="blue", linestyle="dashed")
    #ax4.plot(times, theta_dry_timeseries_nonagg[2], color="blue", linestyle="dotted")

    ax4.set_xlim(xlim_min, xlim_max)
    ax4.set_ylim(ylim_min_nonagg, ylim_max_nonagg)
    ax4.set_yticks(np.linspace(ax4.get_ybound()[0], ax4.get_ybound()[1], 5))
    ax4.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    ax4.tick_params(axis='x', labelsize=tick_fontsize)
    ax4.tick_params(axis='y', labelsize=tick_fontsize)
    ax4.yaxis.get_offset_text().set_fontsize(tick_fontsize)
    ax4.xaxis.set_tick_params(labelbottom=False)
    ax4.grid()
    ax4.axhline(y=0, color='black', linestyle='-', linewidth=1.5)

    """
    5. AGG QV
    """
    ax5 = fig.add_subplot(num_row, num_col, 5)
    ax5.plot(times, qv_wet_timeseries_agg[0], color="red", linestyle="solid", label="Moist")
    ax5.plot(times, qv_wet_timeseries_agg[1], color="red", linestyle="dashed")
    #ax5.plot(times, qv_wet_timeseries_agg[2], color="red", linestyle="dotted")
    ax5.plot(times, qv_dry_timeseries_agg[0], color="blue", linestyle="solid", label="Dry")
    ax5.plot(times, qv_dry_timeseries_agg[1], color="blue", linestyle="dashed")
    #ax5.plot(times, qv_dry_timeseries_agg[2], color="blue", linestyle="dotted")

    ax5.set_xlim(xlim_min, xlim_max)
    ax5.set_ylim(ylim_min_agg, ylim_max_agg)
    ax5.set_yticks(np.linspace(ax5.get_ybound()[0], ax5.get_ybound()[1], 5))
    ax5.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    ax5.tick_params(axis='x', labelsize=tick_fontsize)
    ax5.tick_params(axis='y', labelsize=tick_fontsize)
    ax5.yaxis.get_offset_text().set_fontsize(tick_fontsize)
    ax5.grid()
    ax5.axhline(y=0, color='black', linestyle='-', linewidth=1.5)

    """
    6. NON-AGG QV
    """
    ax6 = fig.add_subplot(num_row, num_col, 6)
    ax6.plot(times, qv_wet_timeseries_nonagg[0], color="red", linestyle="solid", label="Moist")
    ax6.plot(times, qv_wet_timeseries_nonagg[1], color="red", linestyle="dashed")
    #ax6.plot(times, qv_wet_timeseries_nonagg[2], color="red", linestyle="dotted")
    ax6.plot(times, qv_dry_timeseries_nonagg[0], color="blue", linestyle="solid", label="Dry")
    ax6.plot(times, qv_dry_timeseries_nonagg[1], color="blue", linestyle="dashed")
    #ax6.plot(times, qv_dry_timeseries_nonagg[2], color="blue", linestyle="dotted")

    ax6.set_xlim(xlim_min, xlim_max)
    ax6.set_ylim(ylim_min_nonagg, ylim_max_nonagg)
    ax6.set_yticks(np.linspace(ax6.get_ybound()[0], ax6.get_ybound()[1], 5))
    ax6.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    ax6.tick_params(axis='x', labelsize=tick_fontsize)
    ax6.tick_params(axis='y', labelsize=tick_fontsize)
    ax6.yaxis.get_offset_text().set_fontsize(tick_fontsize)
    ax6.grid()
    ax6.axhline(y=0, color='black', linestyle='-', linewidth=1.5)

    # -------------

    ax1.set_title(r"(a) ${\theta}^{'}_{v}$/${\overline{\theta}_{v}}, \alpha = 0$", fontsize=title_fontsize)
    ax3.set_title(r"(b) ${\theta}^{'}$/${\overline{\theta}}, \alpha = 0$", fontsize=title_fontsize)
    ax5.set_title(r"(c) ${\epsilon}q_{v}^{'}$ / (${1 + {\epsilon}\overline{q_{v}}}$), $\alpha = 0$", fontsize=title_fontsize)

    ax2.set_title(r"(d) ${\theta}^{'}_{v}$/${\overline{\theta}_{v}}, \alpha = 1$", fontsize=title_fontsize)
    ax4.set_title(r"(e) ${\theta}^{'}$/${\overline{\theta}}, \alpha = 1$", fontsize=title_fontsize)
    ax6.set_title(r"(f) ${\epsilon}q_{v}^{'}$ / (${1 + {\epsilon}\overline{q_{v}}}$), $\alpha = 1$", fontsize=title_fontsize)

    ax5.set_xlabel("Time [day]", fontsize=label_fontsize)
    ax6.set_xlabel("Time [day]", fontsize=label_fontsize)

    #-----------------------------------------
    fig.subplots_adjust(left=0.05,
                        bottom=0.08,
                        right=0.98,
                        top=0.9,
                        wspace=0.1,
                        hspace=0.35)

    """
    wspace and hspace specify the space reserved between Matplotlib subplots. They are the fractions of axis width and height, respectively.
    left, right, top and bottom parameters specify four sides of the subplots’ positions. They are the fractions of the width and height of the figure.
    top and bottom should add up to 1.0
    """

    #plt.tight_layout()
    plt.show()


def main():

    configs_agg = ["0_L256", "0_L128"] # aggregated
    configs_nonagg = ["1_L256", "1_L128"]  # non-aggregated

    thetav_wet_timeseries_agg = get_timeseries(configs_agg,"thetav_wet")
    thetav_dry_timeseries_agg = get_timeseries(configs_agg,"thetav_dry")
    theta_wet_timeseries_agg = get_timeseries(configs_agg,"theta_wet")
    theta_dry_timeseries_agg = get_timeseries(configs_agg,"theta_dry")
    qv_wet_timeseries_agg = get_timeseries(configs_agg,"qv_wet")
    qv_dry_timeseries_agg = get_timeseries(configs_agg,"qv_dry")

    thetav_wet_timeseries_nonagg = get_timeseries(configs_nonagg,"thetav_wet")
    thetav_dry_timeseries_nonagg = get_timeseries(configs_nonagg,"thetav_dry")
    theta_wet_timeseries_nonagg = get_timeseries(configs_nonagg,"theta_wet")
    theta_dry_timeseries_nonagg = get_timeseries(configs_nonagg,"theta_dry")
    qv_wet_timeseries_nonagg = get_timeseries(configs_nonagg,"qv_wet")
    qv_dry_timeseries_nonagg = get_timeseries(configs_nonagg,"qv_dry")


    plot_timeseries(thetav_wet_timeseries_agg,thetav_dry_timeseries_agg,theta_wet_timeseries_agg,theta_dry_timeseries_agg, qv_wet_timeseries_agg,qv_dry_timeseries_agg,
                             thetav_wet_timeseries_nonagg, thetav_dry_timeseries_nonagg,theta_wet_timeseries_nonagg, theta_dry_timeseries_nonagg,
                             qv_wet_timeseries_nonagg, qv_dry_timeseries_nonagg)


if __name__ == "__main__":

    main()