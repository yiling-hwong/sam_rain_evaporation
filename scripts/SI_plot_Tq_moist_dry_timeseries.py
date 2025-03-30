"""
@author: Yi-Ling HWONG
"""
import matplotlib as mpl
#mpl.rcParams['figure.dpi'] = 300
import matplotlib.pyplot as plt
from applications.utils import *

def get_timeseries(configs,target):

    timeseries_list = []
    for config in configs:
        file_path = "../data/Tq_moist_dry_timeseries/alpha"+config+"_"+target+".csv"
        timeseries = read_csv_to_list(file_path)

        if target == "qv_dom_mean":
            timeseries = [x*1000 for x in timeseries]

        timeseries_list.append(timeseries)

    return timeseries_list

def plot_timeseries(T_wet_timeseries_agg,T_dry_timeseries_agg,T_mean_timeseries_agg,
                    qv_wet_timeseries_agg,qv_dry_timeseries_agg,qv_mean_timeseries_agg,
                    T_wet_timeseries_nonagg, T_dry_timeseries_nonagg, T_mean_timeseries_nonagg,
                    qv_wet_timeseries_nonagg, qv_dry_timeseries_nonagg, qv_mean_timeseries_nonagg):

    num_row = 4
    num_col = 2
    bigframe_pad = 20

    xlim_min = 0  # in days
    xlim_max = 39
    title_fontsize = 16
    label_fontsize = 15
    tick_fontsize = 12
    legend_fontsize = 13


    # --------- Get ylim min and max
    ylim_min_T_agg = 294.5
    ylim_max_T_agg = 296
    ylim_min_T_nonagg = 294
    ylim_max_T_nonagg = 295

    ylim_min_qv_agg = 7
    ylim_max_qv_agg = 15
    ylim_min_qv_nonagg = 12
    ylim_max_qv_nonagg = 15

    # --------- GET TIME AXIS AND ZERO ARRAY
    times = [int(n) for n in range(len(T_wet_timeseries_agg[0]))]

    #------------ TITLE
    fig, big_axes = plt.subplots(figsize=(12, 10), nrows=4, ncols=1, sharex=True)

    for row, big_ax in enumerate(big_axes, start=1):

        if row == 1:
            big_ax.set_title("L128", fontsize=title_fontsize+4, weight="bold",pad=bigframe_pad)
        if row == 3:
            big_ax.set_title("L256", fontsize=title_fontsize+4, weight="bold", pad=bigframe_pad)

        # Turn off axis lines and ticks of the big subplot
        # obs alpha is 0 in RGBA string!
        big_ax.set_xticks([])
        big_ax.set_yticks([])
        # removes the white frame
        big_ax._frameon = False

    """
    1. L128_ALPHA0 T (AGG)
    """
    ax1 = fig.add_subplot(num_row, num_col, 1)
    ax1.plot(times, T_wet_timeseries_agg[0], color="red", linestyle="solid", label="Moist")
    ax1.plot(times, T_dry_timeseries_agg[0], color="blue", linestyle="solid",label="Dry")
    ax1.plot(times, T_mean_timeseries_agg[0], color="black", linestyle="dashed",label="Dom. mean")

    ax1.set_xlim(xlim_min, xlim_max)
    #ax1.set_ylim(ylim_min_T_agg, ylim_max_T_agg)
    ax1.set_yticks(np.linspace(ax1.get_ybound()[0], ax1.get_ybound()[1], 5))
    ax1.tick_params(axis='x', labelsize=tick_fontsize)
    ax1.tick_params(axis='y', labelsize=tick_fontsize)
    ax1.xaxis.set_tick_params(labelbottom=False)
    ax1.yaxis.get_offset_text().set_fontsize(tick_fontsize)
    ax1.set_ylabel(r"$\overline{T}_{PBL}$ [K]", fontsize=label_fontsize)
    ax1.set_title(r"(a) ${\alpha}$=0, agg", fontsize=title_fontsize)

    ax1.grid()

    """
    2. L128_ALPHA1 T (NONAGG)
    """
    ax2 = fig.add_subplot(num_row, num_col, 2)
    ax2.plot(times, T_wet_timeseries_nonagg[0], color="red", linestyle="solid", label="Moist")
    ax2.plot(times, T_dry_timeseries_nonagg[0], color="blue", linestyle="solid",label="Dry")
    ax2.plot(times, T_mean_timeseries_nonagg[0], color="black", linestyle="dashed",label="Dom. mean")

    ax2.set_xlim(xlim_min, xlim_max)
    #ax2.set_ylim(ylim_min_T_nonagg, ylim_max_T_nonagg)
    ax2.set_yticks(np.linspace(ax2.get_ybound()[0], ax2.get_ybound()[1], 5))
    ax2.tick_params(axis='x', labelsize=tick_fontsize)
    ax2.tick_params(axis='y', labelsize=tick_fontsize)
    ax2.xaxis.set_tick_params(labelbottom=False)
    ax2.yaxis.get_offset_text().set_fontsize(tick_fontsize)
    ax2.set_ylabel(r"$\overline{T}_{PBL}$ [K]", fontsize=label_fontsize)
    ax2.grid()
    ax2.legend(loc="lower center", fontsize=legend_fontsize, ncol=3)
    ax2.set_title(r"(e) ${\alpha}$=1, non-agg", fontsize=title_fontsize)

    """
    3. L128_ALPHA1 QV (AGG)
    """
    ax3 = fig.add_subplot(num_row, num_col, 3)
    ax3.plot(times, qv_wet_timeseries_agg[0], color="red", linestyle="solid", label="Moist")
    ax3.plot(times, qv_dry_timeseries_agg[0], color="blue", linestyle="solid",label="Dry")
    ax3.plot(times, qv_mean_timeseries_agg[0],color="black", linestyle="dashed",label="Dom. mean")

    ax3.set_xlim(xlim_min, xlim_max)
    #ax3.set_ylim(ylim_min_qv_agg, ylim_max_qv_agg)
    ax3.set_yticks(np.linspace(ax3.get_ybound()[0], ax3.get_ybound()[1], 5))
    ax3.tick_params(axis='x', labelsize=tick_fontsize)
    ax3.tick_params(axis='y', labelsize=tick_fontsize)
    ax3.xaxis.set_tick_params(labelbottom=False)
    ax3.yaxis.get_offset_text().set_fontsize(tick_fontsize)
    ax3.set_ylabel(r"$\overline{q}_{v,PBL}$ [g kg$^{-1}$]", fontsize=label_fontsize)
    ax3.set_title(r"(b) ${\alpha}$=0, agg", fontsize=title_fontsize)
    ax3.grid()

    """
    4. L128_ALPHA1 QV (NONAGG)
    """
    ax4 = fig.add_subplot(num_row, num_col, 4)
    ax4.plot(times, qv_wet_timeseries_nonagg[0], color="red", linestyle="solid", label="Moist")
    ax4.plot(times, qv_dry_timeseries_nonagg[0], color="blue", linestyle="solid",label="Dry")
    ax4.plot(times, qv_mean_timeseries_nonagg[0],color="black", linestyle="dashed",label="Dom. mean")

    ax4.set_xlim(xlim_min, xlim_max)
    #ax4.set_ylim(ylim_min_qv_nonagg, ylim_max_qv_nonagg)
    ax4.set_yticks(np.linspace(ax4.get_ybound()[0], ax4.get_ybound()[1], 5))
    ax4.tick_params(axis='x', labelsize=tick_fontsize)
    ax4.tick_params(axis='y', labelsize=tick_fontsize)
    ax4.xaxis.set_tick_params(labelbottom=False)
    ax4.yaxis.get_offset_text().set_fontsize(tick_fontsize)
    ax4.set_ylabel(r"$\overline{q}_{v,PBL}$ [g kg$^{-1}$]", fontsize=label_fontsize)
    ax4.set_title(r"(f) ${\alpha}$=1, non-agg", fontsize=title_fontsize)
    ax4.grid()

#######################
# L256
#######################

    """
    5. L256_ALPHA0 T (AGG)
    """
    ax5 = fig.add_subplot(num_row, num_col, 5)
    ax5.plot(times, T_wet_timeseries_agg[1], color="red", linestyle="solid", label="Moist")
    ax5.plot(times, T_dry_timeseries_agg[1], color="blue", linestyle="solid", label="Dry")
    ax5.plot(times, T_mean_timeseries_agg[1], color="black", linestyle="dashed", label="Dom. mean")

    ax5.set_xlim(xlim_min, xlim_max)
    #ax5.set_ylim(ylim_min_T_agg, ylim_max_T_agg)
    ax5.set_yticks(np.linspace(ax5.get_ybound()[0], ax5.get_ybound()[1], 5))
    ax5.tick_params(axis='x', labelsize=tick_fontsize)
    ax5.tick_params(axis='y', labelsize=tick_fontsize)
    ax5.xaxis.set_tick_params(labelbottom=False)
    ax5.yaxis.get_offset_text().set_fontsize(tick_fontsize)
    ax5.set_ylabel(r"$\overline{T}_{PBL}$ [K]", fontsize=label_fontsize)
    ax5.set_title(r"(c) ${\alpha}$=0, agg", fontsize=title_fontsize)
    ax5.grid()

    """
    6. L256_ALPHA1 T (NONAGG)
    """
    ax6 = fig.add_subplot(num_row, num_col, 6)
    ax6.plot(times, T_wet_timeseries_nonagg[1], color="red", linestyle="solid", label="Moist")
    ax6.plot(times, T_dry_timeseries_nonagg[1], color="blue", linestyle="solid", label="Dry")
    ax6.plot(times, T_mean_timeseries_nonagg[1], color="black", linestyle="dashed", label="Dom. mean")

    ax6.set_xlim(xlim_min, xlim_max)
    #ax6.set_ylim(ylim_min_T_nonagg, ylim_max_T_nonagg)
    ax6.set_yticks(np.linspace(ax6.get_ybound()[0], ax6.get_ybound()[1], 5))
    ax6.tick_params(axis='x', labelsize=tick_fontsize)
    ax6.tick_params(axis='y', labelsize=tick_fontsize)
    ax6.xaxis.set_tick_params(labelbottom=False)
    ax6.yaxis.get_offset_text().set_fontsize(tick_fontsize)
    ax6.set_ylabel(r"$\overline{T}_{PBL}$ [K]", fontsize=label_fontsize)
    ax6.set_title(r"(g) ${\alpha}$=1, non-agg", fontsize=title_fontsize)
    ax6.grid()

    """
    7. L256_ALPHA1 QV (AGG)
    """
    ax7 = fig.add_subplot(num_row, num_col, 7)
    ax7.plot(times, qv_wet_timeseries_agg[1], color="red", linestyle="solid", label="Moist")
    ax7.plot(times, qv_dry_timeseries_agg[1], color="blue", linestyle="solid", label="Dry")
    ax7.plot(times, qv_mean_timeseries_agg[1], color="black", linestyle="dashed", label="Dom. mean")

    ax7.set_xlim(xlim_min, xlim_max)
    #ax7.set_ylim(ylim_min_qv_agg, ylim_max_qv_agg)
    ax7.set_yticks(np.linspace(ax7.get_ybound()[0], ax7.get_ybound()[1], 5))
    ax7.tick_params(axis='x', labelsize=tick_fontsize)
    ax7.tick_params(axis='y', labelsize=tick_fontsize)
    ax7.yaxis.get_offset_text().set_fontsize(tick_fontsize)
    ax7.set_ylabel(r"$\overline{q}_{v,PBL}$ [g kg$^{-1}$]", fontsize=label_fontsize)
    ax7.set_title(r"(d) ${\alpha}$=0, agg", fontsize=title_fontsize)
    ax7.grid()

    """
    8. L256_ALPHA1 QV (NONAGG)
    """
    ax8 = fig.add_subplot(num_row, num_col, 8)
    ax8.plot(times, qv_wet_timeseries_nonagg[1], color="red", linestyle="solid", label="Moist")
    ax8.plot(times, qv_dry_timeseries_nonagg[1], color="blue", linestyle="solid", label="Dry")
    ax8.plot(times, qv_mean_timeseries_nonagg[1], color="black", linestyle="dashed", label="Dom. mean")

    ax8.set_xlim(xlim_min, xlim_max)
    #ax8.set_ylim(ylim_min_qv_nonagg, ylim_max_qv_nonagg)
    ax8.set_yticks(np.linspace(ax8.get_ybound()[0], ax8.get_ybound()[1], 5))
    ax8.tick_params(axis='x', labelsize=tick_fontsize)
    ax8.tick_params(axis='y', labelsize=tick_fontsize)
    ax8.set_ylabel(r"$\overline{q}_{v,PBL}$ [g kg$^{-1}$]", fontsize=label_fontsize)
    ax8.set_title(r"(h) ${\alpha}$=1, non-agg", fontsize=title_fontsize)
    ax8.yaxis.get_offset_text().set_fontsize(tick_fontsize)
    ax8.grid()

    ax7.set_xlabel("Time [day]", fontsize=label_fontsize)
    ax8.set_xlabel("Time [day]", fontsize=label_fontsize)

    #-----------------------------------------
    fig.subplots_adjust(left=0.1,
                        bottom=0.07,
                        right=0.98,
                        top=0.93,
                        wspace=0.28,  # 0.02
                        hspace=0.55)  # 0.65

    """
    wspace and hspace specify the space reserved between Matplotlib subplots. They are the fractions of axis width and height, respectively.
    left, right, top and bottom parameters specify four sides of the subplots’ positions. They are the fractions of the width and height of the figure.
    top and bottom should add up to 1.0
    """

    plt.show()


def main():

    configs_agg = ["0_L128", "0_L256"] # aggregated
    configs_nonagg = ["1_L128", "1_L256"]  # non-aggregated

    T_wet_timeseries_agg = get_timeseries(configs_agg,"T_wet")
    T_dry_timeseries_agg = get_timeseries(configs_agg, "T_dry")
    T_mean_timeseries_agg = get_timeseries(configs_agg, "T_dom_mean")
    qv_wet_timeseries_agg = get_timeseries(configs_agg,"qv_wet")
    qv_dry_timeseries_agg = get_timeseries(configs_agg, "qv_dry")
    qv_mean_timeseries_agg = get_timeseries(configs_agg, "qv_dom_mean")

    T_wet_timeseries_nonagg = get_timeseries(configs_nonagg,"T_wet")
    T_dry_timeseries_nonagg = get_timeseries(configs_nonagg, "T_dry")
    T_mean_timeseries_nonagg = get_timeseries(configs_nonagg, "T_dom_mean")
    qv_wet_timeseries_nonagg = get_timeseries(configs_nonagg,"qv_wet")
    qv_dry_timeseries_nonagg = get_timeseries(configs_nonagg, "qv_dry")
    qv_mean_timeseries_nonagg = get_timeseries(configs_nonagg, "qv_dom_mean")


    plot_timeseries(T_wet_timeseries_agg,T_dry_timeseries_agg,T_mean_timeseries_agg,
                    qv_wet_timeseries_agg,qv_dry_timeseries_agg,qv_mean_timeseries_agg,
                    T_wet_timeseries_nonagg, T_dry_timeseries_nonagg, T_mean_timeseries_nonagg,
                    qv_wet_timeseries_nonagg, qv_dry_timeseries_nonagg, qv_mean_timeseries_nonagg)


if __name__ == "__main__":

    main()