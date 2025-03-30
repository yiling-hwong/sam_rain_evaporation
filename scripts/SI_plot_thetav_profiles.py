"""
@author: Yi-Ling HWONG
"""
import matplotlib as mpl
#mpl.rcParams['figure.dpi'] = 300
import matplotlib.pyplot as plt
from applications.utils import *

def get_profiles(config,time_list,target):

    profile_list = []

    for time in time_list:
        file_path = "../data/thetav_profile/" + config + "_" + target + "_day"+str(time)+".csv"
        profile = read_csv_to_list(file_path)

        profile_list.append(profile)

    return profile_list

def get_heights(config):

    file_path = "../data/thetav_profile/" + config + "_z.csv"
    heights = read_csv_to_list(file_path)

    print (len(heights))

    return heights

def plot_profiles(time_list_L128,z_list_L128,thetav_wet_profiles_L128,thetav_dry_profiles_L128,theta_wet_profiles_L128,theta_dry_profiles_L128,qv_wet_profiles_L128,qv_dry_profiles_L128,
                  time_list_L256,z_list_L256,thetav_wet_profiles_L256,thetav_dry_profiles_L256,theta_wet_profiles_L256,theta_dry_profiles_L256,qv_wet_profiles_L256,qv_dry_profiles_L256):

    rows = 4
    cols = 3
    bigframe_pad = 30

    ylim_min = 0
    ylim_max = 4  # 4; max 19
    xlim = 0.002
    xlim_min_thetav = -xlim
    xlim_max_thetav = xlim


    tick_fontsize = 13
    label_fontsize = 14
    title_fontsize = 15
    legend_fontsize = 14
    colors = get_line_colors("YlOrRd_r", len(time_list_L128))  # YlOrRd_r,Dark2,rainbow,jet,turbo,gnuplot,gnuplot2,tab10,Accent
    # https://matplotlib.org/stable/tutorials/colors/colormaps.html

    fig, big_axes = plt.subplots(figsize=(7, 14), nrows=4, ncols=1, sharey=True)

    for row, big_ax in enumerate(big_axes, start=1):

        if row == 1:
            big_ax.set_title("L128, Moist", fontsize=title_fontsize, pad=bigframe_pad)
        if row == 2:
            big_ax.set_title("L128, Dry", fontsize=title_fontsize, pad=bigframe_pad)
        if row == 3:
            big_ax.set_title("L256, Moist", fontsize=title_fontsize, pad=bigframe_pad)
        if row == 4:
            big_ax.set_title("L256, Dry", fontsize=title_fontsize, pad=bigframe_pad)

        # Turn off axis lines and ticks of the big subplot
        # obs alpha is 0 in RGBA string!
        big_ax.set_xticks([])
        big_ax.set_yticks([])
        # removes the white frame
        big_ax._frameon = False

    """
    L128 WET THETA_V
    """
    ax = fig.add_subplot(rows,cols,1)
    for index,time in enumerate(time_list_L128):
        to_plot = thetav_wet_profiles_L128[index]
        ax.plot(to_plot,z_list_L128,color=colors[index],label="day "+str(time))
        plt.axvline(x=0, color="black", linestyle="dashed")
        plt.axhline(y=1, color="black", linestyle="dashed")
        ax.legend(loc="best",fontsize=label_fontsize)
        ax.set_xlim([xlim_min_thetav,xlim_max_thetav])
        ax.set_ylim([ylim_min,ylim_max])
        plt.yticks(fontsize=tick_fontsize)
        plt.xticks(fontsize=tick_fontsize)
        plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
        ax.tick_params(axis='x', labelsize=tick_fontsize)
        ax.tick_params(axis='y', labelsize=tick_fontsize)
        ax.xaxis.get_offset_text().set_fontsize(tick_fontsize)
        plt.ylabel("Height [km]",fontsize=label_fontsize)
        plt.title("(a) " + r"${\theta}^{'}_{v}$/${\overline{\theta}_{v}}$",fontsize = title_fontsize)
        plt.grid()

        ax.get_legend().remove()

    """
    L128 WET THETA
    """
    ax = fig.add_subplot(rows,cols,2)
    for index,time in enumerate(time_list_L128):
        to_plot = theta_wet_profiles_L128[index]
        ax.plot(to_plot,z_list_L128,color=colors[index],label="day "+str(time))
        plt.axvline(x=0, color="black", linestyle="dashed")
        plt.axhline(y=1, color="black", linestyle="dashed")
        ax.legend(loc="best",fontsize=label_fontsize)
        ax.set_xlim([xlim_min_thetav,xlim_max_thetav])
        ax.set_ylim([ylim_min,ylim_max])
        plt.yticks(fontsize=tick_fontsize)
        plt.xticks(fontsize=tick_fontsize)
        plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
        ax.tick_params(axis='x', labelsize=tick_fontsize)
        ax.tick_params(axis='y', labelsize=tick_fontsize)
        ax.xaxis.get_offset_text().set_fontsize(tick_fontsize)
        plt.ylabel("Height [km]",fontsize=label_fontsize)
        plt.title("(b) " + r"${\theta}^{'}$/${\overline{\theta}}$",fontsize = title_fontsize)
        plt.grid()

        ax.get_legend().remove()
        ax.yaxis.set_ticklabels([])
        ax.set(ylabel=None)

    """
    L128 WET QV
    """
    ax = fig.add_subplot(rows,cols,3)
    for index,time in enumerate(time_list_L128):
        to_plot = qv_wet_profiles_L128[index]
        ax.plot(to_plot,z_list_L128,color=colors[index],label="day "+str(time))
        plt.axvline(x=0, color="black", linestyle="dashed")
        plt.axhline(y=1, color="black", linestyle="dashed")
        ax.legend(loc="best",fontsize=label_fontsize)
        ax.set_xlim([xlim_min_thetav,xlim_max_thetav])
        ax.set_ylim([ylim_min,ylim_max])
        plt.yticks(fontsize=tick_fontsize)
        plt.xticks(fontsize=tick_fontsize)
        plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
        ax.tick_params(axis='x', labelsize=tick_fontsize)
        ax.tick_params(axis='y', labelsize=tick_fontsize)
        ax.xaxis.get_offset_text().set_fontsize(tick_fontsize)
        plt.ylabel("Height [km]",fontsize=label_fontsize)
        plt.title("(c) " + r"${\epsilon}q_{v}^{'}$ / (${1 + {\epsilon}\overline{q_{v}}}$)",fontsize = title_fontsize)
        plt.grid()

        ax.get_legend().remove()
        ax.yaxis.set_ticklabels([])
        ax.set(ylabel=None)

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles=handles, labels=labels, bbox_to_anchor=(1.96, 0.52), loc="center right", ncol=1,fontsize=legend_fontsize)

    """
    L128 DRY THETA_V
    """
    ax = fig.add_subplot(rows,cols,4)
    for index,time in enumerate(time_list_L128):
        to_plot = thetav_dry_profiles_L128[index]
        ax.plot(to_plot,z_list_L128,color=colors[index],label="day "+str(time))
        plt.axvline(x=0, color="black", linestyle="dashed")
        plt.axhline(y=1, color="black", linestyle="dashed")
        ax.legend(loc="best",fontsize=label_fontsize)
        ax.set_xlim([xlim_min_thetav,xlim_max_thetav])
        ax.set_ylim([ylim_min,ylim_max])
        plt.yticks(fontsize=tick_fontsize)
        plt.xticks(fontsize=tick_fontsize)
        plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
        ax.tick_params(axis='x', labelsize=tick_fontsize)
        ax.tick_params(axis='y', labelsize=tick_fontsize)
        ax.xaxis.get_offset_text().set_fontsize(tick_fontsize)
        plt.ylabel("Height [km]",fontsize=label_fontsize)
        plt.title("(d) " + r"${\theta}^{'}_{v}$/${\overline{\theta}_{v}}$",fontsize = title_fontsize)
        plt.grid()

        ax.get_legend().remove()

    """
    L128 DRY THETA
    """
    ax = fig.add_subplot(rows,cols,5)
    for index,time in enumerate(time_list_L128):
        to_plot = theta_dry_profiles_L128[index]
        ax.plot(to_plot,z_list_L128,color=colors[index],label="day "+str(time))
        plt.axvline(x=0, color="black", linestyle="dashed")
        plt.axhline(y=1, color="black", linestyle="dashed")
        ax.legend(loc="best",fontsize=label_fontsize)
        ax.set_xlim([xlim_min_thetav,xlim_max_thetav])
        ax.set_ylim([ylim_min,ylim_max])
        plt.yticks(fontsize=tick_fontsize)
        plt.xticks(fontsize=tick_fontsize)
        plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
        ax.tick_params(axis='x', labelsize=tick_fontsize)
        ax.tick_params(axis='y', labelsize=tick_fontsize)
        ax.xaxis.get_offset_text().set_fontsize(tick_fontsize)
        plt.ylabel("Height [km]",fontsize=label_fontsize)
        plt.title("(e) " + r"${\theta}^{'}$/${\overline{\theta}}$",fontsize = title_fontsize)
        plt.grid()

        ax.get_legend().remove()
        ax.yaxis.set_ticklabels([])
        ax.set(ylabel=None)

    """
    L128 DRY QV
    """
    ax = fig.add_subplot(rows,cols,6)
    for index,time in enumerate(time_list_L128):
        to_plot = qv_dry_profiles_L128[index]
        ax.plot(to_plot,z_list_L128,color=colors[index],label="day "+str(time))
        plt.axvline(x=0, color="black", linestyle="dashed")
        plt.axhline(y=1, color="black", linestyle="dashed")
        ax.legend(loc="best",fontsize=label_fontsize)
        ax.set_xlim([xlim_min_thetav,xlim_max_thetav])
        ax.set_ylim([ylim_min,ylim_max])
        plt.yticks(fontsize=tick_fontsize)
        plt.xticks(fontsize=tick_fontsize)
        plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
        ax.tick_params(axis='x', labelsize=tick_fontsize)
        ax.tick_params(axis='y', labelsize=tick_fontsize)
        ax.xaxis.get_offset_text().set_fontsize(tick_fontsize)
        plt.ylabel("Height [km]",fontsize=label_fontsize)
        plt.title("(f) " + r"${\epsilon}q_{v}^{'}$ / (${1 + {\epsilon}\overline{q_{v}}}$)",fontsize = title_fontsize)
        plt.grid()

        ax.get_legend().remove()
        ax.yaxis.set_ticklabels([])
        ax.set(ylabel=None)

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles=handles, labels=labels, bbox_to_anchor=(1.96, 0.52), loc="center right", ncol=1,fontsize=legend_fontsize)


############################

    """
    L256 WET THETA_V
    """
    ax = fig.add_subplot(rows, cols, 7)
    for index, time in enumerate(time_list_L256):
        to_plot = thetav_wet_profiles_L256[index]
        ax.plot(to_plot, z_list_L256, color=colors[index], label="day " + str(time))
        plt.axvline(x=0, color="black", linestyle="dashed")
        plt.axhline(y=1, color="black", linestyle="dashed")
        ax.legend(loc="best", fontsize=label_fontsize)
        ax.set_xlim([xlim_min_thetav, xlim_max_thetav])
        ax.set_ylim([ylim_min, ylim_max])
        plt.yticks(fontsize=tick_fontsize)
        plt.xticks(fontsize=tick_fontsize)
        plt.ticklabel_format(axis="x", style="sci", scilimits=(0, 0))
        ax.tick_params(axis='x', labelsize=tick_fontsize)
        ax.tick_params(axis='y', labelsize=tick_fontsize)
        ax.xaxis.get_offset_text().set_fontsize(tick_fontsize)
        plt.ylabel("Height [km]", fontsize=label_fontsize)
        plt.title("(g) " + r"${\theta}^{'}_{v}$/${\overline{\theta}_{v}}$", fontsize=title_fontsize)
        plt.grid()

        ax.get_legend().remove()

    """
    L256 WET THETA
    """
    ax = fig.add_subplot(rows, cols, 8)
    for index, time in enumerate(time_list_L256):
        to_plot = theta_wet_profiles_L256[index]
        ax.plot(to_plot, z_list_L256, color=colors[index], label="day " + str(time))
        plt.axvline(x=0, color="black", linestyle="dashed")
        plt.axhline(y=1, color="black", linestyle="dashed")
        ax.legend(loc="best", fontsize=label_fontsize)
        ax.set_xlim([xlim_min_thetav, xlim_max_thetav])
        ax.set_ylim([ylim_min, ylim_max])
        plt.yticks(fontsize=tick_fontsize)
        plt.xticks(fontsize=tick_fontsize)
        plt.ticklabel_format(axis="x", style="sci", scilimits=(0, 0))
        ax.tick_params(axis='x', labelsize=tick_fontsize)
        ax.tick_params(axis='y', labelsize=tick_fontsize)
        ax.xaxis.get_offset_text().set_fontsize(tick_fontsize)
        plt.ylabel("Height [km]", fontsize=label_fontsize)
        plt.title("(h) " + r"${\theta}^{'}$/${\overline{\theta}}$", fontsize=title_fontsize)
        plt.grid()

        ax.get_legend().remove()
        ax.yaxis.set_ticklabels([])
        ax.set(ylabel=None)

    """
    L256 WET QV
    """
    ax = fig.add_subplot(rows, cols, 9)
    for index, time in enumerate(time_list_L256):
        to_plot = qv_wet_profiles_L256[index]
        ax.plot(to_plot, z_list_L256, color=colors[index], label="day " + str(time))
        plt.axvline(x=0, color="black", linestyle="dashed")
        plt.axhline(y=1, color="black", linestyle="dashed")
        ax.legend(loc="best", fontsize=label_fontsize)
        ax.set_xlim([xlim_min_thetav, xlim_max_thetav])
        ax.set_ylim([ylim_min, ylim_max])
        plt.yticks(fontsize=tick_fontsize)
        plt.xticks(fontsize=tick_fontsize)
        plt.ticklabel_format(axis="x", style="sci", scilimits=(0, 0))
        ax.tick_params(axis='x', labelsize=tick_fontsize)
        ax.tick_params(axis='y', labelsize=tick_fontsize)
        ax.xaxis.get_offset_text().set_fontsize(tick_fontsize)
        plt.ylabel("Height [km]", fontsize=label_fontsize)
        plt.title("(i) " + r"${\epsilon}q_{v}^{'}$ / (${1 + {\epsilon}\overline{q_{v}}}$)", fontsize=title_fontsize)
        plt.grid()

        ax.get_legend().remove()
        ax.yaxis.set_ticklabels([])
        ax.set(ylabel=None)

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles=handles, labels=labels, bbox_to_anchor=(2.04, 0.52), loc="center right", ncol=1,fontsize=legend_fontsize)


    """
    L256 DRY THETA_V
    """
    ax = fig.add_subplot(rows, cols, 10)
    for index, time in enumerate(time_list_L256):
        to_plot = thetav_dry_profiles_L256[index]
        ax.plot(to_plot, z_list_L256, color=colors[index], label="day " + str(time))
        plt.axvline(x=0, color="black", linestyle="dashed")
        plt.axhline(y=1, color="black", linestyle="dashed")
        ax.legend(loc="best", fontsize=label_fontsize)
        ax.set_xlim([xlim_min_thetav, xlim_max_thetav])
        ax.set_ylim([ylim_min, ylim_max])
        plt.yticks(fontsize=tick_fontsize)
        plt.xticks(fontsize=tick_fontsize)
        plt.ticklabel_format(axis="x", style="sci", scilimits=(0, 0))
        ax.tick_params(axis='x', labelsize=tick_fontsize)
        ax.tick_params(axis='y', labelsize=tick_fontsize)
        ax.xaxis.get_offset_text().set_fontsize(tick_fontsize)
        plt.ylabel("Height [km]", fontsize=label_fontsize)
        plt.title("(j) " + r"${\theta}^{'}_{v}$/${\overline{\theta}_{v}}$", fontsize=title_fontsize)
        plt.grid()

        ax.get_legend().remove()

    """
    L256 DRY THETA
    """
    ax = fig.add_subplot(rows, cols, 11)
    for index, time in enumerate(time_list_L256):
        to_plot = theta_dry_profiles_L256[index]
        ax.plot(to_plot, z_list_L256, color=colors[index], label="day " + str(time))
        plt.axvline(x=0, color="black", linestyle="dashed")
        plt.axhline(y=1, color="black", linestyle="dashed")
        ax.legend(loc="best", fontsize=label_fontsize)
        ax.set_xlim([xlim_min_thetav, xlim_max_thetav])
        ax.set_ylim([ylim_min, ylim_max])
        plt.yticks(fontsize=tick_fontsize)
        plt.xticks(fontsize=tick_fontsize)
        plt.ticklabel_format(axis="x", style="sci", scilimits=(0, 0))
        ax.tick_params(axis='x', labelsize=tick_fontsize)
        ax.tick_params(axis='y', labelsize=tick_fontsize)
        ax.xaxis.get_offset_text().set_fontsize(tick_fontsize)
        plt.ylabel("Height [km]", fontsize=label_fontsize)
        plt.title("(k) " + r"${\theta}^{'}$/${\overline{\theta}}$", fontsize=title_fontsize)
        plt.grid()

        ax.get_legend().remove()
        ax.yaxis.set_ticklabels([])
        ax.set(ylabel=None)

    """
    L256 DRY QV
    """
    ax = fig.add_subplot(rows, cols, 12)
    for index, time in enumerate(time_list_L256):
        to_plot = qv_dry_profiles_L256[index]
        ax.plot(to_plot, z_list_L256, color=colors[index], label="day " + str(time))
        plt.axvline(x=0, color="black", linestyle="dashed")
        plt.axhline(y=1, color="black", linestyle="dashed")
        ax.legend(loc="best", fontsize=label_fontsize)
        ax.set_xlim([xlim_min_thetav, xlim_max_thetav])
        ax.set_ylim([ylim_min, ylim_max])
        plt.yticks(fontsize=tick_fontsize)
        plt.xticks(fontsize=tick_fontsize)
        plt.ticklabel_format(axis="x", style="sci", scilimits=(0, 0))
        ax.tick_params(axis='x', labelsize=tick_fontsize)
        ax.tick_params(axis='y', labelsize=tick_fontsize)
        ax.xaxis.get_offset_text().set_fontsize(tick_fontsize)
        plt.ylabel("Height [km]", fontsize=label_fontsize)
        plt.title("(l) " + r"${\epsilon}q_{v}^{'}$ / (${1 + {\epsilon}\overline{q_{v}}}$)", fontsize=title_fontsize)
        plt.grid()

        ax.get_legend().remove()
        ax.yaxis.set_ticklabels([])
        ax.set(ylabel=None)

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles=handles, labels=labels, bbox_to_anchor=(2.04, 0.52), loc="center right", ncol=1,fontsize=legend_fontsize)


    #-----------------------------------------
    fig.subplots_adjust(left=0.08,
                        bottom=0.04,
                        right=0.77,
                        top=0.95,
                        wspace=0.12,  # 0.02
                        hspace=0.6)  # 0.65

    """
    wspace and hspace specify the space reserved between Matplotlib subplots. They are the fractions of axis width and height, respectively.
    left, right, top and bottom parameters specify four sides of the subplots’ positions. They are the fractions of the width and height of the figure.
    top and bottom should add up to 1.0
    """

    plt.show()


def main():

    config_L128 = "alpha0_L128"
    config_L256 = "alpha0_L256"
    time_list_L128 = [1,3,5,7,9]
    time_list_L256 = [3,5,7,9,12]

    z_list_L128 = get_heights(config_L128)
    thetav_wet_profiles_L128 = get_profiles(config_L128,time_list_L128,"thetav_wet")
    thetav_dry_profiles_L128 = get_profiles(config_L128,time_list_L128,"thetav_dry")
    theta_wet_profiles_L128 = get_profiles(config_L128,time_list_L128,"theta_wet")
    theta_dry_profiles_L128 = get_profiles(config_L128,time_list_L128,"theta_dry")
    qv_wet_profiles_L128 = get_profiles(config_L128,time_list_L128,"qv_wet")
    qv_dry_profiles_L128 = get_profiles(config_L128,time_list_L128,"qv_dry")

    z_list_L256 = get_heights(config_L256)
    thetav_wet_profiles_L256 = get_profiles(config_L256,time_list_L256,"thetav_wet")
    thetav_dry_profiles_L256 = get_profiles(config_L256,time_list_L256,"thetav_dry")
    theta_wet_profiles_L256 = get_profiles(config_L256,time_list_L256,"theta_wet")
    theta_dry_profiles_L256 = get_profiles(config_L256,time_list_L256,"theta_dry")
    qv_wet_profiles_L256 = get_profiles(config_L256,time_list_L256,"qv_wet")
    qv_dry_profiles_L256 = get_profiles(config_L256,time_list_L256,"qv_dry")

    plot_profiles(time_list_L128,z_list_L128,thetav_wet_profiles_L128,thetav_dry_profiles_L128,theta_wet_profiles_L128,theta_dry_profiles_L128,qv_wet_profiles_L128,qv_dry_profiles_L128,
                  time_list_L256,z_list_L256,thetav_wet_profiles_L256,thetav_dry_profiles_L256,theta_wet_profiles_L256,theta_dry_profiles_L256,qv_wet_profiles_L256,qv_dry_profiles_L256)

if __name__ == "__main__":

    main()