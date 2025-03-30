"""
@author: Yi-Ling HWONG
"""
import matplotlib as mpl
#mpl.rcParams['figure.dpi'] = 300
import matplotlib.pyplot as plt
from applications.utils import *


def get_array(config,var):

    file_path = "../data/"+var+"/"+config+".csv"
    arr = np.loadtxt(file_path,delimiter=",")

    return arr

def plot_figure(coldpools_all,pw_all,plot_labels):

    num_row = 4
    num_col = 4
    bigframe_pad = 10

    vmin_cp_L128 = -1
    vmax_cp_L128 = -vmin_cp_L128
    vmin_cp_L256 = -5
    vmax_cp_L256 = -vmin_cp_L256

    vmin_pw_L128 = 28
    vmax_pw_L128 = 48
    vmin_pw_L256 = 28
    vmax_pw_L256 = 48

    title_fontsize = 16
    label_fontsize = 14
    tick_fontsize = 13

    cmap_cp = "seismic"
    cmap_pw = "rainbow"

    fig, big_axes = plt.subplots(figsize=(10, 12), nrows=4, ncols=1, sharey=True)

    for row, big_ax in enumerate(big_axes, start=1):

        if row == 1:
            big_ax.set_title("L128, precipitable water", fontsize=title_fontsize, pad=bigframe_pad)
        if row == 2:
            big_ax.set_title("L128, cold pools", fontsize=title_fontsize, pad=bigframe_pad)
        if row == 3:
            big_ax.set_title("L256, precipitable water", fontsize=title_fontsize, pad=bigframe_pad)
        if row == 4:
            big_ax.set_title("L256, cold pools", fontsize=title_fontsize, pad=bigframe_pad)

        # Turn off axis lines and ticks of the big subplot
        # obs alpha is 0 in RGBA string!
        big_ax.set_xticks([])
        big_ax.set_yticks([])
        # removes the white frame
        big_ax._frameon = False

    """
    1. L128 PW
    """
    for n in range(0, 4):
        ax = fig.add_subplot(num_row, num_col, n + 1)
        im = ax.pcolormesh(pw_all[n], vmin=vmin_pw_L128, vmax=vmax_pw_L128, cmap=cmap_pw)
        plt.title(plot_labels[n], fontsize=label_fontsize)
        plt.gca().set_aspect('equal')

        if n == 0:
            ax.set_xlabel("x [km]", fontsize=label_fontsize)
            ax.set_ylabel("y [km]", fontsize=label_fontsize)
            ax.set_xticks([-0, 50, 100])
            ax.set_yticks([-0, 50, 100])
            ax.tick_params(axis='x', labelsize=tick_fontsize)
            ax.tick_params(axis='y', labelsize=tick_fontsize)
        else:
            ax.xaxis.set_tick_params(labelbottom=False)
            ax.yaxis.set_tick_params(labelleft=False)

    # COLORBAR
    cbar_ax = fig.add_axes([0.9, 0.78, 0.01, 0.15])
      # 1 rows, [left, bottom, width, height]
    cbar = fig.colorbar(im, cax=cbar_ax, extend="both")
    cbar.ax.set_ylabel(r"$PW$ [mm]", fontsize=label_fontsize)
    cbar.ax.tick_params(labelsize=tick_fontsize)

    """
    2. L128 COLD POOLS
    """
    for n in range(0,4):
        ax = fig.add_subplot(num_row, num_col, n+5)
        im = ax.pcolormesh(coldpools_all[n], vmin=vmin_cp_L128, vmax=vmax_cp_L128, cmap=cmap_cp)
        plt.title(plot_labels[n + 4],fontsize=label_fontsize)
        plt.gca().set_aspect('equal')

        if n == 0:
            ax.set_xlabel("x [km]",fontsize=label_fontsize)
            ax.set_ylabel("y [km]",fontsize=label_fontsize)
            ax.set_xticks([-0,50,100])
            ax.set_yticks([-0, 50, 100])
            ax.tick_params(axis='x', labelsize=tick_fontsize)
            ax.tick_params(axis='y', labelsize=tick_fontsize)
        else:
            ax.xaxis.set_tick_params(labelbottom=False)
            ax.yaxis.set_tick_params(labelleft=False)

    # COLORBAR
    cbar_ax = fig.add_axes([0.9, 0.54, 0.01, 0.15])  # 1 rows, [left, bottom, width, height]
    cbar = fig.colorbar(im, cax=cbar_ax, extend="both")
    cbar.ax.set_ylabel(r"${\theta}_v - \overline{{\theta}_v}$ [K]",fontsize=label_fontsize)
    cbar.ax.tick_params(labelsize=tick_fontsize)


    """
    3. L256 PW
    """
    for n in range(4,8):
        ax = fig.add_subplot(num_row, num_col, n+5)
        im = ax.pcolormesh(pw_all[n], vmin=vmin_pw_L256, vmax=vmax_pw_L256, cmap=cmap_pw)
        plt.title(plot_labels[n+4],fontsize=label_fontsize)
        plt.gca().set_aspect('equal')

        if n == 4:
            ax.set_xlabel("x [km]",fontsize=label_fontsize)
            ax.set_ylabel("y [km]",fontsize=label_fontsize)
            ax.set_xticks([-0,100,200])
            ax.set_yticks([-0,100, 200])
            ax.tick_params(axis='x', labelsize=tick_fontsize)
            ax.tick_params(axis='y', labelsize=tick_fontsize)
        else:
            ax.xaxis.set_tick_params(labelbottom=False)
            ax.yaxis.set_tick_params(labelleft=False)

    # COLORBAR
    cbar_ax = fig.add_axes([0.9, 0.3, 0.01, 0.15])
      # 1 rows, [left, bottom, width, height]
    cbar = fig.colorbar(im, cax=cbar_ax, extend="both")
    cbar.ax.set_ylabel(r"$PW$ [mm]",fontsize=label_fontsize)
    cbar.ax.tick_params(labelsize=tick_fontsize)

    """
    4. L256 COLD POOLS
    """
    for n in range(4, 8):
        ax = fig.add_subplot(num_row, num_col, n + 9)
        im = ax.pcolormesh(coldpools_all[n], vmin=vmin_cp_L256, vmax=vmax_cp_L256, cmap=cmap_cp)
        plt.title(plot_labels[n + 8], fontsize=label_fontsize)
        plt.gca().set_aspect('equal')

        if n == 4:
            ax.set_xlabel("x [km]", fontsize=label_fontsize)
            ax.set_ylabel("y [km]", fontsize=label_fontsize)
            ax.set_xticks([-0, 100, 200])
            ax.set_yticks([-0, 100, 200])
            ax.tick_params(axis='x', labelsize=tick_fontsize)
            ax.tick_params(axis='y', labelsize=tick_fontsize)
        else:
            ax.xaxis.set_tick_params(labelbottom=False)
            ax.yaxis.set_tick_params(labelleft=False)

    # COLORBAR
    cbar_ax = fig.add_axes([0.9, 0.055, 0.01, 0.15]) # 1 rows, [left, bottom, width, height]
    cbar = fig.colorbar(im, cax=cbar_ax, extend="both")
    cbar.ax.set_ylabel(r"${\theta}_v - \overline{{\theta}_v}$ [K]", fontsize=label_fontsize)
    cbar.ax.tick_params(labelsize=tick_fontsize)


    #-----------------------------------------
    fig.subplots_adjust(left=0.08,
                        bottom=0.03,
                        right=0.88,
                        top=0.96,
                        wspace=0.12,  # 0.02
                        hspace=0.2)  # 0.65

    """
    wspace and hspace specify the space reserved between Matplotlib subplots. They are the fractions of axis width and height, respectively.
    left, right, top and bottom parameters specify four sides of the subplots’ positions. They are the fractions of the width and height of the figure.
    top and bottom should add up to 1.0
    """

    plt.show()



def main():

    configs = ["alpha0_L128","alpha0.005_L128","alpha0.01_L128","alpha1_L128",
               "alpha0_L256","alpha0.05_L256","alpha0.1_L256","alpha1_L256"]
    plot_labels = [r"(a) ${\alpha}$=0",r"(b) ${\alpha}$=0.005",r"(c) ${\alpha}$=0.01",r"(d) ${\alpha}$=1, CTRL",
                   r"(e) ${\alpha}$=0",r"(f) ${\alpha}$=0.005",r"(g) ${\alpha}$=0.01",r"(h) ${\alpha}$=1, CTRL",
                   r"(i) ${\alpha}$=0", r"(j) ${\alpha}$=0.05", r"(k) ${\alpha}$=0.1", r"(l) ${\alpha}$=1, CTRL",
                   r"(m) ${\alpha}$=0", r"(n) ${\alpha}$=0.05", r"(o) ${\alpha}$=0.1", r"(p) ${\alpha}$=1, CTRL"
                   ]


    coldpools_all = []
    pw_all= []

    for config in configs:

        arr_cp = get_array(config,"cold_pools")
        arr_pw = get_array(config,"pw")

        coldpools_all.append(arr_cp)
        pw_all.append(arr_pw)

    print (len(coldpools_all),len(pw_all))

    plot_figure(coldpools_all,pw_all,plot_labels)


if __name__ == "__main__":

    main()