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
    bigframe_pad = 6

    vmin_cp = -1
    vmax_cp = -vmin_cp

    vmin_pw = 25
    vmax_pw = 45

    title_fontsize = 16
    label_fontsize = 14
    tick_fontsize = 13

    cmap_cp = "seismic"
    cmap_pw = "rainbow"

    fig, big_axes = plt.subplots(figsize=(10, 12), nrows=4, ncols=2, sharey=True)

    for row in range(4):
        for col in range(2):

            #HomoRad_HomoFlx & HomoRad_IntFlx
            if row == 0 and col == 0:
                big_axes[row,col].set_title("HomoRad, HomoFlx", fontsize=title_fontsize, pad=bigframe_pad)
            if row == 0 and col == 1:
                big_axes[row, col].set_title("HomoRad, IntFlx", fontsize=title_fontsize, pad=bigframe_pad)
            if row == 1 and col == 0:
                big_axes[row,col].set_title("HomoRad, HomoFlx", fontsize=title_fontsize, pad=bigframe_pad)
            if row == 1 and col == 1:
                big_axes[row, col].set_title("HomoRad, IntFlx", fontsize=title_fontsize, pad=bigframe_pad)

            # IntRad_HomoFlx & IntRad_IntFlx
            if row == 2 and col == 0:
                big_axes[row,col].set_title("IntRad, HomoFlx", fontsize=title_fontsize, pad=bigframe_pad)
            if row == 2 and col == 1:
                big_axes[row, col].set_title("IntRad, IntFlx", fontsize=title_fontsize, pad=bigframe_pad)
            if row == 3 and col == 0:
                big_axes[row,col].set_title("IntRad, HomoFlx", fontsize=title_fontsize, pad=bigframe_pad)
            if row == 3 and col == 1:
                big_axes[row, col].set_title("IntRad, IntFlx", fontsize=title_fontsize, pad=bigframe_pad)

            # Turn off axis lines and ticks of the big subplot
            # obs alpha is 0 in RGBA string!
            #big_ax.tick_params(labelcolor=(1., 1., 1., 0.0), top='off', bottom='off', left='off', right='off')
            big_axes[row, col].set_xticks([])
            big_axes[row, col].set_yticks([])
            # removes the white frame
            big_axes[row, col]._frameon = False


    """
    1. PW HomoRad_HomoFlx & HomoRad_IntFlx
    """
    for n in range(0, 4):
        ax = fig.add_subplot(num_row, num_col, n + 1)
        im = ax.pcolormesh(pw_all[n], vmin=vmin_pw, vmax=vmax_pw, cmap=cmap_pw)
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
    cbar_ax = fig.add_axes([0.9, 0.80, 0.01, 0.15])
      # 1 rows, [left, bottom, width, height]
    cbar = fig.colorbar(im, cax=cbar_ax, extend="both")
    cbar.ax.set_ylabel(r"$PW$ [mm]", fontsize=label_fontsize)
    cbar.ax.tick_params(labelsize=tick_fontsize)

    """
    2. COLD POOLS HomoRad_HomoFlx & HomoRad_IntFlx
    """
    for n in range(0,4):
        ax = fig.add_subplot(num_row, num_col, n+5)
        im = ax.pcolormesh(coldpools_all[n], vmin=vmin_cp, vmax=vmax_cp, cmap=cmap_cp)
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
    cbar_ax = fig.add_axes([0.9, 0.55, 0.01, 0.15])  # 1 rows, [left, bottom, width, height]
    cbar = fig.colorbar(im, cax=cbar_ax, extend="both")
    cbar.ax.set_ylabel(r"${\theta}_v - \overline{{\theta}_v}$ [K]",fontsize=label_fontsize)
    cbar.ax.tick_params(labelsize=tick_fontsize)


    """
    3. PW IntRad_HomoFlx & IntRad_IntFlx
    """
    for n in range(4,8):
        ax = fig.add_subplot(num_row, num_col, n+5)
        im = ax.pcolormesh(pw_all[n], vmin=vmin_pw, vmax=vmax_pw, cmap=cmap_pw)
        plt.title(plot_labels[n+4],fontsize=label_fontsize)
        plt.gca().set_aspect('equal')

        if n == 4:
            ax.set_xlabel("x [km]",fontsize=label_fontsize)
            ax.set_ylabel("y [km]",fontsize=label_fontsize)
            ax.set_xticks([-0, 50, 100])
            ax.set_yticks([-0, 50, 100])
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
    4. COLD POOLS IntRad_HomoFlx & IntRad_IntFlx
    """
    for n in range(4, 8):
        ax = fig.add_subplot(num_row, num_col, n + 9)
        im = ax.pcolormesh(coldpools_all[n], vmin=vmin_cp, vmax=vmax_cp, cmap=cmap_cp)
        plt.title(plot_labels[n + 8], fontsize=label_fontsize)
        plt.gca().set_aspect('equal')

        if n == 4:
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
    cbar_ax = fig.add_axes([0.9, 0.055, 0.01, 0.15]) # 1 rows, [left, bottom, width, height]
    cbar = fig.colorbar(im, cax=cbar_ax, extend="both")
    cbar.ax.set_ylabel(r"${\theta}_v - \overline{{\theta}_v}$ [K]", fontsize=label_fontsize)
    cbar.ax.tick_params(labelsize=tick_fontsize)


    #-----------------------------------------
    fig.subplots_adjust(left=0.08,
                        bottom=0.03,
                        right=0.88,
                        top=0.97,
                        wspace=0.12,  # 0.02
                        hspace=0.3)  # 0.65

    """
    wspace and hspace specify the space reserved between Matplotlib subplots. They are the fractions of axis width and height, respectively.
    left, right, top and bottom parameters specify four sides of the subplots’ positions. They are the fractions of the width and height of the figure.
    top and bottom should add up to 1.0
    """

    plt.show()



def main():

    configs = ["alpha0.005_L128","alpha0.01_L128","alpha0.01_L128_HomoRad_IntFlx","alpha0.02_L128_HomoRad_IntFlx",
               "alpha0.1_L128_IntRad_HomoFlx","alpha0.125_L128_IntRad_HomoFlx","alpha0.125_L128_IntRad_IntFlx","alpha0.15_L128_IntRad_IntFlx"]
    plot_labels = [r"(a) ${\alpha}$=0.005",r"(b) ${\alpha}$=0.01",r"(c) ${\alpha}$=0.01",r"(d) ${\alpha}$=0.02",
                   r"(e) ${\alpha}$=0.005",r"(f) ${\alpha}$=0.01",r"(g) ${\alpha}$=0.01",r"(h) ${\alpha}$=0.02",
                   r"(i) ${\alpha}$=0.1", r"(j) ${\alpha}$=0.125", r"(k) ${\alpha}$=0.125", r"(l) ${\alpha}$=0.15",
                   r"(m) ${\alpha}$=0.1", r"(n) ${\alpha}$=0.125", r"(o) ${\alpha}$=0.125", r"(p) ${\alpha}$=0.15"
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