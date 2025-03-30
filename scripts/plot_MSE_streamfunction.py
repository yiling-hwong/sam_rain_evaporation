"""
@author: Yi-Ling HWONG
"""
import matplotlib as mpl
#mpl.rcParams['figure.dpi'] = 300
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from applications.utils import *

class MidpointNormalize(colors.Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))

def get_array(config,var):

    file_path = "../data/MSE_streamfunction/"+config+"_"+var+".csv"
    arr = np.loadtxt(file_path,delimiter=",")

    return arr

def plot_figure(z_all,crh_all,mse_all,psi_all,plot_labels):

    """
    Set plot parameters here
    """

    num_row = 2
    num_col = 3
    end_level = 45
    bigframe_pad = 25
    cmap = "jet"

    title_fontsize = 16
    label_fontsize = 14
    tick_fontsize = 12

    num_level_mse = 50
    vmin_mse = 310
    vmax_mse = 340

    #######
    # PLOT
    #######

    fig, big_axes = plt.subplots(figsize=(10, 7), nrows=2, ncols=1, sharey=True)

    for row, big_ax in enumerate(big_axes, start=1):

        if row == 1:
            big_ax.set_title("L128", fontsize=title_fontsize, pad=bigframe_pad)
        if row == 2:
            big_ax.set_title("L256", fontsize=title_fontsize, pad=bigframe_pad)

        # Turn off axis lines and ticks of the big subplot
        # obs alpha is 0 in RGBA string!
        big_ax.set_xticks([])
        big_ax.set_yticks([])
        # removes the white frame
        big_ax._frameon = False

    """
    CONTOUR PLOTS (ALL)
    """

    for n in range(0,6):
        ax = fig.add_subplot(num_row, num_col, n + 1)

        # MSE
        mse_plot = ax.contourf(crh_all[n][1:], z_all[n][:end_level] / 1000, mse_all[n][:end_level, 1:] / 1000, num_level_mse, cmap=cmap, extend="both",
                               vmin=vmin_mse,vmax=vmax_mse)
        ax.set_xlabel('CRH ', color='black',fontsize=label_fontsize)
        ax.tick_params(axis='x', labelsize=tick_fontsize)
        plt.title(plot_labels[n], fontsize=label_fontsize)

        # STREAM FUNCTION CONTOURS

        # get interval for streamfunction

        if n < 3:
            sf_interval = 0.02  # in kg m-2 s-1 (0.02 for smalldom, 0.04 for bigdom)
            sf_min = round(np.amin(-psi_all[n][:end_level, :]), 2)
            sf_max = np.amax(-psi_all[n][:end_level, :])
        else:
            sf_interval = 0.04
            sf_min = round(np.amin(-psi_all[n][:end_level, :]), 2)
            sf_max = np.amax(-psi_all[n][:end_level, :])
        sf_levels = np.arange(sf_min, sf_max + sf_interval, sf_interval)
        print ()
        print("Streamfunction min and max:", sf_min, sf_max)
        print(sf_levels)

        contours=ax.contour(crh_all[n][1:],z_all[n][:end_level]/1000,-psi_all[n][:end_level,:],levels = sf_levels,colors='black',linewidths=0.6,norm=MidpointNormalize(midpoint=0))
        #ax.clabel(contours, inline=True, fontsize=8) #uncomment to show psy values inline

        if n == 0 or n == 3:
            ax.set_ylabel("Height [km]",fontsize=label_fontsize)
            ax.set_yticks([2,4,6,8,10,12,14])
            ax.tick_params(axis='y', labelsize=tick_fontsize)
        else:
            ax.yaxis.set_tick_params(labelleft=False)
            ax.set(ylabel=None)

        if n == 3: #alpha0_L256
            ax.set_xticks([0.3,0.4, 0.5, 0.6, 0.7, 0.8])

            # COLORBAR
            cbar_ax = fig.add_axes([0.915, 0.3, 0.01, 0.4])  # 1 rows, [left, bottom, width, height]
            cbar = fig.colorbar(mse_plot, cax=cbar_ax, extend="both",
                                ticks=np.linspace(np.rint(vmin_mse), np.rint(vmax_mse), 4))
            cbar.ax.set_ylabel("MSE [kJ/kg]", fontsize=label_fontsize)
            cbar.ax.tick_params(labelsize=tick_fontsize)

        if n == 4: #alpha0.05_L256
            ax.set_xticks([0.3, 0.4, 0.5, 0.6, 0.7])

    #-----------------------------------------
    fig.subplots_adjust(left=0.07,
                        bottom=0.1,
                        right=0.9,
                        top=0.9,
                        wspace=0.05,  # 0.02
                        hspace=0.55)  # 0.65

    """
    wspace and hspace specify the space reserved between Matplotlib subplots. They are the fractions of axis width and height, respectively.
    left, right, top and bottom parameters specify four sides of the subplots’ positions. They are the fractions of the width and height of the figure.
    top and bottom should add up to 1.0
    """

    plt.show()


def main():

    configs = ["alpha0_L128", "alpha0.005_L128", "alpha0.01_L128",
               "alpha0_L256", "alpha0.05_L256", "alpha0.1_L256"]
    plot_labels = [r"(a) ${\alpha}$=0",r"(b) ${\alpha}$=0.005",r"(c) ${\alpha}$=0.01",
                   r"(d) ${\alpha}$=0",r"(e) ${\alpha}$=0.05",r"(f) ${\alpha}$=0.1"]

    z_all = []
    crh_all = []
    mse_all = []
    psi_all = [] # stream function

    for config in configs:

        arr_z = get_array(config,"z")
        arr_crh = get_array(config,"crh")
        arr_mse = get_array(config,"mse")
        arr_psi = get_array(config,"streamfunction")

        z_all.append(arr_z)
        crh_all.append(arr_crh)
        mse_all.append(arr_mse)
        psi_all.append(arr_psi)

    plot_figure(z_all,crh_all,mse_all,psi_all,plot_labels)


if __name__ == "__main__":

    main()
