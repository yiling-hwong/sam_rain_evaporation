"""
@author: Yi-Ling HWONG
"""
import matplotlib as mpl
#mpl.rcParams['figure.dpi'] = 300
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
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

def get_array(config,var,plot_moving_average_flag,moving_average_window):

    file_path = "../data/dry_subs_timeseries/"+config+"_"+var+".csv"
    arr = np.loadtxt(file_path,delimiter=",")

    if plot_moving_average_flag == True:
        ma_mode = "same"
        window = np.ones(int(moving_average_window)) / float(moving_average_window)
        arr_ma = np.apply_along_axis(lambda m: np.convolve(m, window, mode=ma_mode), axis=1,arr=arr)  # shape = (53,40)(z,time)
        arr = arr_ma

    return arr

def get_list(config,var):

    file_path = "../data/dry_subs_timeseries/"+config+"_"+var+".csv"
    var_list = read_csv_to_list(file_path)

    if var == "z":
        var_list = [x/1000 for x in var_list]

    return var_list

def plot_figure(w_arrays,qv_arrays,vert_flux_arrays,z_list,time_list):

    title_fontsize = 16
    label_fontsize = 14
    tick_fontsize = 12
    bigframe_pad = 28

    x_lim_min = 1
    x_lim_max = 40  # days
    y_lim_min = 0
    y_lim_max = 5  # height in km, 5 or 3
    y_max_ind = 21 # 20 or 21 for 5km; 8 for 1km; 12 for 2km; 14 for 3km

    #######
    # PLOT
    #######

    fig, big_axes = plt.subplots(figsize=(11, 6), nrows=2, ncols=1, sharey=True)
    gs = gridspec.GridSpec(2, 3)

    for row, big_ax in enumerate(big_axes, start=1):

        # plt.rcParams['text.usetex'] = True
        # plt.rcParams["mathtext.fontset"] = "custom"

        if row == 1:
            #big_ax.set_title(r"L128_$\mathbf{{\alpha}0}$", fontsize=title_fontsize, weight = "bold", pad=bigframe_pad)
            big_ax.set_title(r"L128_${\alpha}0$", fontsize=title_fontsize, pad=bigframe_pad)
        if row == 2:
            #big_ax.set_title(r"L256_$\mathbf{{\alpha}0}$", fontsize=title_fontsize, weight = "bold", pad=bigframe_pad)
            big_ax.set_title(r"L256_${\alpha}0$", fontsize=title_fontsize, pad=bigframe_pad)

        # Turn off axis lines and ticks of the big subplot
        # obs alpha is 0 in RGBA string!
        big_ax.set_xticks([])
        big_ax.set_yticks([])
        # removes the white frame
        big_ax._frameon = False


    """
    1. L128 W
    """
    ax = fig.add_subplot(gs[0, 0])

    arr_to_plot = w_arrays[0][:y_max_ind,:]
    z_to_plot = z_list[0][:y_max_ind]
    cplot = ax.contourf(time_list[0], z_to_plot, arr_to_plot, cmap="seismic", levels=50, extend="both", norm=MidpointNormalize(midpoint=0))  # vmin=minzz, vmax=maxzz)
    ax.set_ylabel('z [km]', color='black',fontsize=label_fontsize)
    ax.set_xlim([x_lim_min, x_lim_max])
    ax.set_ylim([y_lim_min, y_lim_max])
    ax.tick_params(axis='x', labelsize=tick_fontsize)
    ax.tick_params(axis='y', labelsize=tick_fontsize)
    ax.xaxis.set_tick_params(labelbottom=False)

    cb = plt.colorbar(cplot, orientation='vertical', shrink=0.8, aspect=20)
    cb.ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    cb.ax.yaxis.set_offset_position('left')
    cb.ax.set_yticks([-0.012, -0.006, 0, 0.006, 0.012])
    cb.ax.tick_params(axis='y', labelsize=tick_fontsize)
    cb.set_label(label='[m s$^{-1}$]',fontsize=tick_fontsize)
    plt.title("(a) w$_{dry}$",fontsize=title_fontsize)

    """
    2. L128 QV
    """
    ax = fig.add_subplot(gs[0, 1])

    arr_to_plot = qv_arrays[0][:y_max_ind,:]
    z_to_plot = z_list[0][:y_max_ind]
    cplot = ax.contourf(time_list[0], z_to_plot, arr_to_plot, cmap="rainbow", levels=50,extend="both")  # vmin=minzz, vmax=maxzz)
    ax.set_xlim([x_lim_min, x_lim_max])
    ax.set_ylim([y_lim_min, y_lim_max])
    ax.xaxis.set_tick_params(labelbottom=False)
    ax.yaxis.set_tick_params(labelleft=False)

    cb = plt.colorbar(cplot, orientation='vertical', shrink=0.8, aspect=20)
    cb.ax.set_yticks([4,8,12])
    cb.ax.tick_params(axis='y', labelsize=tick_fontsize)
    cb.set_label(label='[g kg$^{-1}$]',fontsize=tick_fontsize)
    plt.title("(b) q$_{v,dry}$",fontsize=title_fontsize)


    """
    3. L128 W_DQ/DZ
    """
    ax = fig.add_subplot(gs[0, 2])

    arr_to_plot = -vert_flux_arrays[0][:y_max_ind,:]
    z_to_plot = z_list[0][:y_max_ind]
    cplot = ax.contourf(time_list[0], z_to_plot, arr_to_plot, cmap="seismic", levels=50,extend="both", norm=MidpointNormalize(midpoint=0))  # vmin=minzz, vmax=maxzz)
    ax.set_xlim([x_lim_min, x_lim_max])
    ax.set_ylim([y_lim_min, y_lim_max])
    ax.xaxis.set_tick_params(labelbottom=False)
    ax.yaxis.set_tick_params(labelleft=False)

    cb = plt.colorbar(cplot, orientation='vertical', shrink=0.8, aspect=20)
    cb.ax.set_yticks([-20,-10,0,10])
    cb.ax.tick_params(axis='y', labelsize=tick_fontsize)
    cb.set_label(label='[g kg$^{-1}$ d$^{-1}$]',fontsize=tick_fontsize)
    #plt.title(r"(c) $\dfrac{\partial(wq_{v,dry})}{\partial{z}}$, vert. flux",fontsize=title_fontsize)
    plt.title(r"(c) $\dfrac{\partial(q_{v,dry})}{\partial{t}}$, vert. flux", fontsize=title_fontsize)

#####################
# L256
#####################

    """
    4. L256 W
    """
    ax = fig.add_subplot(gs[1, 0])

    arr_to_plot = w_arrays[1][:y_max_ind,:]
    z_to_plot = z_list[1][:y_max_ind]
    cplot = ax.contourf(time_list[1], z_to_plot, arr_to_plot, cmap="seismic", levels=50, extend="min", norm=MidpointNormalize(midpoint=0))  # vmin=minzz, vmax=maxzz)
    ax.set_ylabel('z [km]', color='black',fontsize=label_fontsize)
    ax.set_xlabel('Time [day]', color='black',fontsize=label_fontsize)
    ax.set_xlim([x_lim_min, x_lim_max])
    ax.set_ylim([y_lim_min, y_lim_max])
    ax.tick_params(axis='x', labelsize=tick_fontsize)
    ax.tick_params(axis='y', labelsize=tick_fontsize)

    cb = plt.colorbar(cplot, orientation='vertical', shrink=0.8, aspect=20)
    cb.ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    cb.ax.yaxis.set_offset_position('left')
    cb.ax.set_yticks([-0.01, -0.005, 0])
    cb.ax.tick_params(axis='y', labelsize=tick_fontsize)
    cb.set_label(label='[m s$^{-1}$]',fontsize=tick_fontsize)
    plt.title("(d) w$_{dry}$",fontsize=title_fontsize)

    """
    5. L256 QV
    """
    ax = fig.add_subplot(gs[1, 1])


    arr_to_plot = qv_arrays[1][:y_max_ind, :]
    z_to_plot = z_list[1][:y_max_ind]
    cplot = ax.contourf(time_list[1], z_to_plot, arr_to_plot, cmap="rainbow", levels=50,extend="both")  # vmin=minzz, vmax=maxzz)
    ax.set_xlabel('Time [day]', color='black',fontsize=label_fontsize)
    ax.set_xlim([x_lim_min, x_lim_max])
    ax.set_ylim([y_lim_min, y_lim_max])
    ax.tick_params(axis='x', labelsize=tick_fontsize)
    ax.tick_params(axis='y', labelsize=tick_fontsize)
    ax.yaxis.set_tick_params(labelbottom=False)

    cb = plt.colorbar(cplot, orientation='vertical', shrink=0.8, aspect=20)
    cb.ax.set_yticks([4,8,12])
    cb.ax.tick_params(axis='y', labelsize=tick_fontsize)
    cb.set_label(label='[g kg$^{-1}$]',fontsize=tick_fontsize)
    plt.title("(e) q$_{v,dry}$",fontsize=title_fontsize)

    """
    6. L256 W_DQ/DZ
    """
    ax = fig.add_subplot(gs[1, 2])

    arr_to_plot = -vert_flux_arrays[1][:y_max_ind, :]
    z_to_plot = z_list[1][:y_max_ind]
    cplot = ax.contourf(time_list[1], z_to_plot, arr_to_plot, cmap="seismic", levels=50,extend="both",norm=MidpointNormalize(midpoint=0))  # vmin=minzz, vmax=maxzz)
    ax.set_xlabel('Time [day]', color='black',fontsize=label_fontsize)
    ax.set_xlim([x_lim_min, x_lim_max])
    ax.set_ylim([y_lim_min, y_lim_max])
    ax.tick_params(axis='x', labelsize=tick_fontsize)
    ax.tick_params(axis='y', labelsize=tick_fontsize)
    ax.yaxis.set_tick_params(labelbottom=False)

    cb = plt.colorbar(cplot, orientation='vertical', shrink=0.8, aspect=20)
    cb.ax.set_yticks([-10,0,10,20,30])
    cb.ax.tick_params(axis='y', labelsize=tick_fontsize)
    cb.set_label(label='[g kg$^{-1}$ d$^{-1}$]',fontsize=tick_fontsize)
    #plt.title(r"(f) $\dfrac{\partial(wq_{v,dry})}{\partial{z}}$, vert. flux", fontsize=title_fontsize)
    plt.title(r"(f) $\dfrac{\partial(q_{v,dry})}{\partial{t}}$, vert. flux", fontsize=title_fontsize)


    #-----------------------------------------
    fig.subplots_adjust(left=0.05,
                        bottom=0.1,
                        right=0.95,
                        top=0.88,
                        wspace=0.2,  # 0.02
                        hspace=0.45)  # 0.65

    """
    wspace and hspace specify the space reserved between Matplotlib subplots. They are the fractions of axis width and height, respectively.
    left, right, top and bottom parameters specify four sides of the subplots’ positions. They are the fractions of the width and height of the figure.
    top and bottom should add up to 1.0
    """

    plt.show()

def main():

    configs = ["alpha0_L128", "alpha0_L256"]
    plot_moving_average_flag = True
    moving_average_window = 10

    # GET DATA
    w_arrays  = []
    qv_arrays = []
    vert_flux_arrays = [] # dq/dt vertical flux
    z_list = []
    time_list = []
    for config in configs:
        arr = get_array(config,"W_array",plot_moving_average_flag,moving_average_window)
        w_arrays.append(arr)

        arr = get_array(config,"qv_array",plot_moving_average_flag,moving_average_window)
        qv_arrays.append(arr)

        arr = get_array(config,"W_dqdz_array",plot_moving_average_flag,moving_average_window)
        vert_flux_arrays.append(arr)

        z = get_list(config,"z")
        z_list.append(z)

        time = get_list(config,"time")
        time_list.append(time)

    print (len(w_arrays),len(qv_arrays), len(vert_flux_arrays), len(z_list), len(time_list))

    plot_figure(w_arrays,qv_arrays,vert_flux_arrays,z_list,time_list)




if __name__ == "__main__":

    main()