import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
from fonts import *

def error_cal(bin_edges,weights,data):
    errors = []
    bin_centers = []
    
    for bin_index in range(len(bin_edges) - 1):

        # find which data points are inside this bin
        bin_left = bin_edges[bin_index]
        bin_right = bin_edges[bin_index + 1]
        in_bin = np.logical_and(bin_left < data, data <= bin_right)
        

        # filter the weights to only those inside the bin
        weights_in_bin = weights[in_bin]

        # compute the error however you want
        error = np.sqrt(np.sum(weights_in_bin ** 2))
        errors.append(error)

        # save the center of the bins to plot the errorbar in the right place
        bin_center = (bin_right + bin_left) / 2
        bin_centers.append(bin_center)

    errors=np.asarray(errors)
    bin_centers=np.asarray(bin_centers)
    return errors, bin_centers


def onedimension_hist(ax,x,weights,xlog,ylog,bins_start,bins_stop,bins,error,color,**kwargs):
    
    
    if xlog:
        bins=np.logspace(bins_start,bins_stop,bins)
        
    else:
         bins=np.linspace(bins_start,bins_stop,bins)
        
    counts, edges = np.histogram(x,weights=weights,bins=bins)
    ax.plot(edges,np.append(counts,counts[-1]),
                 drawstyle="steps-post",color=color,
                 **kwargs)
    print('Total counts are',np.sum(counts))
    if error:
        
        error,bin_centres = error_cal(edges,weights,x)
        
        ax.errorbar(x=bin_centres, y=counts,
                 yerr=error, color=color,fmt='o', markersize=8,capsize=5)
    if xlog:
        ax.set_xscale('log')
        
        
    if ylog:
        ax.set_yscale('log')
    
    for item in (ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(18)
        item.set_family('serif')

    ax.tick_params(axis='both',which='major',width=3,length=15,direction='in')
    ax.tick_params(axis='both',which='minor',width=1,length=8,direction='in')
    
    ax.grid(True, which="both", ls=":")
    
def weighted_quantile(
    values,
    quantiles,
    sample_weight=None,
    values_sorted=False,
):
    """
    Very close to numpy.percentile, but supports weights. Qantiles should be in [0, 1]!
    
    Parameters
    ----------
    values : array of floats
        Input data.
    quantiles : array of floats
        Quantile values to compute.
    sample_weight : array of floats
        Weights of the input data.
    values_sorted : bool
        Are the input values sorted, or not.

    Returns
    -------
    quantiles : array of floats
        Computed quantiles.
    """

    import numpy as np

    values = np.array(values)
    quantiles = np.array(quantiles)

    if values.size == 0: return(np.nan)

    if sample_weight is None:
        sample_weight = np.ones(len(values))

    sample_weight = np.array(sample_weight)

    assert np.all(quantiles >= 0) and np.all(quantiles <= 1), 'quantiles should be in [0, 1]'   

    if not values_sorted:
        sorter = np.argsort(values)
        values = values[sorter]
        sample_weight = sample_weight[sorter]

    weighted_quantiles = np.cumsum(sample_weight) - 0.5 * sample_weight
    weighted_quantiles /= np.sum(sample_weight)

    return np.interp(quantiles, weighted_quantiles, values) 
    

def plot_hist_errorbar(ax, hist, bins, yerror, **kwargs):
    l = ax.plot(bins,
                np.append(hist, hist[-1]),
                drawstyle="steps-post",
                **kwargs)

    if not np.all(yerror == 0.):
        
        ax.errorbar(np.mean(rolling_window(bins, 2), axis=1), hist, yerr=yerror,
                    ls="none", capsize=4.0, capthick=2,
                    color=l[0].get_color(), alpha=l[0].get_alpha())

    return l


def plot_hist_errorbar_unbinned(ax, samples, bins, weights,error, **kwargs):
    hist,yerror = make_hist_error(samples, bins, weights=weights)
    l = ax.plot(bins,
                np.append(hist, hist[-1]),
                drawstyle="steps-post",
                **kwargs)

    if not np.all(yerror == 0.):
        
        ax.errorbar(np.mean(rolling_window(bins, 2), axis=1), hist, yerr=yerror,
                    ls="none", capsize=4.0, capthick=2,
                    color=l[0].get_color(), alpha=l[0].get_alpha())

    return l

def plot_hist_errorbar_T(ax, hist, bins, error, **kwargs):
    
    l = ax.plot(np.append(hist, hist[-1]),
                bins,
                drawstyle="steps-pre",
                **kwargs)

    if not np.all(error == 0.):
        # skip errorbars that look weird if alpha!=0 is used
        ax.errorbar(hist, np.mean(rolling_window(bins, 2), axis=1), xerr=error,
                    ls="none", capsize=4.0, capthick=2,
                    color=l[0].get_color(), alpha=l[0].get_alpha())

    return l


def plot_hist(ax, hist, bins,weights=None, normed=False, error=False, **kwargs):
#     hist = make_hist_error(samples, bins, weights=weights, normed=False, error=False)
    l = ax.plot(bins,
                np.append(hist, hist[-1]),
                drawstyle="steps-post",
                **kwargs)

    

    return l

def plot_hist_band(ax, hist, bins, yerror, **kwargs):

    p = ax.stairs(hist, bins, baseline=None, **kwargs)
    ax.stairs(
        hist + yerror,
        bins,
        baseline=hist - yerror,
        fill=True,
        alpha=0.25,
        color=p.get_edgecolor()
    )

    return p


def plot_ratio_band(
    ax, hist, bins, yerror, hist_baseline, yerr_baseline=None, **kwargs
):

    ratio = hist / hist_baseline

    if yerr_baseline is not None:
        # propagate error on ratio
        yerror_ratio = ratio * np.sqrt(
            (yerror / hist)**2 + (yerr_baseline / hist_baseline)**2
        )
    else:
        yerror_ratio = yerror / hist_baseline

    p = ax.stairs(ratio, bins, baseline=None, **kwargs)
    ax.stairs(
        ratio + yerror_ratio,
        bins,
        baseline=ratio - yerror_ratio,
        fill=True,
        alpha=0.25,
        color=p.get_edgecolor()
    )

    return p

def plot_ratio_errorbar(ax, hist, bins, yerror, hist_baseline, yerr_baseline=None, **kwargs):
    """
    
    """
    ratio = hist/hist_baseline

    if yerr_baseline is not None:
        # propagate error on ratio
        yerror_ratio = ratio * np.sqrt( (yerror/hist)**2 + (yerr_baseline/hist_baseline)**2 )
    else:
        yerror_ratio = yerror/hist_baseline

    l = plot_hist_errorbar(ax,ratio,bins,yerror_ratio, **kwargs)

    return l


def plot_data_hist_errorbar(ax, hist, bins, yerror, **kwargs):
    """Draw data as points, not histogram"""

    l = ax.errorbar(np.mean(rolling_window(bins, 2), axis=1), hist, yerr=yerror,
                 ls="none", capsize=4.0, capthick=2,
                 **kwargs)

    return l


def plot_data_ratio_errorbar(ax, hist, bins, yerror, hist_baseline, yerr_baseline=None, **kwargs):
    """Draw data as points, not histogram"""
    ratio = hist/hist_baseline

    if yerr_baseline is not None:
        # propagate error on ratio
        yerror_ratio = ratio * np.sqrt( (yerror/hist)**2 + (yerr_baseline/hist_baseline)**2 )
    else:
        yerror_ratio = yerror/hist_baseline

    l = plot_data_hist_errorbar(ax,ratio,bins,yerror_ratio, **kwargs)

    return l


# wrapper for numpy.histogram for calculating uncertainty from weighted entries
def make_hist_error(samples, bins, weights=None, normed=False):
    # set weights if not give
    if weights is None:
        weights = np.ones_like(samples, dtype=int)

    hist, bins = np.histogram(samples, bins, weights=weights)
    if normed:
        norm = 1./np.diff(bins)/hist.sum()
        hist = hist*norm

    # calculate error
    weights_hist, _ = np.histogram(samples, bins, weights=np.power(weights, 2))
    yerror = np.sqrt(weights_hist)

    if normed:
        yerror = yerror*norm

    return hist, yerror


# wrapper for numpy.histogram2d for calculating uncertainty from weighted entries
def make_hist2d_error(x, y, binning, weights=None, normed=False):
    # set weights if not give
    if weights is None:
        weights = np.ones_like(x)

    hist, _, _ = np.histogram2d(x, y, bins=binning, weights=weights)
    if normed:
        raise NotImplmentedError
        # norm = 1./np.diff(bins)/hist.sum()
        # hist = hist*norm

    # calculate error
    weights_hist, _, _ = np.histogram2d(x, y, bins=binning, weights=np.power(weights, 2))
    yerror = np.sqrt(weights_hist)

    # if normed:
    #     yerror = yerror*norm

    return hist, yerror


def initialize_figure(figsize=(12,9), height_ratios=(3,1)):
    # make the figure
    if height_ratios is not None:
        fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw={"height_ratios": height_ratios}, figsize=figsize)
    else:
        fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw={"height_ratios": (1,1)}, figsize=figsize)

    # labels
    ax1.set_ylabel(r"$N_{\mathrm{Events}}$")
    ax2.set_ylabel(r"ratio")

    #ax2.set_ylim(0.5, 2.)
    fig.align_labels()

    return fig, (ax1, ax2)


# save figure in different formats
def savefig(fig, path, name, plot_format=["png", "pdf"], bbox_inches="tight", dpi=200, **kwargs):

    if type(plot_format) == str:
        if not plot_format.startswith("."):
            plot_format = "." + plot_format
        fig.savefig(os.path.join(path, name+plot_format), bbox_inches=bbox_inches, dpi=dpi, **kwargs)

    elif type(plot_format) == list:
        for iformat in plot_format:
            if not iformat.startswith("."):
                iformat = "." + iformat
            fig.savefig(os.path.join(path, name+iformat), bbox_inches=bbox_inches, dpi=dpi, **kwargs)
    else:
        raise NotImplementedError
    return


# adapted from stackoverflow: https://stackoverflow.com/questions/36699155/how-to-get-color-of-most-recent-plotted-line-in-pythons-plt
# can also be used for calculating bin centers
def rolling_window(a, window):
    shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
    strides = a.strides + (a.strides[-1],)
    return np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)



def plot_3D_analysisvariables(
    mu_3d, ssq_3d, data_3d, bins_energy, bins_length,bins_eratio,
    ylim_energy_ratio=(0, 2),
    ylim_length_ratio=(0.5, 1.5),
    ylim_eratio_ratio=(0, 2),
    components=None,
    data_label="data",
    plot_title=None,
    save=None, plot_dir=None
):
    """
    _summary_

    Parameters
    ----------
    mu_3d : np.array
        3d MC bincount
    ssq_3d : np.array
        3d expected MC fluctuations
    data_3d : np.array
        3d data bin count
    bins_energy : np.array
        energy bins
    bins_length : np.array
        length bins
    bins_eratio : np.array
        eratio bins
    ylim_energy_ratio : tuple
        ylim in energy spectrum ratio plot
    ylim_length_ratio : tuple
        ylim in length spectrum ratio plot
    ylim_eratio_ratio : tuple
        ylim in eratio spectrum ratio plot
    components : dict, optional
        dict(
           'settings' -> plot settings per comp,
           'hists' -> mu, ssq per comp for corresponding evaluated parameters
         ), by default None so that no individual components will be drawn
    data_label : str, optional
        label for data, by default "data"
    save : str, optional
        name for plot to save, by default None
    plot_dir : str, optional
        dir for plot to save, by default None
    """

    # energy spectrum

    fig, (ax1, ax2) = initialize_figure(figsize=(8,8))
    import matplotlib.font_manager as font_manager
    font_axis_label = {'family': 'serif',
            'color':  'black',
            'weight': 'normal',
            'size': 15,
            }
    font_title = {'family': 'serif',
            'color':  'black',
            'weight': 'bold',
            'size': 18,
            }
    font_legend = font_manager.FontProperties(family='serif',
                                       weight='normal',
                                       style='normal', size=12)
    if plot_title is not None:
        fig.suptitle(plot_title,fontdict=font_title)
    
    #
    # energy
    #
    hist_mc = mu_3d.sum(axis=(1,2))
    error_mc = np.sqrt(ssq_3d.sum(axis=(1,2)))

    hist_data = data_3d.sum(axis=(1,2))
    error_data = np.sqrt(data_3d.sum(axis=(1,2)))

    hist_ratio_base = hist_mc
    # yerr_ratio_base = yerr_base

    plot_data_hist_errorbar(
        ax1, hist_data, bins_energy, error_data, label=data_label, color='k', marker='o'
    )
    plot_data_ratio_errorbar(
        ax2,
        hist_data,
        bins_energy,
        error_data,
        hist_ratio_base,
        yerr_baseline=None,
        color='k', marker='o'
    )

    plot_hist_errorbar(ax1, hist_mc, bins_energy, error_mc, label="MC sum")
    plot_ratio_errorbar(
        ax2,
        hist_mc,
        bins_energy,
        error_mc,
        hist_ratio_base,
        yerr_baseline=None
    )
    if components is not None:
        for comp, d in components['settings'].items():
            hist = components['hists'][comp]["mu"].sum(axis=(1,2))
            yerror = np.sqrt(components['hists'][comp]["ssq"].sum(axis=(1,2)))
            plot_hist_errorbar(ax1, hist, bins_energy, yerror, **d['plot_settings'])

    # setup figure
    plt.xlabel("Reco energy [GeV]",fontdict=font_axis_label)
    ax1.set_ylabel(r"$N_{\mathrm{Events}}$",fontdict=font_axis_label)
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    #restrict_ylim_data(ax1)
    ax1.set_xlim(bins_energy[0], bins_energy[-1])
    ax1.set_ylim(1e-3,5)
    ax2.set_ylim(ylim_energy_ratio[0], ylim_energy_ratio[1])
    for item in (ax1.get_xticklabels() + ax1.get_yticklabels()):
        item.set_fontsize(12)
        item.set_family('serif')
    for item in (ax2.get_xticklabels() + ax2.get_yticklabels()):
        item.set_fontsize(12)
        item.set_family('serif')
    ax1.tick_params(axis='both',which='major',width=3,length=15,direction='in')
    ax1.tick_params(axis='both',which='minor',width=1,length=8,direction='in')
    ax2.tick_params(axis='both',which='major',width=3,length=15,direction='in')
    ax2.tick_params(axis='both',which='minor',width=1,length=8,direction='in')
    # save and show plot
    ax1.legend(prop=font_legend)
    fig.tight_layout()
    fig.subplots_adjust(hspace=0.1)
    if (save is not None) and (plot_dir is not None):
        savefig(fig, plot_dir, save + '_energy')
    plt.show()

    # length spectrum

    fig, (ax1, ax2) = initialize_figure(figsize=(8,8))

    if plot_title is not None:
        fig.suptitle(plot_title,fontdict=font_title)

    hist_mc = mu_3d.sum(axis=(0,2))
    error_mc = np.sqrt(ssq_3d.sum(axis=(0,2)))

    hist_data = data_3d.sum(axis=(0,2))
    error_data = np.sqrt(data_3d.sum(axis=(0,2)))

    hist_ratio_base = hist_mc
    # yerr_ratio_base = yerr_base

    plot_data_hist_errorbar(
        ax1, hist_data, bins_length, error_data, label=data_label, color='k', marker='o'
    )
    plot_data_ratio_errorbar(
        ax2,
        hist_data,
        bins_length,
        error_data,
        hist_ratio_base,
        yerr_baseline=None,
        color='k', marker='o'
    )

    plot_hist_errorbar(ax1, hist_mc, bins_length, error_mc, label="MC sum")
    plot_ratio_errorbar(
        ax2,
        hist_mc,
        bins_length,
        error_mc,
        hist_ratio_base,
        yerr_baseline=None
    )
    if components is not None:
        for comp, d in components['settings'].items():
            hist = components['hists'][comp]["mu"].sum(axis=(0,2))
            yerror = np.sqrt(components['hists'][comp]["ssq"].sum(axis=(0,2)))
            plot_hist_errorbar(ax1, hist,bins_length , yerror, **d['plot_settings'])

    # setup figure
    plt.xlabel("Reconstructed Length[m]",fontdict=font_axis_label)
    ax1.set_ylabel(r"$N_{\mathrm{Events}}$",fontdict=font_axis_label)
    ax1.set_yscale("log")
    ax1.set_xscale("log")
    ax1.set_xlim(bins_length[0], bins_length[-1])
    ax1.set_ylim(1e-3,5)
    #restrict_ylim_data(ax1)
    ax2.set_ylim(ylim_length_ratio[0], ylim_length_ratio[1])
    for item in (ax1.get_xticklabels() + ax1.get_yticklabels()):
        item.set_fontsize(12)
        item.set_family('serif')
    for item in (ax2.get_xticklabels() + ax2.get_yticklabels()):
        item.set_fontsize(12)
        item.set_family('serif')
    ax1.tick_params(axis='both',which='major',width=2,length=8,direction='in')
    ax1.tick_params(axis='both',which='minor',width=1,length=5,direction='in')
    ax2.tick_params(axis='both',which='major',width=2,length=8,direction='in')
    ax2.tick_params(axis='both',which='minor',width=1,length=5,direction='in')
    # save and show plot
    ax1.legend(prop=font_legend)
    # save and show plot
    
    fig.tight_layout()
    fig.subplots_adjust(hspace=0.1)
    if (save is not None) and (plot_dir is not None):
        savefig(fig, plot_dir, save + '_zenith')
    plt.show()

    # eratio spectrum

    fig, (ax1, ax2) = initialize_figure(figsize=(8,8))

    if plot_title is not None:
        fig.suptitle(plot_title,fontdict=font_title)

    hist_mc = mu_3d.sum(axis=(1,0))
    error_mc = np.sqrt(ssq_3d.sum(axis=(1,0)))

    hist_data = data_3d.sum(axis=(1,0))
    error_data = np.sqrt(data_3d.sum(axis=(1,0)))

    hist_ratio_base = hist_mc
    # yerr_ratio_base = yerr_base

    plot_data_hist_errorbar(
        ax1, hist_data, bins_eratio, error_data, label=data_label, color='k', marker='o'
    )
    plot_data_ratio_errorbar(
        ax2,
        hist_data,
        bins_eratio,
        error_data,
        hist_ratio_base,
        yerr_baseline=None,
        color='k', marker='o'
    )

    plot_hist_errorbar(ax1, hist_mc, bins_eratio, error_mc, label="MC sum")
    plot_ratio_errorbar(
        ax2,
        hist_mc,
        bins_eratio,
        error_mc,
        hist_ratio_base,
        yerr_baseline=None
    )
    if components is not None:
        for comp, d in components['settings'].items():
            hist = components['hists'][comp]["mu"].sum(axis=(1,0))
            yerror = np.sqrt(components['hists'][comp]["ssq"].sum(axis=(1,0)))
            plot_hist_errorbar(ax1, hist,bins_eratio , yerror, **d['plot_settings'])

    # setup figure
    plt.xlabel("Energy Asymmetry",fontdict=font_axis_label)
    ax1.set_ylabel(r"$N_{\mathrm{Events}}$",fontdict=font_axis_label)
    ax1.set_yscale("log")
    ax1.set_ylim(1e-3,5)
    ax1.set_xlim(bins_eratio[-1], bins_eratio[0])
    #restrict_ylim_data(ax1)
    ax2.set_ylim(ylim_eratio_ratio[0], ylim_eratio_ratio[1])
    for item in (ax1.get_xticklabels() + ax1.get_yticklabels()):
        item.set_fontsize(12)
        item.set_family('serif')
    for item in (ax2.get_xticklabels() + ax2.get_yticklabels()):
        item.set_fontsize(12)
        item.set_family('serif')
    ax1.tick_params(axis='both',which='major',width=2,length=8,direction='in')
    ax1.tick_params(axis='both',which='minor',width=1,length=5,direction='in')
    ax2.tick_params(axis='both',which='major',width=2,length=8,direction='in')
    ax2.tick_params(axis='both',which='minor',width=1,length=5,direction='in')
    # save and show plot
    ax1.legend(prop=font_legend)
    # save and show plot
    
    fig.tight_layout()
    fig.subplots_adjust(hspace=0.1)
    if (save is not None) and (plot_dir is not None):
        savefig(fig, plot_dir, save + '_zenith')
    plt.show()

def plot_2D_DC_analysisvariables(
    mu_2d, ssq_2d, data_2d, bins_energy, bins_length,
    ylim_energy_ratio=(0, 2),
    ylim_length_ratio=(0.5, 1.5),
    components=None,
    data_label="data",
    plot_title=None,
    save=None, plot_dir=None
):
    """
    _summary_

    Parameters
    ----------
    mu_2d : np.array
        2d MC bincount
    ssq_2d : np.array
        2d expected MC fluctuations
    data_2d : np.array
        2d data bin count
    bins_energy : np.array
        energy bins
    bins_length : np.array
        length bins
    ylim_energy_ratio : tuple
        ylim in energy spectrum ratio plot
    ylim_length_ratio : tuple
        ylim in length spectrum ratio plot
    components : dict, optional
        dict(
           'settings' -> plot settings per comp,
           'hists' -> mu, ssq per comp for corresponding evaluated parameters
         ), by default None so that no individual components will be drawn
    data_label : str, optional
        label for data, by default "data"
    save : str, optional
        name for plot to save, by default None
    plot_dir : str, optional
        dir for plot to save, by default None
    """

    # energy spectrum

    fig, (ax1, ax2) = initialize_figure(figsize=(8,8))
    import matplotlib.font_manager as font_manager
    font_axis_label = {'family': 'serif',
            'color':  'black',
            'weight': 'normal',
            'size': 15,
            }
    font_title = {'family': 'serif',
            'color':  'black',
            'weight': 'bold',
            'size': 18,
            }
    font_legend = font_manager.FontProperties(family='serif',
                                       weight='normal',
                                       style='normal', size=12)
    if plot_title is not None:
        fig.suptitle(plot_title,fontdict=font_title)
    
    #
    # energy
    #
    hist_mc = mu_2d.sum(axis=1)
    error_mc = np.sqrt(ssq_2d.sum(axis=1))

    hist_data = data_2d.sum(axis=1)
    error_data = np.sqrt(data_2d.sum(axis=1))

    hist_ratio_base = hist_mc
    # yerr_ratio_base = yerr_base

    plot_data_hist_errorbar(
        ax1, hist_data, bins_energy, error_data, label=data_label, color='k', marker='o'
    )
    plot_data_ratio_errorbar(
        ax2,
        hist_data,
        bins_energy,
        error_data,
        hist_ratio_base,
        yerr_baseline=None,
        color='k', marker='o'
    )

    plot_hist_errorbar(ax1, hist_mc, bins_energy, error_mc, label="MC sum")
    plot_ratio_errorbar(
        ax2,
        hist_mc,
        bins_energy,
        error_mc,
        hist_ratio_base,
        yerr_baseline=None
    )
    if components is not None:
        for comp, d in components['settings'].items():
            hist = components['hists'][comp]["mu"].sum(axis=(1))
            yerror = np.sqrt(components['hists'][comp]["ssq"].sum(axis=(1)))
            plot_hist_errorbar(ax1, hist, bins_energy, yerror, **d['plot_settings'])

    # setup figure
    plt.xlabel("Reco energy [GeV]",fontdict=font_axis_label)
    ax1.set_ylabel(r"$N_{\mathrm{Events}}$",fontdict=font_axis_label)
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    #restrict_ylim_data(ax1)
    ax1.set_xlim(bins_energy[0], bins_energy[-1])
    ax1.set_ylim(1e-3,5)
    ax2.set_ylim(ylim_energy_ratio[0], ylim_energy_ratio[1])
    for item in (ax1.get_xticklabels() + ax1.get_yticklabels()):
        item.set_fontsize(12)
        item.set_family('serif')
    for item in (ax2.get_xticklabels() + ax2.get_yticklabels()):
        item.set_fontsize(12)
        item.set_family('serif')
    ax1.tick_params(axis='both',which='major',width=3,length=15,direction='in')
    ax1.tick_params(axis='both',which='minor',width=1,length=8,direction='in')
    ax2.tick_params(axis='both',which='major',width=3,length=15,direction='in')
    ax2.tick_params(axis='both',which='minor',width=1,length=8,direction='in')
    # save and show plot
    ax1.legend(prop=font_legend)
    fig.tight_layout()
    fig.subplots_adjust(hspace=0.1)
    if (save is not None) and (plot_dir is not None):
        savefig(fig, plot_dir, save + '_energy')
    plt.show()

    # length spectrum

    fig, (ax1, ax2) = initialize_figure(figsize=(8,8))

    if plot_title is not None:
        fig.suptitle(plot_title,fontdict=font_title)

    hist_mc = mu_2d.sum(axis=0)
    error_mc = np.sqrt(ssq_2d.sum(axis=0))

    hist_data = data_2d.sum(axis=0)
    error_data = np.sqrt(data_2d.sum(axis=0))
    hist_ratio_base = hist_mc
    # yerr_ratio_base = yerr_base

    plot_data_hist_errorbar(
        ax1, hist_data, bins_length, error_data, label=data_label, color='k', marker='o'
    )
    plot_data_ratio_errorbar(
        ax2,
        hist_data,
        bins_length,
        error_data,
        hist_ratio_base,
        yerr_baseline=None,
        color='k', marker='o'
    )

    plot_hist_errorbar(ax1, hist_mc, bins_length, error_mc, label="MC sum")
    plot_ratio_errorbar(
        ax2,
        hist_mc,
        bins_length,
        error_mc,
        hist_ratio_base,
        yerr_baseline=None
    )
    if components is not None:
        for comp, d in components['settings'].items():
            hist = components['hists'][comp]["mu"].sum(axis=(0))
            yerror = np.sqrt(components['hists'][comp]["ssq"].sum(axis=(0)))
            plot_hist_errorbar(ax1, hist,bins_length , yerror, **d['plot_settings'])

    # setup figure
    plt.xlabel("Reconstructed Length[m]",fontdict=font_axis_label)
    ax1.set_ylabel(r"$N_{\mathrm{Events}}$",fontdict=font_axis_label)
    ax1.set_yscale("log")
    ax1.set_xscale("log")
    ax1.set_xlim(bins_length[0], bins_length[-1])
    ax1.set_ylim(1e-3,5)
    #restrict_ylim_data(ax1)
    ax2.set_ylim(ylim_length_ratio[0], ylim_length_ratio[1])
    for item in (ax1.get_xticklabels() + ax1.get_yticklabels()):
        item.set_fontsize(12)
        item.set_family('serif')
    for item in (ax2.get_xticklabels() + ax2.get_yticklabels()):
        item.set_fontsize(12)
        item.set_family('serif')
    ax1.tick_params(axis='both',which='major',width=2,length=8,direction='in')
    ax1.tick_params(axis='both',which='minor',width=1,length=5,direction='in')
    ax2.tick_params(axis='both',which='major',width=2,length=8,direction='in')
    ax2.tick_params(axis='both',which='minor',width=1,length=5,direction='in')
    # save and show plot
    ax1.legend(prop=font_legend)
    # save and show plot
    
    fig.tight_layout()
    fig.subplots_adjust(hspace=0.1)
    if (save is not None) and (plot_dir is not None):
        savefig(fig, plot_dir, save + '_zenith')
    plt.show()

    



def plot_energy_and_zenith_data_MC(
    mu_2d, ssq_2d, data_2d, bins_energy, bins_zenith,
    ylim_energy_ratio=(0, 2),
    ylim_zenith_ratio=(0.5, 1.5),
    components=None,
    data_label="(pseudo) data",
    plot_title=None,
    save=None, plot_dir=None
):
    """
    _summary_

    Parameters
    ----------
    mu_2d : np.array
        2d MC bincount
    ssq_2d : np.array
        2d expected MC fluctuations
    data_2d : np.array
        2d data bin count
    bins_energy : np.array
        energy bins
    bins_zenith : np.array
        energy bins
    ylim_energy_ratio : tuple
        ylim in energy spectrum ratio plot
    ylim_zenith_ratio : tuple
        ylim in zenith spectrum ratio plot
    components : dict, optional
        dict(
           'settings' -> plot settings per comp,
           'hists' -> mu, ssq per comp for corresponding evaluated parameters
         ), by default None so that no individual components will be drawn
    data_label : str, optional
        label for data, by default "(pseudo) data"
    save : str, optional
        name for plot to save, by default None
    plot_dir : str, optional
        dir for plot to save, by default None
    """

    # energy spectrum

    fig, (ax1, ax2) = initialize_figure(figsize=(8,8))
    import matplotlib.font_manager as font_manager
    font_axis_label = {'family': 'serif',
            'color':  'black',
            'weight': 'normal',
            'size': 15,
            }
    font_title = {'family': 'serif',
            'color':  'black',
            'weight': 'bold',
            'size': 18,
            }
    font_legend = font_manager.FontProperties(family='serif',
                                       weight='normal',
                                       style='normal', size=12)
    if plot_title is not None:
        fig.suptitle(plot_title,fontdict=font_title)

    hist_mc = mu_2d.sum(axis=1)
    error_mc = np.sqrt(ssq_2d.sum(axis=1))

    hist_data = data_2d.sum(axis=1)
    error_data = np.sqrt(data_2d.sum(axis=1))

    hist_ratio_base = hist_mc
    # yerr_ratio_base = yerr_base

    plot_data_hist_errorbar(
        ax1, hist_data, bins_energy, error_data, label=data_label, color='k', marker='o'
    )
    plot_data_ratio_errorbar(
        ax2,
        hist_data,
        bins_energy,
        error_data,
        hist_ratio_base,
        yerr_baseline=None,
        color='k', marker='o'
    )

    plot_hist_errorbar(ax1, hist_mc, bins_energy, error_mc, label="MC sum")
    plot_ratio_errorbar(
        ax2,
        hist_mc,
        bins_energy,
        error_mc,
        hist_ratio_base,
        yerr_baseline=None
    )
    if components is not None:
        for comp, d in components['settings'].items():
            hist = components['hists'][comp]["mu"].sum(axis=1)
            yerror = np.sqrt(components['hists'][comp]["ssq"].sum(axis=1))
            plot_hist_errorbar(ax1, hist, bins_energy, yerror, **d['plot_settings'])

    # setup figure
    plt.xlabel("Reco energy [GeV]",fontdict=font_axis_label)
    ax1.set_ylabel(r"$N_{\mathrm{Events}}$",fontdict=font_axis_label)
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    
    
    ax1.set_xlim(bins_energy[0], bins_energy[-1])
    ax1.set_ylim(1e-3,19)
    ax2.set_ylim(ylim_energy_ratio[0], ylim_energy_ratio[1])
    for item in (ax1.get_xticklabels() + ax1.get_yticklabels()):
        item.set_fontsize(12)
        item.set_family('serif')
    for item in (ax2.get_xticklabels() + ax2.get_yticklabels()):
        item.set_fontsize(12)
        item.set_family('serif')
    ax1.tick_params(axis='both',which='major',width=3,length=15,direction='in')
    ax1.tick_params(axis='both',which='minor',width=1,length=8,direction='in')
    ax2.tick_params(axis='both',which='major',width=3,length=15,direction='in')
    ax2.tick_params(axis='both',which='minor',width=1,length=8,direction='in')
    # save and show plot
    ax1.legend(prop=font_legend)
    fig.tight_layout()
    fig.subplots_adjust(hspace=0.1)
    if (save is not None) and (plot_dir is not None):
        savefig(fig, plot_dir, save + '_energy')
    plt.show()

    # zenith spectrum

    fig, (ax1, ax2) = initialize_figure(figsize=(8,8))

    if plot_title is not None:
        fig.suptitle(plot_title,fontdict=font_title)

    hist_mc = mu_2d.sum(axis=0)
    error_mc = np.sqrt(ssq_2d.sum(axis=0))

    hist_data = data_2d.sum(axis=0)
    error_data = np.sqrt(data_2d.sum(axis=0))

    hist_ratio_base = hist_mc
    # yerr_ratio_base = yerr_base

    plot_data_hist_errorbar(
        ax1, hist_data, np.cos(bins_zenith), error_data, label=data_label, color='k', marker='o'
    )
    plot_data_ratio_errorbar(
        ax2,
        hist_data,
        np.cos(bins_zenith),
        error_data,
        hist_ratio_base,
        yerr_baseline=None,
        color='k', marker='o'
    )

    plot_hist_errorbar(ax1, hist_mc, np.cos(bins_zenith), error_mc, label="MC sum")
    plot_ratio_errorbar(
        ax2,
        hist_mc,
        np.cos(bins_zenith),
        error_mc,
        hist_ratio_base,
        yerr_baseline=None
    )
    if components is not None:
        for comp, d in components['settings'].items():
            hist = components['hists'][comp]["mu"].sum(axis=0)
            yerror = np.sqrt(components['hists'][comp]["ssq"].sum(axis=0))
            plot_hist_errorbar(ax1, hist, np.cos(bins_zenith), yerror, **d['plot_settings'])

    # setup figure
    plt.xlabel("cos(zenith)",fontdict=font_axis_label)
    ax1.set_ylabel(r"$N_{\mathrm{Events}}$",fontdict=font_axis_label)
    ax1.set_yscale("log")
    ax1.set_xlim(np.cos(bins_zenith)[-1], np.cos(bins_zenith)[0])
    ax1.set_ylim(1e-3,19)
    #restrict_ylim_data(ax1)
    ax2.set_ylim(ylim_zenith_ratio[0], ylim_zenith_ratio[1])
    for item in (ax1.get_xticklabels() + ax1.get_yticklabels()):
        item.set_fontsize(12)
        item.set_family('serif')
    for item in (ax2.get_xticklabels() + ax2.get_yticklabels()):
        item.set_fontsize(12)
        item.set_family('serif')
    ax1.tick_params(axis='both',which='major',width=2,length=8,direction='in')
    ax1.tick_params(axis='both',which='minor',width=1,length=5,direction='in')
    ax2.tick_params(axis='both',which='major',width=2,length=8,direction='in')
    ax2.tick_params(axis='both',which='minor',width=1,length=5,direction='in')
    # save and show plot
    ax1.legend(prop=font_legend)
    # save and show plot
    
    fig.tight_layout()
    fig.subplots_adjust(hspace=0.1)
    if (save is not None) and (plot_dir is not None):
        savefig(fig, plot_dir, save + '_zenith')
    plt.show()


def restrict_ylim_data(ax):
    ylim_default = ax.get_ylim()
    if ylim_default[0] < 0.1:
        # clip lower value
        ylim_to_use = (0.1, ylim_default[1])
    else:
        ylim_to_use = ylim_default
    ax.set_ylim(ylim_to_use)
