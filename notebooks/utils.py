import numpy as np
import matplotlib as mat
import matplotlib.font_manager as font_manager
import matplotlib.pyplot as plt

font_axis_label = {'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': 26,
        }
font_title = {'family': 'serif',
        'color':  'black',
        'weight': 'bold',
        'size': 35,
        }
font_legend = font_manager.FontProperties(family='serif',
                                   weight='normal',
                                   style='normal', size=17)

def get_variables( hdf_file ):

    variables = {}

    variables['PrimaryNeutrinoEnergy'] = hdf_file['I3MCWeightDict']['PrimaryNeutrinoEnergy'].values
    variables['PrimaryNeutrinoType'] = hdf_file['I3MCWeightDict']['PrimaryNeutrinoType'].values
    variables['PrimaryNeutrinoAzimuth'] = hdf_file['I3MCWeightDict']['PrimaryNeutrinoAzimuth'].values
    variables['PrimaryNeutrinoZenith'] = hdf_file['I3MCWeightDict']['PrimaryNeutrinoZenith'].values
    variables['InteractionType'] = hdf_file['I3MCWeightDict']['InteractionType'].values

    variables['TrueETot'] = hdf_file['TrueETot'].value.values
    variables['TrueL'] = hdf_file['TrueL'].value.values
    variables['TrueZenith'] = hdf_file['TrueZenith'].value.values
    variables['MCInteractionEventclass'] = hdf_file['MCInteractionEventclass'].value.values

    variables['RecoAzimuth'] = hdf_file['RecoAzimuth'].value.values
    variables['RecoEConfinement'] = hdf_file['RecoEConfinement'].value.values
    variables['RecoERatio'] = hdf_file['RecoERatio'].value.values
    variables['RecoZenith'] = hdf_file['RecoZenith'].value.values
    variables['RecoL'] = hdf_file['RecoL'].value.values
    variables['RecoETot'] = hdf_file['RecoETot'].value.values
    variables['FinalTopology'] = hdf_file['FinalTopology'].value.values
    variables['FinalEventClass'] = hdf_file['FinalEventClass'].value.values
    variables['ConventionalAtmosphericPassingFractions'] = hdf_file['ConventionalAtmosphericPassingFractions'].value.values
    variables['PromptAtmosphericPassingFractions'] = hdf_file['PromptAtmosphericPassingFractions'].value.values

    # # evt gen
    # variables['RecoL_evtgen'] = np.abs(hdf_file['MyEgeneratorOutputFrameKey']['cascade_cascade_00001_distance'].values)
    # variables['RecoL_evtgen_noabs'] = hdf_file['MyEgeneratorOutputFrameKey']['cascade_cascade_00001_distance'].values
    # variables['RecoE1_evtgen'] = hdf_file['MyEgeneratorOutputFrameKey']['cascade_energy'].values 
    # variables['RecoE2_evtgen'] = hdf_file['MyEgeneratorOutputFrameKey']['cascade_cascade_00001_energy'].values 
    # variables['RecoETot_evtgen'] = variables['RecoE1_evtgen'] + variables['RecoE2_evtgen']
    # variables['FinalTopology_evtgen'] = hdf_file['FinalTopology_evtgen'].value.values
    # variables['FinalEventClass_evtgen'] = hdf_file['FinalEventClass_evtgen'].value.values

    return variables

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


def plot_2dHist(x,y,weights,xbins_start,xbins_stop,xbins,ybins_start,ybins_stop,ybins,\
                xlogspace,ylogspace,title,eventcount,ax):

    
    if xlogspace:
        x_bins=np.logspace(xbins_start,xbins_stop,xbins)
        
    else:
        x_bins=np.linspace(xbins_start,xbins_stop,xbins)
        
    
    if ylogspace:
        
        y_bins =np.logspace(ybins_start,ybins_stop,ybins)
    else:
        
        y_bins =np.linspace(ybins_start,ybins_stop,ybins)
    bins = [x_bins,y_bins]
    
    H, xedges, yedges = np.histogram2d(x,y,bins = [x_bins,y_bins],\
                                   weights=weights)
    H /= np.sum(H)
    
    
    #norm=colors.LogNorm()
    h = ax.pcolormesh(xedges, yedges, H.T,norm = mat.colors.LogNorm())
    #h = ax.pcolormesh(xedges, yedges, H.T,norm = mat.colors.Normalize(vmin=vmin, vmax=vmax))
    #h = ax.pcolormesh(xedges, yedges, H.T,norm = mat.colors.Normalize())
    ax.set_xlim(min(x_bins),max(x_bins))
    ax.set_ylim(min(y_bins),max(y_bins))
    ax.set_title(title,fontdict=font_axis_label)
    
    
    
    if xlogspace:
        ax.set_xscale('log')
        
        
    if ylogspace:
        ax.set_yscale('log')
        
        
        
    if eventcount:
        EventCount, _, __ = np.histogram2d(x,y,bins = [x_bins,y_bins])
        print("Total Hist Count is %d"%np.sum(EventCount))
        if xlogspace:
                midbins_x = np.sqrt(x_bins[:-1] * x_bins[1:])
                
        else:
                midbins_x = (x_bins[1:] - x_bins[:-1])/2 + x_bins[:-1]   
                
                
        if ylogspace:
                
                midbins_y = np.sqrt(y_bins[:-1] * y_bins[1:])
        else:
                
                midbins_y = (y_bins[1:] - y_bins[:-1])/2 + y_bins[:-1]   
                
        for i_x in range(len(xedges)-1):
            for i_y  in range(len(yedges)-1):
                ax.text(midbins_x[i_x],
                        midbins_y[i_y],
                        int(EventCount[i_x, i_y]),
                        color='w',
                        ha='center',
                        va='center',
                        fontweight='normal',
                        fontsize = 12,
                       )
   
    for item in (ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(16)
        item.set_family('serif')
    
    clb = plt.colorbar(h)
    
    clb.ax.tick_params(labelsize=16)

    return clb
    
    
    