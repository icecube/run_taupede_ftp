import simweights
import numpy as np
import pickle


spline_file = '/data/ana/Diffuse/NNMFit/MCEq_splines/v1.2.1/MCEq_splines_PRI-Gaisser-H4a_INT-SIBYLL23c_allfluxes.pickle'

def Append_Weights(file):
    
    hdf_file = file['hdf_file']
    nfiles   = file['nfiles']
    
    # conventional            
    flux_keys_conv =  ['conv_antinumu','conv_numu','conv_antinue','conv_nue','conv_antinutau','conv_nutau']
    spline_object_conv = SplineHandler(spline_file, flux_keys_conv)
    conv_flux = spline_object_conv.return_weight
    
    generator_conv = lambda pdgid, energy, cos_zen: conv_flux(pdgid, energy, cos_zen)
    
    weighter = simweights.NuGenWeighter(hdf_file,nfiles=nfiles)
    file['variables']['Weights_Conventional'] = weighter.get_weights(generator_conv)
    file['variables']['Weights_Conventional_PassedVeto'] = file['variables']['Weights_Conventional']*file['variables']['ConventionalAtmosphericPassingFractions']
            
    # prompt
    flux_keys_pr =  ['pr_antinumu','pr_numu','pr_antinue','pr_nue','pr_antinutau','pr_nutau']
    spline_object_pr = SplineHandler(spline_file, flux_keys_pr)
    pr_flux = spline_object_pr.return_weight
    
    generator_pr = lambda pdgid, energy, cos_zen: pr_flux(pdgid, energy, cos_zen)
    
    weighter = simweights.NuGenWeighter(hdf_file,nfiles=nfiles)
    file['variables']['Weights_Prompt'] = weighter.get_weights(generator_pr)
    file['variables']['Weights_Prompt_PassedVeto'] = file['variables']['Weights_Prompt']*file['variables']['PromptAtmosphericPassingFractions']

    # combine atmospheric
    file['variables']['Weights_Atmospheric'] = file['variables']['Weights_Conventional_PassedVeto'] + file['variables']['Weights_Prompt_PassedVeto']

    # astro
    def AstroFluxModel(pdgid, energy, cos_zen):
        flux = 0.5*(2.12*1e-18)*(energy/1e5)**-2.87
        return flux

    weighter = simweights.NuGenWeighter(hdf_file,nfiles=nfiles)
    file['variables']['Weights_Astro'] = weighter.get_weights(AstroFluxModel)

    return file




class SplineHandler(object):
    """
    Class implementing the flux weight calculation from a
    spline file created before
    (adjusted for MCEq)
    """
    IS_SYS_PARAM = False
    def __init__(self, spline_file, flux_keys, barr_key=None):
        self.spline_file = spline_file
        self.flux_keys = flux_keys
        self.barr_key = barr_key
        self.Ecut = 5e8 ## force highE weights to zero
        if self.barr_key is not None:
            self.spline_dict = self._load_pickle(spline_file)[1][self.barr_key]
            self.mag = 0
            self.spline_in_log = False
        else:
            self.spline_dict = self._load_pickle(spline_file)
            self.mag = self.spline_dict[0]["Emultiplier"]
            self.spline_in_log = True
        self._pid_dict = {"conv_numu" : 14,
                          "conv_antinumu" : -14,
                          "conv_nue" : 12,
                          "conv_antinue" : -12,
                          "conv_nutau" : 16,
                          "conv_antinutau" : -16,
                          "k_numu" : 14,
                          "k_antinumu" : -14,
                          "k_nue" : 12,
                          "k_antinue" : -12,
                          "pi_numu" : 14,
                          "pi_antinumu" : -14,
                          "pi_nue" : 12,
                          "pi_antinue" : -12,
                          "numu" : 14,
                          "antinumu" : -14,
                          "nue" : 12,
                          "antinue" : -12,
                          "conv_numuMSIS00_ICSouthPoleJanuary": 14,
                          "conv_antinumuMSIS00_ICSouthPoleJanuary": -14,
                          "conv_numuMSIS00_ICSouthPoleJuly": 14,
                          "conv_antinumuMSIS00_ICSouthPoleJuly": -14,
                          "pr_antinumu" : -14,
                          "pr_numu" : 14,
                          "pr_nue" : 12,
                          "pr_antinue" : -12,
                          "pr_antinutau" : -16,
                          "pr_nutau": 16}
    
    def resolve_pid(self, flux_key):
        
        return self._pid_dict[flux_key]
    
    def _load_pickle(self, pickle_file):
        """
        Returns the content of a pickle file.
        Compatible with python2 AND python3 pickle files.
        """
        try:
            with open(pickle_file, 'rb') as f:
                pickle_data = pickle.load(f)
        except UnicodeDecodeError as e:
            with open(pickle_file, 'rb') as f:
                pickle_data = pickle.load(f, encoding='latin1')
        except Exception as e:
            print('Unable to load data ', pickle_file, ':', e)
            raise
        return pickle_data
    
    def return_weight(self, pid_ints, energys, cosZs):
        """
        Return weight from spline. Correct for the E**mag factor that was
        applied during creationp.
        Args: _particleID, coszenith, energy
        """
        theta_deg = 180./np.pi*np.arccos(cosZs)
        logenergy = np.log10(energys)
        weights = np.zeros_like(cosZs)
        #logger.debug("Calculating MCEq weights from spline %s",
        #             self.spline_file)
        
        for flux_key in self.flux_keys:
            pid_idcs = np.argwhere(pid_ints == self.resolve_pid(flux_key))
            if self.barr_key is None:
                weights[pid_idcs] = 10**self.spline_dict[1][flux_key](theta_deg[pid_idcs],
                                                                      logenergy[pid_idcs],
                                                                      grid=False)
            else:
                #special treatment for barr-splines (were built slightly diff.)
                weights[pid_idcs] = self.spline_dict[flux_key](
                    theta_deg[pid_idcs],
                    logenergy[pid_idcs],
                    grid=False)
            ##hard fix to remove the highE madness of MCEq gradients
            #logger.warning("Forcing {} atmospheric weights for super highE weights to zero for numerical stability.".format(flux_key))
            weights[np.argwhere(logenergy>np.log10(self.Ecut))] = 0.
            ## check for NaN
            ##logger.warning("Found {} events with NaN weights.".format(len(weights[np.argwhere(np.isnan(weights))])))
            #weights[np.argwhere(np.isnan(weights))] = 0.
            ##correct for the E**mag factor from MCEq
            weights[pid_idcs] /= energys[pid_idcs]**self.mag
        return weights



def onedimension_hist(ax,x,weights,xlog,ylog,bins_start,bins_stop,bins,error,color,**kwargs):
    
    
    if xlog:
        bins=np.logspace(bins_start,bins_stop,bins)
        
    else:
         bins=np.linspace(bins_start,bins_stop,bins)
        
    counts, edges = np.histogram(x,weights=weights,bins=bins)
    print('hist count is',np.sum(counts))
    ax.plot(edges,np.append(counts,counts[-1]),
                 drawstyle="steps-post",color=color,
                 **kwargs)
    
    if error:
        
        error,bin_centres = error_cal(edges,weights,x)
        
        ax.errorbar(x=bin_centres, y=counts,
                 yerr=error, color=color,fmt='o', markersize=8,capsize=5)
    if xlog:
        plt.semilogx()
    if ylog:
        plt.semilogy()
    
