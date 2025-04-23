from icecube import dataio, icetray, dataclasses, simclasses
from icecube import millipede, MuonGun
from icecube.dataclasses import I3Double, I3Particle, I3Direction, I3Position, I3VectorI3Particle, I3Constants, I3String, I3UInt64
from icecube.dataclasses import I3RecoPulse, I3RecoPulseSeriesMap, I3RecoPulseSeriesMapMask
from icecube.icetray import I3Units, I3Frame, I3ConditionalModule, traysegment
from icecube import simclasses
import numpy as np
from I3Tray import *
from scipy.optimize import curve_fit
from scipy.stats import chi2
from scipy.special import erfc

'''
Whole section taken from R.Naab 
to be consistent with passing fraction calculations developed for Track sample in GlobalFit
'''
################################################################
    ############## Passing Fraction Calculations ##############
################################################################
def add_primary(frame):
    if "I3MCTree_preMuonProp" in frame:###switch to premuon prop
        def sanitize(particle):
            if particle is None:
                return dataclasses.I3Particle()
            else:
                return particle

        mcTree = frame["I3MCTree_preMuonProp"]
        primary = None
        neutrino = None
        for p in mcTree:
            if mcTree.depth(p) != 0: continue

            if p.is_neutrino:
                if neutrino is None or p.energy > neutrino.energy:
                    neutrino = p

            if primary is None or p.energy > primary.energy:
                primary = p
        if "MCPrimary" not in frame:
            frame["MCPrimary"] = sanitize(primary)



'''
Snippet to calculate the penetrating depth (to be used for passing fraction calculation later).
This is the depth at which a muon travelling colinear with the primary neutrino would have entered the detector.
Taken from R.Naab
'''
def penetrating_depth(frame, gcd, depth_name_suffix=''):
    ##  add penetrating depth dependence to the self veto probability calculation
    
    
    from icecube import MuonGun
    p = frame['MCPrimary']

    # previously used MuonGun.Cylinder(1000, 500) surface, which could result in a few events having a "nan" surface.intersection
    surface = MuonGun.Cylinder(1000, 500)
    d = surface.intersection(p.pos, p.dir)
    getDepth=p.pos + d.first*p.dir
    impactDepth = MuonGun.depth((getDepth).z)*1.e3
    frame["penetrating_depth"+depth_name_suffix+"_old"] = dataclasses.I3Double(impactDepth)

    # new calculation:
    surface_hex = MuonGun.ExtrudedPolygon.from_file(gcd)
    d = surface_hex.intersection(p.pos, p.dir)
    getDepth=p.pos + d.first*p.dir
    impactDepth = MuonGun.depth((getDepth).z)*1.e3
    frame["penetrating_depth"+depth_name_suffix] = dataclasses.I3Double(impactDepth)



import pickle
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

def PassingFraction(frame):
    import sys
    sys.path.append("/data/user/nlad/PassingFractions/")
    from spline_evaluator import Spline_Evaluator
    spline_hdl = Spline_Evaluator()
    pid = np.asarray(frame['I3MCWeightDict']['PrimaryNeutrinoType'])
    Depth = np.asarray(frame['penetrating_depth'].value)
    zenith = np.asarray(frame['I3MCWeightDict']['PrimaryNeutrinoZenith'])
    energy = np.asarray(frame['I3MCWeightDict']['PrimaryNeutrinoEnergy'])
    logofconv_pf = spline_hdl.evaluate_all(pid,Depth,zenith,energy)
    logofprompt_pf = spline_hdl.evaluate_all(pid,Depth,zenith,energy,fluxtype='pr')
   
    conv_pf=np.power(10,logofconv_pf.tolist())
    prompt_pf=np.power(10,logofprompt_pf.tolist())
    frame.Put('ConventionalAtmosphericPassingFractions', I3Double(conv_pf))
    frame.Put('PromptAtmosphericPassingFractions',I3Double(prompt_pf))
    
    return True



