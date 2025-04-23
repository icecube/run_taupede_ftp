# icecube imports
from icecube import dataio, icetray, dataclasses
from icecube import phys_services, photonics_service, millipede, VHESelfVeto
from icecube.photonics_service import I3PhotoSplineService
from icecube.dataclasses import I3Double, I3Particle, I3Direction, I3Position, I3VectorI3Particle, I3Constants, I3VectorOMKey
from icecube.dataclasses import I3RecoPulse, I3RecoPulseSeriesMap, I3RecoPulseSeriesMapMask, I3TimeWindow, I3TimeWindowSeriesMap
from icecube.icetray import I3Units, I3Frame, I3ConditionalModule, traysegment

from I3Tray import I3Tray
from icecube import STTools, DomTools
from icecube.STTools.seededRT.configuration_services import I3DOMLinkSeededRTConfigurationService
from icecube import linefit, lilliput, clast, cscd_llh
from icecube.lilliput.segments import I3SinglePandelFitter, I3IterativePandelFitter

# python system imports
import sys, os, datetime
from glob import glob
from optparse import OptionParser
import numpy as n

# for level 2
from icecube import STTools, DomTools
from icecube.STTools.seededRT.configuration_services import I3DOMLinkSeededRTConfigurationService
from icecube import linefit, lilliput, clast, cscd_llh


@traysegment
def Level2ReconstructionWrapper(tray, name, Pulses='SplitInIcePulses'):

    """
    SRT cleaning: standard level 2
    """
    seededRTConfig = I3DOMLinkSeededRTConfigurationService(
        ic_ic_RTRadius              = 150.0*I3Units.m,
        ic_ic_RTTime                = 1000.0*I3Units.ns,
        treat_string_36_as_deepcore = False,
        useDustlayerCorrection      = False,
        allowSelfCoincidence        = True)
    tray.Add('I3SeededRTCleaning_RecoPulseMask_Module', 'seededrt',
        InputHitSeriesMapName  = Pulses,
        OutputHitSeriesMapName = 'SRT' + Pulses,
        STConfigService        = seededRTConfig,
        SeedProcedure          = 'HLCCoreHits',
        NHitsThreshold         = 2,
        MaxNIterations         = 3,
        Streams                = [I3Frame.Physics])

    """
    offline muon reconstruction of LineFit and SPEFit2: taken from standard level 2 filterscripts
    """
    tray.Add(linefit.simple, 'LineFit',
        inputResponse = 'SRT' + Pulses,
        fitName = 'LineFit',
        If = lambda frame : not 'LineFit' in frame)
    tray.Add(I3SinglePandelFitter, 'SPEFitSingle',
        pulses = 'SRT' + Pulses,
        fitname="SPEFitSingle",
        seeds = ['LineFit'],
        If = lambda frame : not 'SPEFit2' in frame)
    tray.Add(I3IterativePandelFitter, 'SPEFit2',
        pulses = 'SRT' + Pulses,
        n_iterations = 2,
        fitname="SPE2Fit",
        seeds = ['SPEFitSingle'],
        If = lambda frame : not 'SPEFit2' in frame)

    """
    offline cascade reconstruction of CascadeLast and CascadeLlhVertexFit: taken from standard level 2 filterscripts
    """
    tray.Add('I3CLastModule', 'CascadeLast_L2',
        Name = 'CascadeLast_L2',
        InputReadout = Pulses,
        If = lambda frame : not 'CascadeLast_L2' in frame)
    tray.Add('I3CscdLlhModule', 'CascadeLlhVertexFit_L2',
        InputType = 'RecoPulse', # ! Use reco pulses
        RecoSeries = Pulses, # ! Name of input pulse series
        FirstLE = True, # Default
        SeedWithOrigin = False, # Default
        SeedKey = 'CascadeLast_L2', # ! Seed fit - CLast reco
        MinHits = 8, # ! Require 8 hits
        AmpWeightPower = 0.0, # Default
        ResultName = 'CascadeLlhVertexFit_L2', # ! Name of fit result
        Minimizer = 'Powell', # ! Set the minimizer to use
        PDF = 'UPandel', # ! Set the pdf to use
        ParamT = '1.0, 0.0, 0.0, false',   # ! Setup parameters
        ParamX = '1.0, 0.0, 0.0, false',   # ! Setup parameters
        ParamY = '1.0, 0.0, 0.0, false',   # ! Setup parameters
        ParamZ = '1.0, 0.0, 0.0, false',   # ! Setup parameters
        If = lambda frame : not 'CascadeLlhVertexFit_L2' in frame)

    """
    clean up
    """
    deletekeys = ['LineFit_HuberFit', 'LineFit_Pulses_delay_cleaned', 'LineFit_debiasedPulses', 'LineFit_linefit_final_rusage',
                  'SPEFitSingle', 'SPEFitSingleFitParams']
    tray.Add('Delete', Keys=deletekeys)
