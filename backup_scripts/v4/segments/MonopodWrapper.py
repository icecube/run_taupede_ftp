# icecube imports
from icecube import dataio, icetray, dataclasses
from icecube import phys_services, photonics_service, millipede, VHESelfVeto
from icecube.photonics_service import I3PhotoSplineService
from icecube.dataclasses import I3Double, I3Particle, I3Direction, I3Position, I3VectorI3Particle, I3Constants, I3VectorOMKey
from icecube.dataclasses import I3RecoPulse, I3RecoPulseSeriesMap, I3RecoPulseSeriesMapMask, I3TimeWindow, I3TimeWindowSeriesMap
from icecube.icetray import I3Units, I3Frame, I3ConditionalModule, traysegment
from icecube.millipede import MonopodFit, MuMillipedeFit, TaupedeFit, HighEnergyExclusions, MillipedeFitParams
#from icecube.simprod.segments.PropagateMuons import PropagateMuons
from I3Tray import I3Tray

# python system imports
import sys, os, datetime
from glob import glob
from optparse import OptionParser
import numpy as n




@traysegment
def MonopodWrapper(tray, name, Seed='CombinedCascadeSeed_L3', Iterations=4, PhotonsPerBin=15, **millipede_params):

    """
    the one iteration amplitude monopod fit
    """
    tray.Add(MonopodFit, name + 'Seed',
        Seed = Seed,
        Iterations = 1,
        PhotonsPerBin = -1,
        **millipede_params)
        
    """
    the multiple iteration timed monopod fit
    """
    tray.Add(MonopodFit, name,
        Seed = name + 'Seed',
        Iterations = Iterations,
        PhotonsPerBin = PhotonsPerBin,
        **millipede_params)
        
    """
    clean up
    """
    tray.Add('Delete', keys=[Seed, name + 'Seed', name + 'SeedFitParams'])
