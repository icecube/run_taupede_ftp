# icecube imports
from icecube import dataio, icetray, dataclasses
from icecube import phys_services, photonics_service, millipede, VHESelfVeto
from icecube.photonics_service import I3PhotoSplineService
from icecube.dataclasses import I3Double, I3Particle, I3Direction, I3Position, I3VectorI3Particle, I3Constants, I3VectorOMKey
from icecube.dataclasses import I3RecoPulse, I3RecoPulseSeriesMap, I3RecoPulseSeriesMapMask, I3TimeWindow, I3TimeWindowSeriesMap
from icecube.icetray import I3Units, I3Frame, I3ConditionalModule, traysegment
from I3Tray import I3Tray
# python system imports
import sys, os, datetime
from glob import glob
from optparse import OptionParser
import numpy as n


@traysegment
def SelfVetoWrapper(tray, name ):
    
    
    pulses = 'SplitInIcePulses'

   
    #No need for Bad String remover
    
    """
    full self veto
    """
    tray.AddModule('HomogenizedQTot', 
                   'qtot_total', 
                   Pulses=pulses,
                   Output="QTot")
    # run the veto modules
    tray.AddModule('I3LCPulseCleaning', 
                   'cleaning', 
                   OutputHLC='HLCPulses', 
                   OutputSLC='', 
                   Input=pulses)
    tray.AddModule('VHESelfVeto', 
                   'selfveto', 
                   Pulses='HLCPulses',
                   Geometry="I3Geometry",
                   OutputBool = 'VHESelfVeto')
    tray.AddModule('HomogenizedQTot', 
                   'qtot_causal', 
                   Pulses=pulses, 
                   Output='CausalQTot', 
                   VertexTime='VHESelfVetoVertexTime')
    
    
        
    
    
