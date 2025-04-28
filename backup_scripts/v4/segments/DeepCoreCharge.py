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

DeepCore_Strings = [79,80,81,82,83,84,85,86]

def ChargeDepositions(frame):
    DeepCore_Charge_map = {}
    TotalCharge_map = {}

    DC_charge = []
    Total_Charge = []

    pulsemap = I3RecoPulseSeriesMap.from_frame(frame, "SplitInIcePulses")
    for pmt in pulsemap.keys():
        if pmt.string in DeepCore_Strings:
            DC_charge.append(sum([pulse.charge for pulse in pulsemap[pmt]]))
            DeepCore_Charge_map[pmt] = sum([pulse.charge for pulse in pulsemap[pmt]])
        Total_Charge.append(sum([pulse.charge for pulse in pulsemap[pmt]]))
        TotalCharge_map[pmt]= sum([pulse.charge for pulse in pulsemap[pmt]])
    Q_DeepCore = sum(DC_charge)
    Q_Total = sum(Total_Charge)

    #frame['DeepCore_Charge_Map'] = dataclasses.I3MapStringDouble(DeepCore_Charge_map)
    #frame['TotalCharge_Map'] = dataclasses.I3MapStringDouble(TotalCharge_map)
    frame['DeepCore_Charge'] = dataclasses.I3Double(Q_DeepCore)
    frame['TotalCharge'] = dataclasses.I3Double(Q_Total)
                                                
    return True                                           
    