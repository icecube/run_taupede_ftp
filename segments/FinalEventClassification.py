# IceCube imports
from icecube import dataio, icetray, dataclasses, simclasses
from icecube import millipede, MuonGun
from icecube.dataclasses import I3Double, I3Particle, I3Direction, I3Position, I3VectorI3Particle, I3Constants, I3String, I3UInt64
from icecube.icetray import I3Units, I3Frame, I3ConditionalModule, traysegment
from icecube import simclasses
from icecube.hdfwriter import I3HDFWriter
from I3Tray import I3Tray
import numpy as np
from I3Tray import *



'''
These are not really 'cuts' to drop events, but rather to classify events into one of the three topologies
'''

# cut values
qtotcut = 6000         # standard HESE 6000 pe CausalQTot cut
etotcutmin = 60000     # standard HESE 60 TeV lower energy cut
lthreshold=10
etotcutmax = 10000000     # extended HESE 10 PeV upper energy cut
econfinementcut = 0.99 # corresponds to a maximum energy deposition outside the confinement region of 1%
eratiocutmin = -0.98   # corresponds to a minimum energy deposition into either cascade of 1%
eratiocutmax = 0.3     # cut based on S/B estimator with a conservative decision to keep background statistics high

def gethesemask(selfvetoarr):
    mask = (selfvetoarr == 0)
    return mask

def geteventclassmask(eventclassarr, eventclassval):
    mask = (eventclassarr == eventclassval)
    return mask

def gettopologymask(larr, lcutval=lthreshold):
    mask = (larr >= lcutval)
    return mask

def geteconfinementmask(econfinementarr, econfinementcutval=econfinementcut):
    mask = (econfinementarr >= econfinementcutval)
    return mask

def geteratiomask(eratioarr, eratiocutminval=eratiocutmin, eratiocutmaxval=eratiocutmax):
    mask = (eratioarr >= eratiocutminval) and (eratioarr <= eratiocutmaxval)
    return mask
    
    
def getsinglefilter(eventclassarr, larr, econfinementarr, eratioarr):
   
    # single stream
    singlestream = geteventclassmask(eventclassarr, 1)
   
    # double stream
    eventclassmask = geteventclassmask(eventclassarr, 2)
    nottopologymask = not(gettopologymask(larr))
    econfinementmask = geteconfinementmask(econfinementarr, econfinementcut)
    noteratiomask = not(geteratiomask(eratioarr, eratiocutmin, eratiocutmax))   
    doublestream = eventclassmask and (nottopologymask or (econfinementmask and noteratiomask))
   
    combinedmask = singlestream or doublestream
 
    return combinedmask
 
def getdoublefilter(eventclassarr, larr, econfinementarr, eratioarr):
 
    # double stream
    eventclassmask = geteventclassmask(eventclassarr, 2)
    topologymask = gettopologymask(larr)
    econfinementmask = geteconfinementmask(econfinementarr, econfinementcut)
    eratiomask = geteratiomask(eratioarr, eratiocutmin, eratiocutmax)
    doublestream = eventclassmask and topologymask and econfinementmask and eratiomask
 
    return doublestream
 
def gettrackfilter(eventclassarr, larr, econfinementarr, eratioarr):
 
    # track stream
    trackstream = geteventclassmask(eventclassarr, 3)
   
    # double stream
    eventclassmask = geteventclassmask(eventclassarr, 2)
    topologymask = gettopologymask(larr)
    noteconfinementmask = not(geteconfinementmask(econfinementarr, econfinementcut))
    doublestream = eventclassmask and topologymask and noteconfinementmask
   
    combinedmask = trackstream or doublestream
    
    return combinedmask

def checkfinaltopology(frame):
    eventclassarr = frame['HESEEventclass'].value
    larr = frame['RecoL'].value
    econfinementarr = frame['RecoEConfinement'].value
    eratioarr = frame['RecoERatio'].value
    FinalEventClass = 0
    single = getsinglefilter(eventclassarr=eventclassarr, larr=larr, econfinementarr=econfinementarr, eratioarr=eratioarr)
    double = getdoublefilter(eventclassarr=eventclassarr, larr=larr, econfinementarr=econfinementarr, eratioarr=eratioarr)
    track = gettrackfilter(eventclassarr=eventclassarr, larr=larr, econfinementarr=econfinementarr, eratioarr=eratioarr)
    
    if single and  (not(double or track)):
        FinalEventClass = 1
    elif double and (not(track or single)) :
        FinalEventClass = 2
    elif track and (not(single or double)):
        FinalEventClass = 3
    else:
        print('something is wrong, no PID assigned')

    frame['FinalTopology'] = dataclasses.I3Double(FinalEventClass)


def checkfinaltopology_evtgen(frame):
    eventclassarr = frame['HESEEventclass'].value
    larr = np.abs(frame['MyEgeneratorOutputFrameKey']['cascade_cascade_00001_distance'])
    econfinementarr = frame['RecoEConfinement'].value
    eratioarr = frame['RecoERatio'].value
    FinalEventClass = 0
    single = getsinglefilter(eventclassarr=eventclassarr, larr=larr, econfinementarr=econfinementarr, eratioarr=eratioarr)
    double = getdoublefilter(eventclassarr=eventclassarr, larr=larr, econfinementarr=econfinementarr, eratioarr=eratioarr)
    track = gettrackfilter(eventclassarr=eventclassarr, larr=larr, econfinementarr=econfinementarr, eratioarr=eratioarr)
    
    if single and  (not(double or track)):
        FinalEventClass = 1
    elif double and (not(track or single)) :
        FinalEventClass = 2
    elif track and (not(single or double)):
        FinalEventClass = 3
    else:
        print('something is wrong, no PID assigned')

    frame['FinalTopology_evtgen'] = dataclasses.I3Double(FinalEventClass)
