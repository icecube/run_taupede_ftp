from I3Tray import *
from icecube import icetray, phys_services, dataio, dataclasses,MuonGun
from icecube.simprod import segments
from icecube.simclasses import I3MMCTrack
from icecube import MuonGun, simclasses
from icecube.icetray import traysegment, I3Module
import sys, math
import copy
import glob
import numpy as np
from icecube.icetray import traysegment, I3Module
import os
@traysegment
def Get_MuonWeight(tray,name,infile,flux_model,prefix,nfile):
    
    def harvest_generators(infile):
        """
        Harvest serialized generator configurations from a set of I3 files.
        """
        from icecube.icetray.i3logging import log_info as log
        generator = None
        print('nfile is {0}'.format(nfile))
        #for fname in infile:
        for fname in infile:

            
                f = dataio.I3File(fname)
                print(fname)
                
                fr = f.pop_frame(icetray.I3Frame.Stream('S'))
                f.close()
                if fr is not None:
                    
                    for k in fr.keys():
                        v = fr[k]
                        if isinstance(v, MuonGun.GenerationProbability):
                            log('%s: found "%s" (%s)' % (f, k, type(v).__name__), unit="MuonGun")
                            
                            if generator is None:
                                generator = v*nfile
                            else:
                                generator += v*nfile

        return generator
      
    model = MuonGun.load_model(flux_model)
    generator = harvest_generators(infile)
    
    tray.AddModule('I3MuonGun::WeightCalculatorModule', 'MuonWeight'+prefix, Model=model,
        Generator=generator)
