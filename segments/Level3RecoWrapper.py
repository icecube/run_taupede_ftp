# icecube imports
from icecube import dataio, icetray, dataclasses
from icecube import phys_services, photonics_service, millipede, VHESelfVeto
from icecube.photonics_service import I3PhotoSplineService
from icecube.dataclasses import I3Double, I3Particle, I3Direction, I3Position, I3VectorI3Particle, I3Constants, I3VectorOMKey
from icecube.dataclasses import I3RecoPulse, I3RecoPulseSeriesMap, I3RecoPulseSeriesMapMask, I3TimeWindow, I3TimeWindowSeriesMap
from icecube.icetray import I3Units, I3Frame, I3ConditionalModule, traysegment
from I3Tray import I3Tray
from icecube.level3_filter_cascade.level3_Recos import CascadeLlhVertexFit, SPEFit
from icecube import STTools, DomTools
from icecube.STTools.seededRT.configuration_services import I3DOMLinkSeededRTConfigurationService

@traysegment
def Level3ReconstructionWrapper(tray, name,pulses):
    
    tray.Add('I3OMSelection<I3RecoPulseSeries>', 'omselection',
        InputResponse = 'SRT' + pulses,
        OmittedStrings = [79,80,81,82,83,84,85,86],
        OutputOMSelection = 'SRT' + pulses + '_BadOMSelectionString',
        OutputResponse = 'SRT' + pulses + '_IC_Singles')
    
    
    """
    CascadeLlhVertexFit: standard level 3 fit on SRT pulses without DeepCore
    """
    tray.Add(CascadeLlhVertexFit, 'CascadeLlhVertexFit_L3',
        Pulses = 'SRT' + pulses + '_IC_Singles')
    
    """
    SPEFit: standard level 3 fit on SRT pulses without DeepCore (first guesses are SPEFit2 and LineFit from level 2)
    """
    tray.Add(SPEFit, 'SPEFit16',
        Pulses = 'SRT' + pulses + '_IC_Singles',
        Iterations = 16)

    
    
    """
    make the cascade level 3 seed: take the best combination out of all level 2 and level 3 fits to build a seed
    """
    def addlevel3seed(frame, Output):
    
        # the seed particle
        seed = I3Particle()
        seed.pos = I3Position(0, 0, 0)
        seed.dir = I3Direction(0, 0)
        seed.time = 0
        seed.energy = 0.
        seed.length = 0.
        seed.speed = I3Constants.c
        seed.fit_status = I3Particle.OK
        seed.shape = I3Particle.Cascade
        seed.location_type = I3Particle.InIce
        
        # possible solutions (ordered, construct seed in any case, even if level 2 + 3 recos failed)
        vertexfits = ['CascadeLlhVertexFit_L3', 'CascadeLlhVertexFit_L2', 'CascadeLast_L2']
        directionfits = ['CscdL3_SPEFit16', 'SPEFit2', 'LineFit']
        
        # vertex + time
        for vertexfitname in vertexfits:
            if not vertexfitname in frame:
                continue
            vertexfit = frame[vertexfitname]
            if vertexfit.fit_status == I3Particle.OK and vertexfit.pos.r >= 0 and vertexfit.time >= 0:
                seed.pos = vertexfit.pos
                seed.time = vertexfit.time
                break
        
        # direction
        for directionfitname in directionfits:
            if not directionfitname in frame:
                continue
            directionfit = frame[directionfitname]
            if directionfit.fit_status == I3Particle.OK and directionfit.dir.theta >= 0 and directionfit.dir.phi >= 0:
                seed.dir = directionfit.dir
                break
        
        # save it
        frame.Put(Output, seed)

    tray.Add(addlevel3seed, Output=name)
    


