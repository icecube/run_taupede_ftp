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
def MCInfoWrapper(tray, name, Pulses='SplitInIcePulses'):
    """
    add useful monte carlo information
    """
    #tray.Add('TauGeneratesMuon', If=lambda f : 'I3MCTree' in f)

    #tray.Add('VertexInFiducialVolume', If=lambda f : 'I3MCTree' in f)

    def VertexIsFarOutsideTheDetector(frame):
        """
        Stores an I3Bool, containing True if the MC vertex is from far outside the detector
        """
        if not 'I3MCTree' in frame:
	        return True
        extra_distance = 50.*I3Units.m

        # cylinder boundary values for IC86:
        minz = -513.*I3Units.m - extra_distance
        maxz = 525.*I3Units.m + extra_distance
        maxrho = 602.*I3Units.m + extra_distance

        vertexPos = frame['VertexPosition'] # from 'VertexInFiducialVolume'
        z = vertexPos.z
        rho = n.sqrt(vertexPos.x**2 + vertexPos.y**2)

        if z < minz or z > maxz or rho > maxrho:
            isOutside = True
        else:
            isOutside = False
        frame['VertexIsFarOutsideTheDetector'] = icetray.I3Bool(isOutside)
   	
    tray.Add(VertexIsFarOutsideTheDetector)


    def GetNeutrino(frame):
        if not 'I3MCTree' in frame:
            return True
        def sanitize(particle):
            if particle is None:
                return I3Particle()
            else:
                return particle
        mcTree = frame['I3MCTree']
        primary = None
        neutrino = None
        for p in mcTree:
            if mcTree.depth(p) != 0: continue
            if p.is_neutrino:
                if neutrino is None or p.energy > neutrino.energy:
                    neutrino = p
            if primary is None or p.energy > primary.energy:
                primary = p
        frame['MostEnergeticPrimary'] = sanitize(primary)
        frame['MostEnergeticNeutrino'] = sanitize(neutrino)
        frame['MostEnergeticInIce'] = sanitize(dataclasses.get_most_energetic_inice(mcTree))
        frame['MostEnergeticMuon'] = sanitize(dataclasses.get_most_energetic_muon(mcTree))
        frame['MostEnergeticCascade'] = sanitize(dataclasses.get_most_energetic_cascade(mcTree))
    tray.Add(GetNeutrino, Streams=[I3Frame.Physics])

    def GetMCTrack(frame):
        if not 'I3MCTree' in frame:
            return True
        mcTree = frame['I3MCTree']
        if frame.Has("MCTrack"):
            frame.Delete("MCTrack")
        trackParticle = None
        cascadeParticle = None
        numCascades = 0
        neutrino = None
        for p in mcTree:
            depth = mcTree.depth(p)
            if depth == 0:
                if neutrino is None or p.energy > neutrino.energy:
                    neutrino = p
            if depth != 1: continue # depth==0 is the root (assume it is the primary neutrino)
            if p.type in [I3Particle.ParticleType.MuPlus,
                        I3Particle.ParticleType.MuMinus,
                        I3Particle.ParticleType.TauPlus,
                        I3Particle.ParticleType.TauMinus]:
                if trackParticle is not None:
                    continue # ignore multiple leptons (CORSIKA)
					#raise RuntimeError("got multiple leptons from a single neutrino.")
                trackParticle = p
            else:
                if cascadeParticle is None or p.energy > cascadeParticle.energy:
                    cascadeParticle = p
                numCascades += 1
        theTrack = None
        if trackParticle is not None:
            theTrack = trackParticle
        else:
            if numCascades == 0: theTrack = None
            if numCascades == 1: theTrack = cascadeParticle
            if neutrino is None:
                raise RuntimeError("Internal error. Cascades found, but no neutrino in MCTree.")
            theTrack = neutrino
        if theTrack is None:
            raise RuntimeError("no MC track could be found in MCTree")

        # shift the vertex to the point of closest approach to the origin (0,0,0)
        a = - (theTrack.pos.x*theTrack.dir.x + theTrack.pos.y*theTrack.dir.y + theTrack.pos.z*theTrack.dir.z)
        newPos = I3Position(theTrack.pos.x + theTrack.dir.x * a,
                            theTrack.pos.y + theTrack.dir.y * a,
                            theTrack.pos.z + theTrack.dir.z * a)
        newTime = theTrack.time + a/I3Constants.c

        # generate a "reconstructed" particle from the MCTrack
        outputTrack = I3Particle()  
        outputTrack.shape = I3Particle.ParticleShape.InfiniteTrack
        outputTrack.pos = newPos
        outputTrack.dir = theTrack.dir
        outputTrack.time = newTime
        outputTrack.fit_status = I3Particle.FitStatus.OK
        outputTrack.location_type = I3Particle.LocationType.InIce

        frame['MCTrack'] = outputTrack
    tray.Add(GetMCTrack, Streams=[I3Frame.Physics])

    """
    find detector interaction point
    """
    def _FindDetectorVolumeIntersections(frame, recoParticle, geometry):
        intersectionPoints = VHESelfVeto.IntersectionsWithInstrumentedVolume(geometry, recoParticle)
        intersectionTimes = []
        for intersectionPoint in intersectionPoints:
            vecX = intersectionPoint.x - recoParticle.pos.x
            vecY = intersectionPoint.y - recoParticle.pos.y
            vecZ = intersectionPoint.z - recoParticle.pos.z
            prod = vecX*recoParticle.dir.x + vecY*recoParticle.dir.y + vecZ*recoParticle.dir.z
            dist = n.sqrt(vecX**2 + vecY**2 + vecZ**2)
            if prod < 0.: dist *= -1.
            if abs(prod-dist) > 1e-3*icetray.I3Units.m:
                raise RuntimeError("intersection points are not on track")
            intersectionTimes.append(dist/I3Constants.c + recoParticle.time)
        min_index = 0
        if len(intersectionTimes)==0:
            frame["IntersectionPoint"] = recoParticle.pos
        else:
            for i in range(len(intersectionTimes)):
            # trying to catch special cases (not intersecting, only one intersection...)
                if i==0:
                    min_time = intersectionTimes[i]
                else:
                    if intersectionTimes[i]<min_time:
                        min_index = i
                        min_time = intersectionTimes[i]
            frame["IntersectionPoint"] = I3Position(intersectionPoints[min_index])
        sortedTimes = sorted(intersectionTimes)
        return sortedTimes

    def FindDetectorVolumeIntersections(frame, TrackName='', OutputTimeWindow=None, TimePadding=0.):
        if frame.Has('IntersectionPoint'):
            return True
        if not TrackName in frame:
            return True
        if OutputTimeWindow is not None:
            twName = OutputTimeWindow
        else:
            twName = TrackName + "TimeRange"

        theTrack = frame[TrackName]
        geometry = frame["I3Geometry"]

        times = _FindDetectorVolumeIntersections(frame, theTrack, geometry)

        # trying to catch special cases
        if len(times) == 0:
             # raise RuntimeError("track does not intersect the detector volume")
             frame[twName] = dataclasses.I3TimeWindow()
        elif len(times) == 1:
             raise RuntimeError("tracks with only one intersection are not supported")
        else:
             tWindow = dataclasses.I3TimeWindow(times[0]-TimePadding, times[-1]+TimePadding)
             frame[twName] = tWindow

    tray.Add(FindDetectorVolumeIntersections,
            TimePadding = 60.*I3Units.m/I3Constants.c,
            TrackName = 'MCTrack',
            OutputTimeWindow = 'ContainedTimeRange')



    def rawweightcalculator(frame):
        interactiontype = frame['I3MCWeightDict']['InteractionType']
        timerange = frame['ContainedTimeRange']
        if not frame.Has('MCTrack'):
            mctrack = get_mctrack(frame)
        else:
            mctrack = frame['MCTrack']
        interactiondepth = mctrack.shift_time_track(timerange.start-mctrack.time).z # calc depth of neutrino intersection

        frame.Put('MCInteractionType', I3Double(interactiontype))
        frame.Put('MCInteractionDepth', I3Double(interactiondepth))

    tray.Add(rawweightcalculator)
