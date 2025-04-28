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
from .ContainementCheck import iscontained, ispointinpolygon, getclosestdistance
# python system imports
import sys, os, datetime
from glob import glob
from optparse import OptionParser
import numpy as n

@traysegment
def MillipedeWrapper(tray, name,innerboundary, outerboundary,outeredge_x,outeredge_y,
                     Seeds=['HESEMonopodFit', 'HESETaupedeFit', 'SPEFit16'], PhotonsPerBin=5, ShowerSpacing=5, **millipede_params):
    
    """
    helper method to calculate the opening angle between two directions (input angles given in radians, return angle in degrees)
    """
    def getopeningangle(z1, a1, z2, a2):
        dot = n.sin(z1)*n.cos(a1)*n.sin(z2)*n.cos(a2)+n.sin(z1)*n.sin(a1)*n.sin(z2)*n.sin(a2)+n.cos(z1)*n.cos(z2)
        if dot < -1 : dot = -1. # treat rounding errors
        if dot > 1 : dot = 1. # treat rounding errors
        return n.arccos(dot)*180/n.pi
    
    """
    add a seed that mumillipede likes
    """
    def addseed(frame, Seeds):
        for seed in Seeds:
            seedparticle = I3Particle(frame[seed])
            # in the rare case that the seed particle is outside the outer boundary, shift it so it is inside and millipede has something ro reconstruct
            seedparticle.shape = I3Particle.InfiniteTrack
            backuppos = I3Position(seedparticle.pos)
            shiftdistance = 0.
            if seedparticle.pos.x > outerboundary:
                shiftdistance = (outerboundary+innerboundary)/(2*n.sin(seedparticle.dir.theta)*n.cos(seedparticle.dir.phi))
                backuppos.x = (outerboundary+innerboundary)/2.
            elif seedparticle.pos.x < -outerboundary:
                shiftdistance = -(outerboundary+innerboundary)/(2*n.sin(seedparticle.dir.theta)*n.cos(seedparticle.dir.phi))
                backuppos.x = -(outerboundary+innerboundary)/2.
            if seedparticle.pos.y > outerboundary:
                shiftdistance = (outerboundary+innerboundary)/(2*n.sin(seedparticle.dir.theta)*n.sin(seedparticle.dir.phi))
                backuppos.y = (outerboundary+innerboundary)/2.
            elif seedparticle.pos.y < -outerboundary:
                shiftdistance = -(outerboundary+innerboundary)/(2*n.sin(seedparticle.dir.theta)*n.sin(seedparticle.dir.phi))
                backuppos.y = -(outerboundary+innerboundary)/2.
            if seedparticle.pos.z > outerboundary:
                shiftdistance = (outerboundary+innerboundary)/(2*n.cos(seedparticle.dir.theta))
                backuppos.z = (outerboundary+innerboundary)/2.
            elif seedparticle.pos.z < -outerboundary:
                shiftdistance = -(outerboundary+innerboundary)/(2*n.cos(seedparticle.dir.theta))
                backuppos.z = -(outerboundary+innerboundary)/2.
            newpos = seedparticle.shift_along_track(shiftdistance)
            if n.abs(newpos.x) > outerboundary or n.abs(newpos.y) > outerboundary or n.abs(newpos.z) > outerboundary:
                # brute force set to intermediate boundary
                seedparticle.time = seedparticle.time-(seedparticle.pos-backuppos).magnitude/I3Constants.c
                seedparticle.pos = backuppos
            else:
                # smooth shift to intermediate boundary
                seedparticle.time = seedparticle.time+shiftdistance/I3Constants.c
                seedparticle.pos = newpos
            seedparticle.fit_status = I3Particle.OK
            seedparticle.location_type = I3Particle.InIce
            seedparticle.type = I3Particle.EMinus
            seedparticle.speed = I3Constants.c
            frame.Put(seed + 'MillipedeSeed', seedparticle)
    tray.Add(addseed, Seeds=Seeds)
    
    """
    add reconstructed energy depositions for all seeds
    """
    for seed in Seeds:
        tray.Add('MuMillipede',
            SeedTrack = seed + 'MillipedeSeed',
            Output = seed + 'MillipedeFit',
            Boundary = outerboundary,
            PhotonsPerBin = PhotonsPerBin,
            ShowerSpacing = ShowerSpacing,
            MuonSpacing=0,
            **millipede_params)

    """
    find best fit
    """
    def findbestfit(frame, Seeds, Output, Fitstatus, Eventclass):
    
        # find best fit
        fitmap = []
        for seed in Seeds:
            seedtag = seed + 'MillipedeSeed'
            fittag = seed + 'MillipedeFit'
            fitparamstag = seed + 'MillipedeFitFitParams'
            if not fittag in frame or not fitparamstag in frame:
                continue
            fitmap.append((seedtag, fittag, frame[fitparamstag].rlogl))
        dtype = [('seedtag', 'S60'), ('fittag', 'S60'), ('rlogl', float)]
        fitmap = n.array(fitmap, dtype=dtype)
        sortedfitmap = n.sort(fitmap, order='rlogl')
        bestseedtag = sortedfitmap[0]['seedtag'].decode('utf-8')
        bestfittag = sortedfitmap[0]['fittag'].decode('utf-8')
        
        if n.all(fitmap['rlogl'] == fitmap['rlogl'][0]):
            bestseedtag = Seeds[1] + 'MillipedeSeed'
            bestfittag = Seeds[1] + 'MillipedeFit'
        else:
            bestseedtag = sortedfitmap[0]['seedtag'].decode('utf-8')
            bestfittag = sortedfitmap[0]['fittag'].decode('utf-8')
        
        
        # assign correct seeds
        singleseedtag = Seeds[0]
        doubleseedtag = Seeds[1]
        trackseedtag = Seeds[2]
        
        # adapt the manual fit status depending on the opening angle between taupede and the millipede best fit
        testfit = frame[doubleseedtag + 'MillipedeFit'][0]
        bestfit = frame[bestfittag][0]
        openingangle = getopeningangle(testfit.dir.zenith, testfit.dir.azimuth, bestfit.dir.zenith, bestfit.dir.azimuth)
        print(f'Opening angle: {openingangle} for seed {bestfittag}')
        if openingangle > 30.: # limit to 30 degree opening angle
            frame.Delete(Fitstatus)
            frame.Put(Fitstatus, I3Double(I3Particle.FitStatus.InsufficientQuality))

        # add hese event class
        eventclass = 0
        print(frame[Fitstatus])
        if not Fitstatus in frame:
            frame.Put(Eventclass, I3Double(eventclass))
            print('No FitStatus In Frame')
            return
        if frame[Fitstatus].value == 0:
            eventclass = 2
        else:
            cascadefittag = singleseedtag + 'MillipedeFitFitParams'
            trackfittag = trackseedtag + 'MillipedeFitFitParams'
            if cascadefittag in frame and trackfittag in frame:
                if frame[cascadefittag].rlogl <= frame[trackfittag].rlogl:
                    eventclass = 1
                else:
                    eventclass = 3
            else:
                frame.Put(Eventclass, I3Double(eventclass))
                return
        print('FinalEvent Class',I3Double(eventclass))
        
        
        print(bestseedtag)
        
  
        # reset shape of millipede fit particle to default
        frame[bestseedtag].shape = I3Particle.ParticleShape.InfiniteTrack
        frame[bestseedtag].speed = I3Constants.c
        frame[bestseedtag].length = n.nan
        frame[bestseedtag].energy = n.nan
        frame[bestseedtag].type = I3Particle.ParticleType.unknown

        # add frame content
        frame.Put(Eventclass, I3Double(eventclass))
        frame.Rename(bestseedtag, Output)
        frame.Rename(bestfittag, Output + 'Particles')
        frame.Rename(bestfittag + 'FitParams', Output + 'FitParams')
        
        
    tray.Add(findbestfit, Seeds=Seeds, Output=name, Fitstatus='HESETaupedeFitManualFitStatus', Eventclass='HESEEventclass')

    """
    calculate truncated deposited energy
    """
    def adddepositedenergy(frame, Seed):
        print(seed)
        if not frame.Has(Seed + 'Particles'):
            print(f'frame does not contain key {Seed + "Particles"}. '
                  f'Skipping deposited energy calc.')
            return
        trackvector, trunctrackvector = frame[Seed + 'Particles'], I3VectorI3Particle()
        etot, truncetot = 0, 0
        for sec in trackvector:
            if iscontained(sec.pos.x, sec.pos.y, sec.pos.z, innerboundary,outeredge_x, outeredge_y): # soft containment
                truncetot += sec.energy
                trunctrackvector.append(sec)
            etot += sec.energy
        frame.Put(Seed + 'TruncatedParticles', trunctrackvector)
        frame.Put(Seed + 'TruncatedDepositedEnergy', I3Double(truncetot))
        frame.Put(Seed + 'DepositedEnergy', I3Double(etot))
    tray.Add(adddepositedenergy, Seed=name)

    """
    clean up
    """
    deletekeys = []
    for seed in Seeds:
        deletekeys += [seed + 'MillipedeSeed', seed + 'MillipedeFit', seed + 'MillipedeFitFitParams']
    tray.Add('Delete', keys=deletekeys)
