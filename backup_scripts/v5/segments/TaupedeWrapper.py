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
def TaupedeWrapper(tray, name,innerboundary, outerboundary,outeredge_x,outeredge_y,Seed='HESEMonopodFit', Iterations=4, PhotonsPerBin=15,**millipede_params):
    ethreshold = 1e3 
    """
    determine a good fit using the particle properties instead of the rlogl
    """
    def isgoodfit(cascade1, cascade2):
        goodfit = True
        if cascade1.energy < ethreshold or cascade2.energy < ethreshold:
            goodfit = False
        if n.abs(cascade1.pos.x) > outerboundary or n.abs(cascade1.pos.y) > outerboundary or n.abs(cascade1.pos.z) > outerboundary:
            goodfit = False
        if n.abs(cascade2.pos.x) > outerboundary or n.abs(cascade2.pos.y) > outerboundary or n.abs(cascade2.pos.z) > outerboundary:
            goodfit = False
        return goodfit

    scales = [10., 25., 50., 100.]
    deletekeys = []

    """
    add seeds
    """
    def addbruteforceseeds(frame, Seed, Output):
        taupedeseedbase = I3Particle(frame[Seed])
        taupedeseedbase.speed = I3Constants.c
        taupedeseedbase.fit_status = I3Particle.OK
        taupedeseedbase.shape = I3Particle.InfiniteTrack
        taupedeseedbase.location_type = I3Particle.InIce
        for scale in scales:

            # center seed
            taupedeseedcenter = I3Particle(taupedeseedbase)
            taupedeseedcenter.length = scale
            frame.Put(Output+'_scale%i_center' % scale, taupedeseedcenter)

            # forward seed
            taupedeseedforward = I3Particle(taupedeseedbase)
            taupedeseedforward.pos = taupedeseedbase.shift_along_track(scale)
            taupedeseedforward.time = taupedeseedbase.time + scale/I3Constants.c
            taupedeseedforward.length = scale
            frame.Put(Output+'_scale%i_forward' % scale, taupedeseedforward)
    
            # backward seed
            taupedeseedbackward = I3Particle(taupedeseedbase)
            taupedeseedbackward.pos = taupedeseedbase.shift_along_track(-scale)
            taupedeseedbackward.time = taupedeseedbase.time - scale/I3Constants.c
            taupedeseedbackward.length = scale
            frame.Put(Output+'_scale%i_backward' % scale, taupedeseedbackward)
    tray.Add(addbruteforceseeds, Seed=Seed, Output='TaupedeBruteForceWrapperSeed')

    for scale in scales:
        for tag in ['forward', 'backward', 'center']:
            seedtag = 'TaupedeBruteForceWrapperSeed_scale%i_%s' % (scale, tag)
            recotag = 'TaupedeBruteForceWrapperReco_scale%i_%s' % (scale, tag)
            if tag == 'forward':
                lengthBounds = [0, 800]
            elif tag == 'backward':
                lengthBounds = [-800, 0]
            elif tag == 'center':
                lengthBounds = [-800, 800]
            tray.Add(TaupedeFit, recotag,
                Seed = seedtag,
                PhotonsPerBin = -1,
                LengthBounds = lengthBounds,
                StepL = 1,
                StepT = 0,
                Iterations = 4,
                **millipede_params)
            deletekeys.append(seedtag)
            deletekeys.append(recotag)
            deletekeys.append(recotag+'Particles')
            deletekeys.append(recotag+'FitParams')

    """
    find the three best fits as new seeds
    """
    def findbestseeds(frame, Seed):
        fitmap = []
        for scale in scales:
            for tag in ['forward', 'backward', 'center']:
                fittag = '%s_scale%i_%s' % (Seed, scale, tag)
                fitparticlestag = '%s_scale%i_%sParticles' % (Seed, scale, tag)
                fitparamstag  = '%s_scale%i_%sFitParams' % (Seed, scale, tag)
                if not fitparamstag in frame:
                    continue
                rlogl = frame[fitparamstag].rlogl
                goodfit = isgoodfit(frame[fitparticlestag][0], frame[fitparticlestag][1])
                fitmap.append((fittag, rlogl, goodfit))

        dtype = [('fittag', 'S60'), ('rlogl', float), ('goodfit', bool)]
        fitmap = n.array(fitmap, dtype=dtype)
        fitmask = fitmap['goodfit']
        bestfittags = [] # save three best results, preferably with goodfit tag otherwise global minimum rlogl
        for fitentry in n.sort(fitmap[fitmask], order='rlogl'):
            if len(bestfittags) < 3:
                bestfittags.append(fitentry['fittag'])
        for fitentry in n.sort(fitmap[n.invert(fitmask)], order='rlogl'):
            if len(bestfittags) < 3:
                bestfittags.append(fitentry['fittag'])
        for i in range(len(bestfittags)):
            bestfittag = bestfittags[i].decode("utf-8") 
            if frame[bestfittag+'Particles'][0].time < frame[bestfittag+'Particles'][1].time:
                bestseed = I3Particle(frame[bestfittag])
            else:
                bestseed = I3Particle(frame[bestfittag+'Particles'][1])
                bestseed.length = n.abs(frame[bestfittag].length)
                bestseed.shape = I3Particle.InfiniteTrack
                bestseed.fit_status = I3Particle.OK
                bestseed.location_type = I3Particle.InIce
            frame.Put(Seed+str(i), bestseed)
    tray.Add(findbestseeds, Seed='TaupedeBruteForceWrapperReco')

    for i in range(3):
        tray.Add(TaupedeFit, name+str(i),
            Seed='TaupedeBruteForceWrapperReco'+str(i),
            LengthBounds = [-800, 800],
            StepL = 1,
            Iterations = Iterations,
            PhotonsPerBin = PhotonsPerBin,
            **millipede_params)
        deletekeys.append('TaupedeBruteForceWrapperReco' + str(i))
        deletekeys.append('TaupedeBruteForceWrapperReco' + str(i) + 'Particles')
        deletekeys.append('TaupedeBruteForceWrapperReco' + str(i) + 'FitParams')

    """
    find the best fit
    """
    def findbestfit(frame, Seed, Output):
        fitmap = []
        for i in range(3):
            fittag = Seed + str(i)
            fitparticlestag = Seed + str(i) + 'Particles'
            fitparamstag = Seed + str(i) + 'FitParams'
            if not fitparamstag in frame:
                continue
            rlogl = frame[fitparamstag].rlogl
            goodfit = isgoodfit(frame[fitparticlestag][0], frame[fitparticlestag][1])
            fitmap.append((fittag, rlogl, goodfit))
        dtype = [('fittag', 'S60'), ('rlogl', float), ('goodfit', bool)]
        fitmap = n.array(fitmap, dtype=dtype)
        fitmask = fitmap['goodfit']
        filteredfitmap = fitmap[fitmask]
        if len(filteredfitmap) > 0:
            sortedfitmap = n.sort(filteredfitmap, order='rlogl')
        else:
            sortedfitmap = n.sort(fitmap, order='rlogl')
        bestfittag = sortedfitmap[0]['fittag'].decode('utf-8')
        if frame[bestfittag+'Particles'][0].time < frame[bestfittag+'Particles'][1].time:
            bestfit = I3Particle(frame[bestfittag])
            bestfitparticles = I3VectorI3Particle([frame[bestfittag+'Particles'][0], frame[bestfittag+'Particles'][1]])
        else:
            bestfit = I3Particle(frame[bestfittag+'Particles'][1])
            bestfit.length = n.abs(frame[bestfittag].length)
            bestfit.shape = I3Particle.InfiniteTrack
            bestfit.fit_status = I3Particle.OK
            bestfit.location_type = I3Particle.InIce
            cascade1 = I3Particle(frame[bestfittag+'Particles'][1])
            cascade2 = I3Particle(frame[bestfittag+'Particles'][0])
            cascade1.length = bestfit.length
            cascade2.length = 0.
            cascade1.fit_status = I3Particle.OK
            cascade2.fit_status = I3Particle.NotSet
            cascade1.location_type = I3Particle.InIce
            cascade2.location_type = I3Particle.Anywhere
            bestfitparticles = I3VectorI3Particle([cascade1, cascade2])

        frame.Put(Output, bestfit)
        frame.Put(Output + 'Particles', bestfitparticles)
        frame.Put(Output + 'FitParams', frame[bestfittag + 'FitParams'])
    tray.Add(findbestfit, Seed=name, Output=name)

    """
    add manual fit status information of the seed
    """
    def addfitstatus(frame, Seed):

        fitstatus = I3Particle.FitStatus.OK
        outtag = Seed + 'ManualFitStatus'

        # check if fit result in frame
        if Seed not in frame or Seed + 'Particles' not in frame:
            fitstatus = I3Particle.FitStatus.MissingSeed
            frame.Put(outtag, I3Double(fitstatus))
            return True

        # check for general weirdness in the fit result
        fitstatus = float(frame[Seed].fit_status)
        cascade1, cascade2 = frame[Seed + 'Particles']
        if fitstatus != 0 or (cascade1.energy + cascade2.energy == 0) or cascade1.pos.r < 0 or cascade2.pos.r < 0:
            fitstatus = I3Particle.FitStatus.FailedToConverge
            frame.Put(outtag, I3Double(fitstatus))
            return True

        # check for soft containment and energy threshold (E1, E2 > 1TeV)
        goodfit = isgoodfit(cascade1, cascade2)
        contained1 = iscontained(cascade1.pos.x, cascade1.pos.y, cascade1.pos.z, 
                                 innerboundary,outeredge_x, outeredge_y)
        contained2 = iscontained(cascade2.pos.x, cascade2.pos.y, cascade2.pos.z, 
                                 innerboundary,outeredge_x, outeredge_y)
        if not goodfit or not contained1 or not contained2:
            fitstatus = I3Particle.FitStatus.InsufficientQuality
       
        frame.Put(outtag, I3Double(fitstatus))
    tray.Add(addfitstatus, Seed=name)

    """
    clean up
    """
    for i in range(3):
        deletekeys.append(name + str(i))
        deletekeys.append(name + str(i) + 'Particles')
        deletekeys.append(name + str(i) + 'FitParams')
    tray.Add('Delete', Keys=deletekeys)
