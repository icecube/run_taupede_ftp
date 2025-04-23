# IceCube imports
from icecube import dataio, icetray, dataclasses
from icecube import millipede, MuonGun
from icecube.dataclasses import I3Double, I3Particle, I3Direction, I3Position, I3VectorI3Particle, I3Constants
from icecube.dataclasses import I3RecoPulse, I3RecoPulseSeriesMap, I3RecoPulseSeriesMapMask
from icecube.icetray import I3Units, I3Frame, I3ConditionalModule, traysegment
from I3Tray import I3Tray
from .ContainementCheck import iscontained, ispointinpolygon, getclosestdistance
# Python system imports
import sys, datetime, os
from glob import glob
from optparse import OptionParser
import numpy as n



"""
calculation of the energy ratio
"""
def getenergyratio(double_particles):
    energythreshold=1e3
    
    cascade1, cascade2 = double_particles
    eratio = (cascade1.energy-cascade2.energy)/(cascade1.energy+cascade2.energy) if cascade1.energy+cascade2.energy >= energythreshold else 1.
    return eratio

"""
calculation of the energy confinement
"""
def getenergyconfinement(double_particles, track_particles, key):
    energythreshold=1e3
    distancethreshold = 40
    # the double cascade
    cascade1, cascade2 = double_particles

    # due to a falsely configured muon spacing in millipede, need to split the cascade and track segments and then rebin
    if key == 'reco':
        particle_buffer = []
        track_dx_buffer, track_dy_buffer, track_dz_buffer, track_de_buffer = [], [], [], []
        for particle in track_particles:
            if particle.shape == I3Particle.Cascade:
                track_dx_buffer.append(particle.pos.x)
                track_dy_buffer.append(particle.pos.y)
                track_dz_buffer.append(particle.pos.z)
                track_de_buffer.append(particle.energy)
            else:
                particle_buffer.append(particle)
        track_dx_buffer = n.asarray(track_dx_buffer)
        track_dy_buffer = n.asarray(track_dy_buffer)
        track_dz_buffer = n.asarray(track_dz_buffer)
        track_de_buffer = n.asarray(track_de_buffer)
        for particle in particle_buffer:
            argmin = n.argmin(n.sqrt((track_dx_buffer - particle.pos.x)**2
                                   + (track_dy_buffer - particle.pos.y)**2
                                   + (track_dz_buffer - particle.pos.z)**2))
            track_de_buffer[argmin] += particle.energy
    elif key == 'true':
        track_dx_buffer, track_dy_buffer, track_dz_buffer, track_de_buffer = [], [], [], []
        for particle in track_particles:
            track_dx_buffer.append(particle.pos.x)
            track_dy_buffer.append(particle.pos.y)
            track_dz_buffer.append(particle.pos.z)
            track_de_buffer.append(particle.energy)
        track_dx_buffer = n.asarray(track_dx_buffer)
        track_dy_buffer = n.asarray(track_dy_buffer)
        track_dz_buffer = n.asarray(track_dz_buffer)
        track_de_buffer = n.asarray(track_de_buffer)

    # find energy depositions in the vicinity of the double cascade vertices
    mask1 = ( n.sqrt((track_dx_buffer-cascade1.pos.x)**2
                   + (track_dy_buffer-cascade1.pos.y)**2
                   + (track_dz_buffer-cascade1.pos.z)**2) < distancethreshold )
    mask2 = ( n.sqrt((track_dx_buffer-cascade2.pos.x)**2
                   + (track_dy_buffer-cascade2.pos.y)**2
                   + (track_dz_buffer-cascade2.pos.z)**2) < distancethreshold )
    # calculate the energy confinement
    esum = n.sum(track_de_buffer)
    edep = n.sum(track_de_buffer[mask1 | mask2])
    econfinement = edep/esum if esum >= energythreshold else 0.
    return econfinement


"""
calculate reco observables
"""

def calculaterecoobservables(frame,innerboundary,outeredge_x, outeredge_y):
    energythreshold=1e3
    # the single vertex containment
    cascade = frame['HESEMonopodFit']
    contained = iscontained(cascade.pos.x, cascade.pos.y, cascade.pos.z, innerboundary,outeredge_x, outeredge_y)

    # the double vertex properties and containment
    cascade1, cascade2 = frame['HESETaupedeFitParticles']
    contained1 = iscontained(cascade1.pos.x, cascade1.pos.y, cascade1.pos.z, innerboundary,outeredge_x, outeredge_y)
    contained2 = iscontained(cascade2.pos.x, cascade2.pos.y, cascade2.pos.z, innerboundary,outeredge_x, outeredge_y)
    length = n.log10((cascade1.pos-cascade2.pos).magnitude)
    reco_l = (cascade1.pos-cascade2.pos).magnitude
    eratio = getenergyratio([cascade1, cascade2])
    e1 = n.log10(cascade1.energy)
    e2 = n.log10(cascade2.energy)
    
    # the direction is taken from the HESEMillipedeFit (i.e. best fit out of three hypotheses)
    zenith = frame['HESEMillipedeFit'].dir.zenith
    azimuth = frame['HESEMillipedeFit'].dir.azimuth
    
    # if truncated deposited energy is below threshold switch over to non-truncated energy deposition
    etot = n.log10(frame['HESEMillipedeFitTruncatedDepositedEnergy'].value)
    reco_e = frame['HESEMillipedeFitTruncatedDepositedEnergy'].value
    if frame.Has('energy_reco'):
	      frame.Delete('energy_reco')
    usetruncated = True
    if etot < n.log10(energythreshold):
        etot = n.log10(frame['HESEMillipedeFitDepositedEnergy'].value)
        reco_e = (frame['HESEMillipedeFitDepositedEnergy'].value)
        usetruncated = False
        
    # the calculation of the energy confinement
    if usetruncated:
        econfinement = getenergyconfinement([cascade1, cascade2], frame['HESEMillipedeFitTruncatedParticles'], 'reco')
    else:
        econfinement = getenergyconfinement([cascade1, cascade2], frame['HESEMillipedeFitParticles'], 'reco')
    # sanitize values
    if not n.isfinite(length):
        length = 0 # limit to 1 m
        reco_l = 1
    if not n.isfinite(e1):
        e1 = 0 # limit to 1 GeV
    if not n.isfinite(e2):
        e2 = 0 # limit to 1 GeV
    if not n.isfinite(etot):
        etot = 0 # limit to 1 GeV
        reco_e = 1


    # limit length to > 1m and < ~2km and smear the boundary region
    if length <= 0:
        length = n.minimum(0.0 + n.abs(n.random.normal(loc=0., scale=0.3)),0.99) 
        reco_l = 10**length
    elif length >= 3.3:
         length = n.maximum(3.3 - n.abs(n.random.normal(loc=0., scale=0.05)),3.2)
         reco_l = 10**length

    # save it to the frame
    frame.Put('RecoContainedSingle', icetray.I3Bool(contained))
    frame.Put('RecoContained1', icetray.I3Bool(contained1))
    frame.Put('RecoContained2', icetray.I3Bool(contained2))
    frame.Put('RecoLogL', I3Double(length))
    frame.Put('RecoERatio', I3Double(eratio))
    frame.Put('RecoEConfinement', I3Double(econfinement))
    frame.Put('RecoLogE1', I3Double(e1))
    frame.Put('RecoLogE2', I3Double(e2))
    frame.Put('RecoLogETot', I3Double(etot))
    frame.Put('RecoZenith', I3Double(zenith))
    frame.Put('RecoAzimuth', I3Double(azimuth))
    frame.Put('HESETaupedeFit1', I3Particle(cascade1))
    frame.Put('HESETaupedeFit2', I3Particle(cascade2))
    frame.Put('RecoETot', I3Double(reco_e))
    frame.Put('RecoL', I3Double(reco_l))
    frame.Put('HESEMonopodFit_x', I3Double(cascade.pos.x))
    frame.Put('HESEMonopodFit_y', I3Double(cascade.pos.y))
    frame.Put('HESEMonopodFit_z', I3Double(cascade.pos.z))

    frame.Put('HESETaupedeFit1_x', I3Double(cascade1.pos.x))
    frame.Put('HESETaupedeFit1_y', I3Double(cascade1.pos.y))
    frame.Put('HESETaupedeFit1_z', I3Double(cascade1.pos.z))

    frame.Put('HESETaupedeFit2_x', I3Double(cascade2.pos.x))
    frame.Put('HESETaupedeFit2_y', I3Double(cascade2.pos.y))
    frame.Put('HESETaupedeFit2_z', I3Double(cascade2.pos.z))

    # add fit params information
    for key in ['HESEMonopodFit', 'HESETaupedeFit', 'HESEMillipedeFit']:
        if key in frame:
            frame.Put(key + 'Logl', I3Double(frame[key + 'FitParams'].logl))
            frame.Put(key + 'LoglNdof', I3Double(frame[key + 'FitParams'].ndof))
            frame.Put(key + 'Chi2', I3Double(frame[key + 'FitParams'].chi_squared))
            frame.Put(key + 'Chi2Ndof', I3Double(frame[key + 'FitParams'].chi_squared_dof))
            
            
    return True
