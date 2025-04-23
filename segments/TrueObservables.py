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
calculate true observables
"""
def calculatetrueobservables(frame,innerboundary,outeredge_x, outeredge_y,name_suffix=''):
    energythreshold=1e3
    # the single vertex containment
    if not frame.Has('MostEnergeticCascade'):
        mcTree = frame["I3MCTree"]
        cascade = dataclasses.get_most_energetic_cascade(mcTree)
        print('woops, most energetic cascade wasn`t added')
    else:
        cascade = frame['MostEnergeticCascade']
    contained = iscontained(cascade.pos.x, cascade.pos.y, cascade.pos.z, innerboundary,outeredge_x, outeredge_y)

    # the double vertex properties and containment
    cascade1, cascade2 = frame['MCInteractionDoubleBang'+name_suffix]
    contained1 = iscontained(cascade1.pos.x, cascade1.pos.y, cascade1.pos.z, innerboundary,outeredge_x, outeredge_y)
    contained2 = iscontained(cascade2.pos.x, cascade2.pos.y, cascade2.pos.z, innerboundary,outeredge_x, outeredge_y)
    length = n.log10((cascade1.pos-cascade2.pos).magnitude)
    true_l = ((cascade1.pos-cascade2.pos).magnitude)
    eratio = getenergyratio([cascade1, cascade2])
    e1 = n.log10(cascade1.energy)
    e2 = n.log10(cascade2.energy)
    true_e1 = cascade1.energy
    true_e2 = cascade2.energy
    zenith = cascade1.dir.zenith
    azimuth = cascade1.dir.azimuth
    
    # if truncated deposited energy is below threshold switch over to non-truncated energy deposition
    etot = n.log10(frame['MCTrueTruncatedDepositedEnergy'+name_suffix].value)
    true_etot = (frame['MCTrueTruncatedDepositedEnergy'+name_suffix].value)
    usetruncated = True
    if etot < n.log10(energythreshold):
        etot = n.log10(frame['MCTrueDepositedEnergy'+name_suffix].value)
        true_etot = (frame['MCTrueDepositedEnergy'+name_suffix].value)
        usetruncated = False
    # the calculation of the energy confinement
    if usetruncated:
        econfinement = getenergyconfinement([cascade1, cascade2], frame['MCTrueTruncateddEdx'+name_suffix], 'true')
    else:
        econfinement = getenergyconfinement([cascade1, cascade2], frame['MCTruedEdx'+name_suffix], 'true')
    # sanitize values
    if not n.isfinite(length):
        length = 0 # limit to 1 m
        true_l = 1
    if not n.isfinite(e1):
        e1 = 0 # limit to 1 GeV
    if not n.isfinite(e2):
        e2 = 0 # limit to 1 GeV
    if not n.isfinite(etot):
        etot = 0 # limit to 1 GeV
        true_etot = 1
    
    # limit length to > 1m and < ~2km and smear the boundary region
    if length <= 0:
        length = n.minimum(0.0 + n.abs(n.random.normal(loc=0., scale=0.3)), 0.99)
        true_l = 10**length
    elif length >= 3.3:
        length = n.maximum(3.3 - n.abs(n.random.normal(loc=0., scale=0.05)), 3.2)
        true_l = 10**length


    # save it to the frame
    frame.Put('TrueContainedSingle'+name_suffix, icetray.I3Bool(contained))
    frame.Put('TrueContained1'+name_suffix, icetray.I3Bool(contained1))
    frame.Put('TrueContained2'+name_suffix, icetray.I3Bool(contained2))
    frame.Put('TrueLogL'+name_suffix, I3Double(length))
    frame.Put('TrueERatio'+name_suffix, I3Double(eratio))
    frame.Put('TrueEConfinement'+name_suffix, I3Double(econfinement))
    frame.Put('TrueLogE1'+name_suffix, I3Double(e1))
    frame.Put('TrueLogE2'+name_suffix, I3Double(e2))
    frame.Put('TrueLogETot'+name_suffix, I3Double(etot))
    frame.Put('TrueZenith'+name_suffix, I3Double(zenith))
    frame.Put('TrueAzimuth'+name_suffix, I3Double(azimuth))
    frame.Put('TrueE1'+name_suffix,I3Double(true_e1))
    frame.Put('TrueE2'+name_suffix,I3Double(true_e2))
    frame.Put('TrueETot'+name_suffix,I3Double(true_etot))
    frame.Put('TrueL'+name_suffix,I3Double(true_l))
    frame.Put('MCInteractionDoubleBang1'+name_suffix, I3Particle(cascade1))
    frame.Put('MCInteractionDoubleBang2'+name_suffix, I3Particle(cascade2))
    return True
    
    
