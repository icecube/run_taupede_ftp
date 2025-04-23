# IceCube imports
from icecube import dataio, icetray, dataclasses
from icecube import millipede, MuonGun
from icecube.dataclasses import I3Double, I3Particle, I3Direction, I3Position, I3VectorI3Particle, I3Constants
from icecube.dataclasses import I3RecoPulse, I3RecoPulseSeriesMap, I3RecoPulseSeriesMapMask
from icecube.icetray import I3Units, I3Frame, I3ConditionalModule, traysegment
from I3Tray import I3Tray

# Python system imports
import sys, datetime, os
from glob import glob
from optparse import OptionParser
import numpy as n


"""
calculation of the energy ratio
"""

def getenergyratio(double_particles):
    cascade1, cascade2 = double_particles
    eratio = (cascade1.energy-cascade2.energy)/(cascade1.energy+cascade2.energy) if cascade1.energy+cascade2.energy >= energythreshold else 1.
    return eratio

  
"""
calculation of the energy confinement
"""
def getenergyconfinement(double_particles, track_particles, key):

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
 
