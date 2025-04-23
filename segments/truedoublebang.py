
# icecube imports
from icecube import dataio, icetray, dataclasses, simclasses
from icecube import phys_services, photonics_service, millipede, VHESelfVeto
from icecube.photonics_service import I3PhotoSplineService
from icecube.dataclasses import I3Double, I3Particle, I3Direction, I3Position, I3VectorI3Particle, I3Constants, I3VectorOMKey
from icecube.dataclasses import I3RecoPulse, I3RecoPulseSeriesMap, I3RecoPulseSeriesMapMask, I3TimeWindow, I3TimeWindowSeriesMap
from icecube.icetray import I3Units, I3Frame, I3ConditionalModule, traysegment
from icecube.millipede import MonopodFit, MuMillipedeFit, TaupedeFit, HighEnergyExclusions, MillipedeFitParams
from .ContainementCheck import iscontained, ispointinpolygon, getclosestdistance
#from icecube.simprod.segments.PropagateMuons import PropagateMuons
from I3Tray import I3Tray
import numpy as np
from icecube.sim_services import ShowerParameters



invisibletypes = [I3Particle.NuE, I3Particle.NuEBar,
                     I3Particle.NuMu, I3Particle.NuMuBar,
                     I3Particle.NuTau, I3Particle.NuTauBar,
                    I3Particle.MuPlus,I3Particle.MuMinus,
                    I3Particle.TauPlus,I3Particle.TauMinus]

def showerparameters(cascade,mctree):
        energy=list()
        shower_max = list()
        particle = None
        if cascade.type in invisibletypes:
            for particle in mctree.get_daughters(cascade):
                if (np.round((particle.pos-cascade.pos).magnitude, 2) ==
                        np.round(cascade.length, 2)):
                    if particle.type not in invisibletypes:
                       if ShowerParameters(particle.type,particle.energy).b==0:
                          energy.append(0)
                          shower_max.append(0)
                       else:
                            energy.append(particle.energy*\
                                      ((ShowerParameters(particle.type,particle.energy).emScale)))
                            shower_max.append(((ShowerParameters(particle.type,particle.energy).a\
                                      -1)/ShowerParameters(particle.type,particle.energy).b))
                    
        shower_max = np.array(shower_max)
        energy = np.array(energy)
        ShowerMaxCombined = (np.dot(shower_max,energy)/np.sum(energy))
        e  = np.sum(energy)

        return ShowerMaxCombined, e, particle
    
"""
reset correct particle IDs to refer from VertexOutgoingLepton
and VertexOutgoingHadron to the I3MCTree
"""
def flattentree(mctree, parent, children):
    daughters = mctree.get_daughters(parent)
    for daughter in daughters:
        children.append(daughter)
        if len(mctree.get_daughters(daughter)) > 0:
            flattentree(mctree, daughter, children)

"""
get flattened mctree containing particles that deposit light in the detector
"""
def flattensecondaries(mctree, particle, secondaries):
    daughters = mctree.get_daughters(particle)
    for d in daughters:
        if len(mctree.get_daughters(d)) == 0:

             if  ((d.shape!=I3Particle.Dark) and (d.type not in invisibletypes)):
                  secondaries.append(d)
        else:
                  flattensecondaries(mctree, d, secondaries)
                
             

def addtruedoublebang(frame,
                          outeredge_x,
                          outeredge_y,
                     innerboundary,outerboundary):
        # get frame objects
        mctree = frame['I3MCTree']
        eventclass = frame['MCInteractionEventclass']
        
        if ('VertexOutgoingLepton' in frame and
                        'VertexOutgoingHadron' in frame):
             print('theres lepton and hadron present')
        elif ('VertexOutgoingHadronOne' in frame and
                        'VertexOutgoingHadronTwo' in frame):
               print('both hadrons')
        else:
               print('I dont know whats happening, help!')
        # single cascade
        if eventclass == 1:
            if ('VertexOutgoingLepton' in frame and
                        'VertexOutgoingHadron' in frame):
                    hadron = frame['VertexOutgoingHadron']
                    lepton = frame['VertexOutgoingLepton']

                    # the single cascade
                    showermax1,energy1,particle = showerparameters(hadron,mctree)
                    showermax2,energy2,particle = showerparameters(lepton,mctree)
                    cascade1 = I3Particle(hadron)
                    cascade1.energy = energy1 + energy2
                    cascade1.shape = I3Particle.InfiniteTrack
                    showermax = ((showermax1 * energy1 + showermax2 * energy2) /
                                 (energy1 + energy2))
                    cascade1.pos = cascade1.shift_along_track(showermax)
                    cascade1.time = hadron.time + showermax/I3Constants.c
                    cascade1.length = 0.
                    cascade1.speed = 0.
                    cascade1.shape = I3Particle.Cascade
                    cascade1.location_type = I3Particle.InIce

                    # the second cascade is empty
                    cascade2 = I3Particle()
                    cascade2.pos = cascade1.pos
                    cascade2.dir = cascade1.dir
                    cascade2.time = cascade1.time
                    cascade2.energy = 0.
                    cascade2.length = 0.
                    cascade2.speed = 0.
                    cascade2.shape = I3Particle.Cascade
                    cascade2.location_type = I3Particle.InIce

            elif ('VertexOutgoingHadronOne' in frame and
                        'VertexOutgoingHadronTwo' in frame):
                    hadron1 = frame['VertexOutgoingHadronOne']
                    hadron2 = frame['VertexOutgoingHadronTwo']

                    # the single cascade
                    showermax1,energy1,particle = showerparameters(hadron1,mctree)
                    showermax2,energy2,particle = showerparameters(hadron2,mctree)
                    cascade1 = I3Particle(hadron1)
                    cascade1.energy = energy1 + energy2
                    cascade1.shape = I3Particle.InfiniteTrack
                    showermax = ((showermax1 * energy1 + showermax2 * energy2) /
                                 (energy1 + energy2))
                    cascade1.pos = cascade1.shift_along_track(showermax)
                    cascade1.time = hadron1.time + showermax/I3Constants.c
                    cascade1.length = 0.
                    cascade1.speed = 0.
                    cascade1.shape = I3Particle.Cascade
                    cascade1.location_type = I3Particle.InIce

                    # the second cascade is empty
                    cascade2 = I3Particle()
                    cascade2.pos = cascade1.pos
                    cascade2.dir = cascade1.dir
                    cascade2.time = cascade1.time
                    cascade2.energy = 0.
                    cascade2.length = 0.
                    cascade2.speed = 0.
                    cascade2.shape = I3Particle.Cascade
                    cascade2.location_type = I3Particle.InIce

            else:
                    cascade=None

                    if 'VertexOutgoingHadron' in frame:
                        cascade = frame['VertexOutgoingHadron']
                        
                    elif 'VertexOutgoingLepton' in frame:
                        cascade = frame['VertexOutgoingLepton']
                        

                    if cascade.type in invisibletypes:
                        showermax,e,particle = showerparameters(cascade,mctree)
                        cascade1 = I3Particle(particle)
                        cascade1.energy = e
                        cascade1.shape = I3Particle.InfiniteTrack
                        cascade1.pos = cascade1.shift_along_track(showermax)
                        cascade1.time = particle.time + showermax/I3Constants.c
                        cascade1.length = 0.
                        cascade1.speed = 0.
                        cascade1.shape = I3Particle.Cascade
                        cascade1.location_type = I3Particle.InIce
                    else:
                        cascade1 = I3Particle(cascade)
                        cascade1.energy = cascade.energy*\
                                          ((ShowerParameters(cascade.type,cascade.energy).emScale))
                        cascade1.shape = I3Particle.InfiniteTrack
                        showermax = ((ShowerParameters(cascade.type,cascade.energy).a-1)\
                                    /ShowerParameters(cascade.type,cascade.energy).b)
                        cascade1.pos = cascade1.shift_along_track(showermax)
                        cascade1.time = cascade.time + showermax/I3Constants.c
                        cascade1.length = 0.
                        cascade1.speed = 0.
                        cascade1.shape = I3Particle.Cascade
                        cascade1.location_type = I3Particle.InIce
            
                    

                    # the second cascade is empty
                    cascade2 = I3Particle()
                    cascade2.pos = cascade1.pos
                    cascade2.dir = cascade1.dir
                    cascade2.time = cascade1.time
                    cascade2.energy = 0.
                    cascade2.length = 0.
                    cascade2.speed = 0.
                    cascade2.shape = I3Particle.Cascade
                    cascade2.location_type = I3Particle.InIce

            mcdoublebang = I3VectorI3Particle([cascade1, cascade2])

            vertexshift1 = 0.
            vertexshift2 = 0.

        # double cascade
        elif eventclass == 2:
            hadron = frame['VertexOutgoingHadron']
            lepton = frame['VertexOutgoingLepton']

            # the hadronic cascade
            cascade1 = I3Particle(hadron)
            cascade1.energy = hadron.energy*((ShowerParameters(hadron.type,hadron.energy).emScale))
            ShowerMax = ((ShowerParameters(hadron.type,hadron.energy).a-1)\
                             /ShowerParameters(hadron.type,hadron.energy).b)
            cascade1.pos = cascade1.shift_along_track(ShowerMax)
            cascade1.time = hadron.time + ShowerMax/I3Constants.c
            cascade1.length = 0.
            cascade1.speed = 0.
            cascade1.shape = I3Particle.Cascade
            cascade1.location_type = I3Particle.InIce

            # the decay cascade
            showermax,e,particle = showerparameters(lepton,mctree)
            cascade2 = I3Particle(particle)
            cascade2.energy = e
            cascade2.shape = I3Particle.InfiniteTrack
            cascade2.pos = cascade2.shift_along_track(showermax)
            cascade2.time = particle.time + showermax/I3Constants.c
            cascade2.length = 0.
            cascade2.speed = 0.
            cascade2.shape = I3Particle.Cascade
            cascade2.location_type = I3Particle.InIce
            # since we shift to the shower maximum it can happen
            # that the cascades are switched in time but we keep the ordering
            # constant to not confuse the energy ratio calculation
            # of which the interaction and decay cascades are
            mcdoublebang = I3VectorI3Particle([cascade1, cascade2])
            vertexshift1 = 0.
            vertexshift2 = 0.

        # track
        elif eventclass == 3:
            if 'VertexOutgoingLepton' in frame: # NuGen
                leptons = [frame['VertexOutgoingLepton']]
            else: # MuonGun or CORSIKA (not well defined)
                leptons = mctree.get_daughters(mctree.primaries[0])

            # flatten secondaries
            secondaries = []
            for lepton in leptons:
                flattensecondaries(mctree, lepton, secondaries)
            lepton = leptons[0] # single muons for MuonGun, arbitrary for CORSIKA

            # calculate the energy weighted mean position,
            # max position and total energy deposition
            maxparticle = I3Particle()
            maxparticle.energy = 0.
            meanpos = I3Position(0, 0, 0)
            totalenergy = 0
            for sec in secondaries:
                # containment to inner boundary
                if not iscontained(sec.pos.x,
                                   sec.pos.y,
                                   sec.pos.z,
                                   innerboundary,
                                   outeredge_x=outeredge_x,
                                   outeredge_y=outeredge_y):
                    continue
                if sec.energy > maxparticle.energy:
                    maxparticle = I3Particle(sec)
                    maxparticle.energy = sec.energy
                meanpos += sec.pos*sec.energy
                totalenergy += sec.energy
            # fallback solution for clipping events is to increase
            # boundary to outer boundary
            if totalenergy == 0:
                for sec in secondaries:
                    # containment to outer boundary
                    if not iscontained(sec.pos.x,
                                       sec.pos.y,
                                       sec.pos.z,
                                       outerboundary,
                                       outeredge_x,
                                       outeredge_y):
                        continue
                    if sec.energy > maxparticle.energy:
                        maxparticle = I3Particle(sec)
                        maxparticle.energy = sec.energy
                    meanpos += sec.pos*sec.energy
                    totalenergy += sec.energy
            # last fallback solution is to take the lepton information
            if totalenergy == 0:
                totalenergy = lepton.energy
                meanpos = lepton.pos
                meanlength = lepton.length
            else:
                meanpos /= totalenergy
                meanlength = (lepton.pos-meanpos).magnitude

            # for atmospheric muons and for GSR -> tau events there
            # is no hadronic cascade, in that case take maximum energy deposition
            if 'VertexOutgoingHadron' in frame:
                hadron = frame['VertexOutgoingHadron']
            else:
                hadron = maxparticle
                totalenergy -= hadron.energy # need to conserve energy

            # the hadronic cascade
            cascade1 = I3Particle(hadron)
            if 'VertexOutgoingHadron' in frame:
                cascade1 = I3Particle(hadron)
                cascade1.energy = hadron.energy*((ShowerParameters(hadron.type,hadron.energy).emScale))
                ShowerMax = ((ShowerParameters(hadron.type,hadron.energy).a-1)\
                                 /ShowerParameters(hadron.type,hadron.energy).b)
                cascade1.pos = cascade1.shift_along_track(ShowerMax)
                cascade1.time = hadron.time + ShowerMax/I3Constants.c
                cascade1.length = 0.
                cascade1.speed = 0.
                cascade1.shape = I3Particle.Cascade
                cascade1.location_type = I3Particle.InIce


                # the second cascade is the energy weighted average of the track
                cascade2 = I3Particle(lepton)
                cascade2.pos = meanpos
                cascade2.time = lepton.time + meanlength/I3Constants.c
                cascade2.energy = totalenergy
                cascade2.length = 0.
                cascade2.speed = 0.
                cascade2.shape = I3Particle.Cascade
                cascade2.location_type = I3Particle.InIce

                mcdoublebang = I3VectorI3Particle([cascade1, cascade2])
                vertexshift1 = 0.
                vertexshift2 = 0.

    
                
    

        # save it into the frame
        frame.Put('MCInteractionDoubleBang', mcdoublebang)
        frame.Put('VertexShift1', I3Double(vertexshift1))
        frame.Put('VertexShift2', I3Double(vertexshift2))

    
