# icecube imports
from icecube import dataio, icetray, dataclasses
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

             

def addtrueenergydepositions(frame,
                             outeredge_x,
                             outeredge_y,
                            innerboundary,outerboundary):
        def reconstructioneventclassification(frame):
    
            # get frame objects
            eventclass = frame['MCInteractionEventclass']
            cascade1, cascade2 = frame['MCInteractionDoubleBang']
            length = (cascade1.pos-cascade2.pos).magnitude
            goodevent = (cascade1.energy >= ethreshold) and (cascade2.energy >= ethreshold) # soft energy requirement
            contained1 = iscontained(cascade1.pos.x, cascade1.pos.y, cascade1.pos.z, innerboundary,outeredge_x, outeredge_y) # soft inner boundary containment
            contained2 = iscontained(cascade2.pos.x, cascade2.pos.y, cascade2.pos.z, innerboundary,outeredge_x, outeredge_y) # soft inner boundary containment

            recoeventclass = 0
            if eventclass == 1:
                recoeventclass = 1
            elif eventclass == 2:
                if not goodevent or not contained1 or not contained2 or (length < 20):
                    recoeventclass = 1 # treat as single cascade
                else:
                    recoeventclass = 2
            elif eventclass == 3:
                if not goodevent or not contained1 or not contained2 or (length < 20):
                    recoeventclass = 1 # treat as single cascade
                else:
                    recoeventclass = 3

            # save it into the frame
            frame.Put('MCReconstructionEventclass', I3Double(recoeventclass))
            
        #tray.Add(reconstructioneventclassification)
        
        
        # get all interaction particles
        mctree = frame['I3MCTree']
        parents = []
        if 'VertexOutgoingHadron' in frame: # NuGen
            parents.append(frame['VertexOutgoingHadron'])
            print('Parent is Hadron')
        if 'VertexOutgoingLepton' in frame: # NuGen
            parents.append(frame['VertexOutgoingLepton'])
            print('Parent is Lepton')
        if 'VertexOutgoingHadronOne' in frame and 'VertexOutgoingHadronTwo' in frame:
            parents.append(frame['VertexOutgoingHadronOne'])
            parents.append(frame['VertexOutgoingHadronTwo'])
            print('Two Hadrons')
        if ('VertexOutgoingHadron' not in frame and
                'VertexOutgoingLepton' not in frame):
            # MuonGun or CORSIKA
            bundle = mctree.get_daughters(mctree.primaries[0])
            # single muons in MuonGun, bundles in CORSIKA
            for muon in bundle:
                parents.append(muon)
        print(parents)
        # get all energy depositions
        secondaries = []
        for parent in parents:
            flattensecondaries(mctree, parent, secondaries)
        print(secondaries)
        # sort them in time
        timesecmap = []
        for index in range(len(secondaries)):
            sec = secondaries[index]
            timesecmap.append([sec.time, index])
        sortedsecondaries = []
        for time, index in sorted(timesecmap):
            sortedsecondaries.append(secondaries[index])
        
        # truncate dE/dx
        seedparticle = I3Particle(sortedsecondaries[len(sortedsecondaries)//2])
        outputsecondaries, truncoutputsecondaries = [], []
        etot, truncetot = 0, 0
        for sec in sortedsecondaries:
            if iscontained(sec.pos.x,
                           sec.pos.y,
                           sec.pos.z,
                           innerboundary,                          
                           outeredge_x,
                           outeredge_y): # soft containment
                truncoutputsecondaries.append(sec)
                truncetot += sec.energy
            outputsecondaries.append(sec)
            etot += sec.energy
            # use the particle with the closest approach to the detector
            # center as seed for millipede
            if sec.pos.magnitude < seedparticle.pos.magnitude:
                seedparticle = I3Particle(sec)

        # in the rare case that the seed particle is outside the outer
        # boundary, shift it so it is inside and millipede has
        # something ro reconstruct
        seedparticle.shape = I3Particle.InfiniteTrack
        backuppos = I3Position(seedparticle.pos)
        shiftdistance = 0.
        if seedparticle.pos.x > outerboundary:
            shiftdistance = (outerboundary+innerboundary)/(2*np.sin(seedparticle.dir.theta)*np.cos(seedparticle.dir.phi))
            backuppos.x = (outerboundary+innerboundary)/2.
        elif seedparticle.pos.x < -outerboundary:
            shiftdistance = -(outerboundary+innerboundary)/(2*np.sin(seedparticle.dir.theta)*np.cos(seedparticle.dir.phi))
            backuppos.x = -(outerboundary+innerboundary)/2.
        if seedparticle.pos.y > outerboundary:
            shiftdistance = (outerboundary+innerboundary)/(2*np.sin(seedparticle.dir.theta)*np.sin(seedparticle.dir.phi))
            backuppos.y = (outerboundary+innerboundary)/2.
        elif seedparticle.pos.y < -outerboundary:
            shiftdistance = -(outerboundary+innerboundary)/(2*np.sin(seedparticle.dir.theta)*np.sin(seedparticle.dir.phi))
            backuppos.y = -(outerboundary+innerboundary)/2.
        if seedparticle.pos.z > outerboundary:
            shiftdistance = (outerboundary+innerboundary)/(2*np.cos(seedparticle.dir.theta))
            backuppos.z = (outerboundary+innerboundary)/2.
        elif seedparticle.pos.z < -outerboundary:
            shiftdistance = -(outerboundary+innerboundary)/(2*np.cos(seedparticle.dir.theta))
            backuppos.z = -(outerboundary+innerboundary)/2.
        newpos = seedparticle.shift_along_track(shiftdistance)
        if np.abs(newpos.x) > outerboundary or np.abs(newpos.y) > outerboundary or np.abs(newpos.z) > outerboundary:
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

        frame.Put('MCTrueTruncateddEdx', I3VectorI3Particle(truncoutputsecondaries))
        frame.Put('MCTrueTruncatedDepositedEnergy', I3Double(truncetot))
        frame.Put('MCTruedEdx', I3VectorI3Particle(outputsecondaries))
        frame.Put('MCTrueDepositedEnergy', I3Double(etot))
        frame.Put('MCTruedEdxSeedTrack', seedparticle)
