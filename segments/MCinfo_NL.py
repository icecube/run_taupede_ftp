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




def get_interactiontype_from_frame(frame):
    if frame.Has('I3MCWeightDict'):
        interactiontype = frame['I3MCWeightDict']['InteractionType']
    else:
        tree = frame['I3MCTree']
        for part in tree:
            if not part.is_neutrino:
                continue

            outgoing_lepton = None
            outgoing_neutrino = None
            outgoing_hadron = None
        
            first = True
            for daughter in tree.get_daughters(part):
                if daughter.location_type_string != 'InIce':
                    continue

                if daughter.type_string in ['EMinus', 'EPlus',
                                            'MuMinus', 'MuPlus',
                                            'TauMinus', 'TauPlus']:
                    outgoing_lepton = daughter
                elif daughter.type_string in ['NuE', 'NuEBar',
                                              'NuMu', 'NuMuBar',
                                              'NuTau', 'NuTauBar']:
                    outgoing_neutrino = daughter
                else:
                    outgoing_hadron = daughter

            # interactiontype = 1 : CC
            # interactiontype = 2 : NC
            # interactiontype = 3 : GR
            if outgoing_neutrino is not None and outgoing_hadron is not None:
                interactiontype = 2
            elif outgoing_lepton is not None and outgoing_hadron is not None:
                interactiontype = 1
            else:
                interactiontype = 3
            break

    return interactiontype

@icetray.traysegment
def mcinfo(tray, name,
              innerboundary, outerboundary,
              outeredge_x, outeredge_y,
              dataset='NuTau',
              ethreshold=1e3,
              PhotonsPerBin=15,
              ShowerSpacing=5,name_suffix='',
              **millipede_params):
    
    

    """
    invisible particles w.r.t. to I3MCTree
    the shape is not actually set correctly, and mostly set to Null,
    so e.g. a TauPlus with shape Null will contribute to the
    energy calculation if not taken out manually.
    this is taken care of by not going through the tree
    from the primary but by taking the secondary particles
    "VertexOutgoingLepton" and "VertexOutgoingHadron" as starting points
    in the method addtrueenergydepositions()
    """
    invisibletypes = [I3Particle.NuE, I3Particle.NuEBar,
                     I3Particle.NuMu, I3Particle.NuMuBar,
                     I3Particle.NuTau, I3Particle.NuTauBar,
                    I3Particle.MuPlus,I3Particle.MuMinus,
                    I3Particle.TauPlus,I3Particle.TauMinus,I3Particle.unknown]
    
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
                
             

    """
    compare if two particles are the same without using the particle id
    """
    def issameparticle(particle1, particle2):
        isequal = True
        if (particle1.pos-particle2.pos).magnitude > 1e-3:
            isequal = False
        if np.abs(particle1.dir.zenith-particle2.dir.zenith) > 1e-3:
            isequal = False
        if np.abs(particle1.dir.azimuth-particle2.dir.azimuth) > 1e-3:
            isequal = False
        if np.abs(particle1.energy-particle2.energy) > 1e-3:
            isequal = False
        if np.abs(particle1.time-particle2.time) > 1e-3:
            isequal = False
        if particle1.type != particle2.type:
            isequal = False
        return isequal

    """
    reset the particle id
    """
    def resetparticle(frame):
        mctree = frame['I3MCTree']
        primary = mctree.primaries[0]
        flattenedtree = []
        flattentree(mctree, primary, flattenedtree)
        for particletag in ['VertexOutgoingLepton', 'VertexOutgoingHadron']:
            if not particletag in frame:
               
               continue
            originalparticle = frame[particletag]
            correctparticle = None
            for particle in flattenedtree:
                if issameparticle(originalparticle, particle):
                    correctparticle = I3Particle(particle)
                    break
            
            frame.Delete(particletag)
            
            frame.Put(particletag, correctparticle)
    
    
    tray.Add(resetparticle)

    """
    simple event classification using the interaction type and event topology,
    glashow resonance events can occur in all three neutrino types,
    in NuE naturally, in NuMu and NuTau through neutrino regeneration and decay.
    1: single cascade
    2: double cascade
    3: track
    """
    def interactioneventclassification(frame,
                                       dataset,
                                       outeredge_x,
                                       outeredge_y):
        print('starting interaction classification')
        possible_datasets = ['MuonGun', 'CORSIKA', 'NuE', 'NuMu', 'NuTau']
        if dataset not in possible_datasets:
            raise ValueError(f'dataset not in {possible_datasets}!')
        eventclass = 0
        taulength = 0
        tauenergy = 0
        tauneutrinoenergy = 0
        if dataset in ['MuonGun', 'CORSIKA']:
            eventclass = 3
        elif dataset == 'NuE':
            interactiontype = get_interactiontype_from_frame(frame)
            if interactiontype == 1: # select cc interaction
                eventclass = 1
            elif interactiontype == 2: # select nc interaction
                eventclass = 1
            elif interactiontype == 3: # select glashow resonance
                if 'VertexOutgoingHadron' in frame: # hadronic decay
                    eventclass = 1
                elif 'VertexOutgoingHadronOne' or 'VertexOutgoingHadronTwo' in frame:
                    eventclass = 1
                elif 'VertexOutgoingLepton' in frame:
                    if (frame['VertexOutgoingLepton'].type in 
                            [I3Particle.EMinus, I3Particle.EPlus]):
                        eventclass = 1
                    elif (frame['VertexOutgoingLepton'].type in
                          [I3Particle.MuMinus, I3Particle.MuPlus]):
                        eventclass = 3
                    elif (frame['VertexOutgoingLepton'].type in
                          [I3Particle.TauPlus, I3Particle.TauMinus]):
                        # select events where the detected lepton is a tau
                        mctree = frame['I3MCTree']
                        tau = I3Particle(frame['VertexOutgoingLepton'])
                        tau.shape = I3Particle.InfiniteTrack
                        decaypos = tau.shift_along_track(tau.length)
                        taulength = tau.length
                        tauenergy = tau.energy
                        tauneutrinoenergy = mctree.parent(tau).energy
                        if not frame['TauGeneratesMuon'].value:
                        # select hadronic/electronic decay channel
                        # tau needs to decay within outer boundary
                        # otherwise it is not a double cascade
                        # but considered a single cascade
                        # since there is no hadronic neutrino interaction
                        # this is considered a single cascade
                        # (it is technically a lollipop event)
                            if not iscontained(decaypos.x,
                                               decaypos.y,
                                               decaypos.z,
                                               outerboundary,
                                               outeredge_x,
                                               outeredge_y):
                                eventclass = 3
                            else:
                                eventclass = 1
                        else:
                            eventclass = 3
                            # this is considered a track
                            # (it is technically a sugar daddy)
        elif dataset == 'NuMu':
            interactiontype = get_interactiontype_from_frame(frame)
            #print(interactiontype)
            if interactiontype == 1: # select cc interaction
                eventclass = 3
            elif interactiontype == 2: # select nc interaction
                eventclass = 1
            elif interactiontype == 3: # select glashow resonance
                if 'VertexOutgoingHadron' in frame: # hadronic decay
                    eventclass = 1
                elif 'VertexOutgoingHadronOne' or 'VertexOutgoingHadronTwo' in frame:
                    eventclass = 1
                elif 'VertexOutgoingLepton' in frame:
                    if (frame['VertexOutgoingLepton'].type in
                            [I3Particle.EMinus, I3Particle.EPlus]):
                        eventclass = 1
                    elif (frame['VertexOutgoingLepton'].type in
                          [I3Particle.MuMinus, I3Particle.MuPlus]):
                         print('yes, i was here')
                         print(frame['VertexOutgoingLepton'])
                         eventclass = 3
                    elif (frame['VertexOutgoingLepton'].type in
                          [I3Particle.TauPlus, I3Particle.TauMinus]):
                        # select events where the detected lepton is a tau
                        mctree = frame['I3MCTree']
                        tau = I3Particle(frame['VertexOutgoingLepton'])
                        tau.shape = I3Particle.InfiniteTrack
                        decaypos = tau.shift_along_track(tau.length)
                        taulength = tau.length
                        tauenergy = tau.energy
                        tauneutrinoenergy = mctree.parent(tau).energy
                        if not frame['TauGeneratesMuon'].value:
                        # select hadronic/electronic decay channel
                        # tau needs to decay within outer boundary,
                        # otherwise it is not a double cascade but
                        # considered a single cascade
                        # since there is no hadronic neutrino interaction
                        # this is considered a single cascade
                        # (it is technically a lollipop event)
                            if not iscontained(decaypos.x,
                                               decaypos.y,
                                               decaypos.z,
                                               outerboundary,
                                               outeredge_x,
                                               outeredge_y):
                                eventclass = 3
                            else:
                                eventclass = 1
                        else:
                            eventclass = 3
                            # this is considered a track
                            # (it is technically a sugar daddy)
        elif dataset == 'NuTau':
            interactiontype = get_interactiontype_from_frame(frame)
            if interactiontype == 1: # select cc interaction
                if (frame['VertexOutgoingLepton'].type in
                        [I3Particle.EMinus, I3Particle.EPlus]):
                    eventclass = 1
                elif (frame['VertexOutgoingLepton'].type in
                      [I3Particle.MuMinus, I3Particle.MuPlus]):
                    eventclass = 3
                elif (frame['VertexOutgoingLepton'].type in
                      [I3Particle.TauPlus, I3Particle.TauMinus]):
                    # select events where the detected lepton is in fact a tau
                    mctree = frame['I3MCTree']
                    tau = I3Particle(frame['VertexOutgoingLepton'])
                    tau.shape = I3Particle.InfiniteTrack
                    decaypos = tau.shift_along_track(tau.length)
                    taulength = tau.length
                    tauenergy = tau.energy
                    tauneutrinoenergy = mctree.parent(tau).energy
                    if not frame['TauGeneratesMuon'].value:
                    # select hadronic/electronic decay channel
                    # tau needs to decay within outer boundary,
                    # otherwise it is not a double cascade
                    # but considered a track
                        if not iscontained(decaypos.x,
                                           decaypos.y,
                                           decaypos.z,
                                           outerboundary,
                                           outeredge_x,
                                           outeredge_y):
                            eventclass = 3
                    # NEW: this is to minimize job failures from taus that decay
                    # within abs(x)/abs(y) = 650, but outside simulated volume
                    # (no decay cascade simulated)
                        elif not iscontained(decaypos.x,
                                             decaypos.y,
                                             decaypos.z,
                                             innerboundary,
                                             outeredge_x,
                                             outeredge_y):
                            closestdistance = getclosestdistance(
                                decaypos.x, decaypos.y,
                                outeredge_x, outeredge_y)
                            if closestdistance > 200:
                                print("######## reclassifying event ",
                                      frame["I3EventHeader"].run_id,
                                      frame["I3EventHeader"].event_id,
                                      "with potentially unsimulated tau decay ",
                                      closestdistance, "m outside detector at ",
                                      decaypos.x, "/", decaypos.y,
                                      " as eventclass 3 ########")
                                eventclass = 3
                            else:
                                print("tau decay is ",
                                      closestdistance,
                                      "m outside detector, keep fingers crossed")
                                eventclass = 2
                        else:
                            eventclass = 2
                    else:
                        eventclass = 3
            elif interactiontype == 2: # select nc interaction
                eventclass = 1
            elif interactiontype == 3: # select glashow resonance
                if 'VertexOutgoingHadron' in frame: # hadronic decay
                    eventclass = 1
                elif 'VertexOutgoingHadronOne' or 'VertexOutgoingHadronTwo' in frame:
                    eventclass = 1
                elif 'VertexOutgoingLepton' in frame:
                    if (frame['VertexOutgoingLepton'].type in
                            [I3Particle.EMinus, I3Particle.EPlus]):
                        eventclass = 1
                    elif (frame['VertexOutgoingLepton'].type in
                          [I3Particle.MuMinus, I3Particle.MuPlus]):
                        eventclass = 3
                    elif (frame['VertexOutgoingLepton'].type in
                          [I3Particle.TauPlus, I3Particle.TauMinus]):
                        # select events where the detected lepton is a tau
                        mctree = frame['I3MCTree']
                        tau = I3Particle(frame['VertexOutgoingLepton'])
                        tau.shape = I3Particle.InfiniteTrack
                        decaypos = tau.shift_along_track(tau.length)
                        taulength = tau.length
                        tauenergy = tau.energy
                        tauneutrinoenergy = mctree.parent(tau).energy
                        if not frame['TauGeneratesMuon'].value:
                        # select hadronic/electronic decay channel
                        # tau needs to decay within outer boundary,
                        # otherwise it is not a double cascade but
                        # considered a single cascade
                        # since there is no hadronic neutrino interaction
                        # this is considered a single cascade
                        # (it is technically a lollipop event)
                            if not iscontained(decaypos.x,
                                               decaypos.y,
                                               decaypos.z,
                                               outerboundary,
                                               outeredge_x,
                                               outeredge_y):
                                eventclass = 3
                            else:
                                eventclass = 1
                        else:
                            eventclass = 3
                            # this is considered a track
                            # (it is technically a sugar daddy)

        frame.Put('MCInteractionEventclass'+name_suffix, I3Double(eventclass))
        frame.Put('TauDecayLength'+name_suffix, I3Double(taulength))
        frame.Put('TauEnergy'+name_suffix, I3Double(tauenergy))
        frame.Put('TauNeutrinoEnergy'+name_suffix, I3Double(tauneutrinoenergy))

    tray.AddModule(
        interactioneventclassification,
        dataset=dataset,
        outeredge_x=outeredge_x,
        outeredge_y=outeredge_y)

    """
    add monte carlo double bang information:
    single cascade  :
        first cascade = hadronic + em shower,
        second cascade = zero
    double cascade  :
        first cascade = hadronic shower,
        second cascade = tau decay cascade
    track           :
        first cascade = hadronic shower,
        second cascade = mean of energy weighted track secondaries
    """
    def addtruedoublebang(frame,
                          outeredge_x,
                          outeredge_y):
        # get frame objects
        mctree = frame['I3MCTree']
        eventclass = frame['MCInteractionEventclass'+name_suffix]
        if ('VertexOutgoingLepton' in frame and
                        'VertexOutgoingHadron' in frame):
             print('theres lepton and hadron present')
             
        elif ('VertexOutgoingHadronOne' in frame and
                        'VertexOutgoingHadronTwo' in frame):
               print('both hadrons')
        else:
               print('I dont know whats happening, help!')
        print('eventclass is {0}'.format(eventclass))
        # single cascade
        if eventclass == 1:
            if ('VertexOutgoingLepton' in frame and
                        'VertexOutgoingHadron' in frame):
                    print('Lepton and Hadron')
                    hadron = frame['VertexOutgoingHadron']
                    #print(hadron)
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
            #print('lepton is',lepton)
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
            #print(particle)
            cascade2 = I3Particle(particle)
            cascade2.energy = e
            cascade2.shape = I3Particle.InfiniteTrack
            cascade2.pos = cascade2.shift_along_track(showermax)
            cascade2.time = particle.time + showermax/I3Constants.c
            cascade2.length = 0.
            cascade2.speed = 0.
            cascade2.shape = I3Particle.Cascade
            cascade2.location_type = I3Particle.InIce
            print('True Energy is',e)
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
                print('getting leptons')
            # flatten secondaries
            secondaries = []
            for lepton in leptons:
                flattensecondaries(mctree, lepton, secondaries)
            lepton = leptons[0] # single muons for MuonGun, arbitrary for CORSIKA
            print(lepton)
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
            if 'VertexOutgoingHadron'+name_suffix in frame:
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
        frame.Put('MCInteractionDoubleBang'+name_suffix, mcdoublebang)
        frame.Put('VertexShift1'+name_suffix, I3Double(vertexshift1))
        frame.Put('VertexShift2'+name_suffix, I3Double(vertexshift2))

    tray.AddModule(
        addtruedoublebang, 'truedb',
        outeredge_x=outeredge_x,
        outeredge_y=outeredge_y)

    
    
    """
    more complex event classification using the interaction type and event topology,
    a single cascade always stays a single cascade,
    a double cascade can become a single cascade if the decay length is very short
    or if the decay cascade is not contained,
    a track can become a single cascade if its length is very short
    in NuE naturally, in NuMu and NuTau through neutrino regeneration and decay.
    1: single cascade
    2: double cascade
    3: track
    """
    def reconstructioneventclassification(frame,
                                          ethreshold,
                                          outeredge_x,
                                          outeredge_y):
        # get frame objects
        eventclass = frame['MCInteractionEventclass'+name_suffix]
        cascade1, cascade2 = frame['MCInteractionDoubleBang'+name_suffix]
        length = (cascade1.pos-cascade2.pos).magnitude
        goodevent = ((cascade1.energy >= ethreshold) and
                     (cascade2.energy >= ethreshold))
        # soft energy requirement
        contained1 = iscontained(
            cascade1.pos.x,
            cascade1.pos.y,
            cascade1.pos.z,
            innerboundary,
            outeredge_x,
            outeredge_y)
        # soft inner boundary containment
        contained2 = iscontained(
            cascade2.pos.x,
            cascade2.pos.y,
            cascade2.pos.z,
            innerboundary,
            outeredge_x,
            outeredge_y)
        # soft inner boundary containment

        recoeventclass = 0
        if eventclass == 1:
            recoeventclass = 1
        elif eventclass == 2:
            if not goodevent or not contained1 or not contained2 or (length < 10):
                recoeventclass = 1 # treat as single cascade
            else:
                recoeventclass = 2
        elif eventclass == 3:
            if not goodevent or not contained1 or not contained2 or (length < 10):
                recoeventclass = 1 # treat as single cascade
            else:
                recoeventclass = 3

        # save it into the frame
        frame.Put('MCReconstructionEventclass'+name_suffix, I3Double(recoeventclass))

    tray.AddModule(
        reconstructioneventclassification,
        'reco classification',
        ethreshold=ethreshold,
        outeredge_x=outeredge_x,
        outeredge_y=outeredge_y)    

    """
    add deposited monte carlo and reco energy depositions
    true: sum up over all contained energy depositions
    reco: run amplitude millipede dEdx
    """
    def addtrueenergydepositions(frame,
                                 outeredge_x,
                                 outeredge_y):
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
        #print(parents)
        # get all energy depositions
        secondaries = []
        for parent in parents:
            flattensecondaries(mctree, parent, secondaries)
        if len(secondaries)==0:
            secondaries=parents
        #print(secondaries)
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
                           outeredge_y): #flattensecondaries soft containment
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

        frame.Put('MCTrueTruncateddEdx'+name_suffix, I3VectorI3Particle(truncoutputsecondaries))
        frame.Put('MCTrueTruncatedDepositedEnergy'+name_suffix, I3Double(truncetot))
        frame.Put('MCTruedEdx'+name_suffix, I3VectorI3Particle(outputsecondaries))
        frame.Put('MCTrueDepositedEnergy'+name_suffix, I3Double(etot))
        frame.Put('MCTruedEdxSeedTrack'+name_suffix, seedparticle)

    tray.AddModule(
        addtrueenergydepositions,
        'add true e deposited',
        
        outeredge_x=outeredge_x,
        outeredge_y=outeredge_y)
