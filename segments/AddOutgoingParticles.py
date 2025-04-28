from icecube import dataclasses, icetray, simclasses, recclasses
from icecube.dataclasses import I3Particle, I3Constants, I3Double
from icecube.dataclasses import I3VectorI3Particle, I3Position, I3Direction
from icecube.icetray import I3Units, I3Frame
from icecube.sim_services import ShowerParameters
import numpy as np

def mctreeinfo(frame):

    tree = frame['I3MCTree_preMuonProp']
    positions = []
    times = []
    outgoing_leptons = []
    outgoing_hadrons = []
    neutrinos = [I3Particle.NuE, I3Particle.NuMu, I3Particle.NuTau,I3Particle.NuEBar, I3Particle.NuMuBar, I3Particle.NuTauBar]
    hadrons = [I3Particle.Hadrons]
    for t in tree:
        
        if frame['I3MCWeightDict']['InteractionType']==-1:
            print('No In-ice neutrino interaction happened')
            continue
            
        if t.is_neutrino:
            continue
        if t.location_type_string != 'InIce':
            continue
            
        parent = tree.parent(t)
        
        if t.type_string in ['EMinus', 'EPlus',
                            'MuMinus', 'MuPlus',
                            'TauMinus', 'TauPlus'] and parent.type in neutrinos:
            outgoing_leptons.append(t)
            positions.append(t.pos)
            times.append(t.time)
            
        elif (t.type in hadrons) and (parent.type in neutrinos):
                
                outgoing_hadrons.append(t)
                positions.append(t.pos)
                times.append(t.time)
       
    frame['NumVerticesInFiducialVolume'] = icetray.I3Int(len(positions))
    frame['VertexPosition'] = positions[0]
    
        
    if (len(outgoing_leptons) != 0) and (len(outgoing_hadrons) !=0):
        frame['VertexOutgoingLepton'] = outgoing_leptons[0]
        frame['VertexOutgoingHadron']= outgoing_hadrons[0]
   
    elif (len(outgoing_leptons) != 0) and (len(outgoing_hadrons) == 0):
        frame['VertexOutgoingLepton'] = outgoing_leptons[0]
        
    elif (len(outgoing_leptons) == 0) and (len(outgoing_hadrons) >1) :
        
        frame['VertexOutgoingHadronOne']= outgoing_hadrons[0]
        frame['VertexOutgoingHadronTwo']= outgoing_hadrons[1]
        
    else:
        frame['VertexOutgoingHadron']= outgoing_hadrons[0]
        
    
    frame['VertexTime'] = dataclasses.I3Double(times[0]) 
 
    return True
