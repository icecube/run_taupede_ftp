from icecube import icetray, phys_services, dataio, dataclasses
from icecube.simprod import segments
from icecube.dataclasses import I3Double, I3Particle, I3Direction, I3Position, I3VectorI3Particle, I3Constants, I3VectorOMKey
from icecube.simclasses import I3MMCTrack
from icecube import MuonGun, simclasses
from icecube.icetray import traysegment, I3Module

@traysegment
def Propagate(tray,name):
    tray.Add("Rename", "rename",
                Keys =  ['MMCTrackList', 'MMCTrackList_orig'])###the tray wouls already contain mmctracklist
    
    randomService = phys_services.I3SPRNGRandomService(
    seed = 10000,
    nstreams = 200000000,
    streamnum = 100014318)
    tray.context['I3RandomService'] = randomService
    # randomService = phys_services.I3GSLRandomService(314159)

    tray.Add(segments.PropagateMuons, 'PropagateMuons',
         RandomService=randomService,
         # RNGStateName=I3MCTree_preMuonProp_RNGState,
         SaveState=True,
         InputMCTreeName="I3MCTree_preMuonProp",
         OutputMCTreeName="I3MCTree")
    
@traysegment
def Get_MuonWeight(tray,name,grid=False,infile=None):
    def harvest_generators(fname,grid):
        """
        Harvest serialized generator configurations from a set of I3 files.
        """
        from icecube.icetray.i3logging import log_info as log
        generator = None
        if grid==True:
            from icecube.dataio import I3FileStagerFile

            fsf = I3FileStagerFile.GridFTPStager()
            local_infile = fsf.GetReadablePath(fname)
            print(local_infile)
            f = dataio.I3File(str(local_infile), 'r')
        else:
            f = dataio.I3File(fname)
        fr = f.pop_frame(icetray.I3Frame.Stream('S'))
        f.close()
        if fr is not None:
            for k in fr.keys():
                v = fr[k]
                if isinstance(v, MuonGun.GenerationProbability):
                    log('%s: found "%s" (%s)' % (fname, k, type(v).__name__), unit="MuonGun")
                    if generator is None:
                        generator = v
                    else:
                        generator += v
        return generator

    model = MuonGun.load_model('GaisserH4a_atmod12_SIBYLL')
    generator = harvest_generators(infile,grid)
    tray.AddModule('I3MuonGun::WeightCalculatorModule', 'MuonWeight', Model=model,
        Generator=generator)
    
@traysegment
def Get_CorsikaWeight(tray,name,flux=None):
    from icecube.weighting.fluxes import GaisserH4a, Hoerandel5,Honda2004
    from icecube.weighting import weighting, get_weighted_primary
    import numpy as np
    corsika_dids = [20904]##this is the dataset mese uses (from m.silva)
    simprod_nfiles={"21002":9979,"21217":21851,"21218":11991,"21219":15376,
                    "21220":9999,"21221":10000,
                    # "20904":740335, 
                    "20904":50000,###not using the entire set?
                    "20891": 497612, "20881": 99904,
                    "20852":99094, "20849":9948, "20848":99782, "20789":99998,
                    "20788":99998, "20787":99743,
                    "scat5":9992,"scat-5":9994,"abs5":9993,"abs-5":9988,
                    "domeff90":9831,"domeff95":9829,"domeff105":9813,"domeff110":9813,
                    "p0=-2.0_p1=-0.2":9822,"p0=-2.0_p1=0.0":9825,"p0=-2.0_p1=0.2":9828,
                    "p0=-1.0_p1=-0.2":9840,"p0=-1.0_p1=0.0":9839,"p0=-1.0_p1=0.2":9830,
                    "p0=0.0_p1=-0.2":9830,"p0=0.0_p1=0.0":9834,"p0=0.0_p1=0.2":9847,
                    "p0=1.0_p1=-0.2":9831,"p0=1.0_p1=0.0":9832,"p0=1.0_p1=0.2":9817
                    }
    def weighter_corsika(frame,flux=None):
        if not frame.Has("MCPrimary"):
            get_weighted_primary(frame, MCPrimary="MCPrimary")
        MCPrimary = frame["MCPrimary"]
        energy = MCPrimary.energy
        ptype = MCPrimary.type
        placeholder = True
        for DID in corsika_dids:
            if(placeholder):
                generator = weighting.from_simprod(DID) * simprod_nfiles["%5.0f"%DID];
                placeholder = False
            else:
                generator += weighting.from_simprod(DID) * simprod_nfiles["%5.0f"%DID];
        if flux is 'GaisserH4a':
            flux=GaisserH4a()
            weights = flux(energy, ptype) / generator(energy, ptype)
            if np.isnan(weights): weights=0.0;
            frame['Weight_GaisserH4a']=dataclasses.I3Double(weights)
        elif flux is 'Hoerandel5':
            flux=Hoerandel5()
            weights = flux(energy, ptype) / generator(energy, ptype)
            if np.isnan(weights): weights=0.0;
            frame['Weight_Hoerandel5']=dataclasses.I3Double(weights)
        elif flux is 'Honda2004':
            flux=Honda2004()
            weights = flux(energy, ptype) / generator(energy, ptype)
            if np.isnan(weights): weights=0.0;
            frame['Weight_Honda2004']=dataclasses.I3Double(weights)
        #do a little cleanup now
        flux = 0; generator =0; weights = 0; MCPrimary = 0;
    tray.AddModule(weighter_corsika, "WeighterCorsika_"+name,flux=flux)

@traysegment
def Get_CorsikaSimWeight(tray,name):
    import simweights
    import numpy as np
    from icecube import dataio

    weight_keys = [
        "CylinderLength",
        "CylinderRadius",
        "EnergyPrimaryMax",
        "EnergyPrimaryMin",
        "NEvents",
        "OverSampling",
        "ParticleType",
        "PrimaryEnergy",
        "PrimarySpectralIndex",
        "PrimaryType",
        "ThetaMax",
        "ThetaMin",
        "Weight",
    ]

    particle_keys = ["type", "energy", "zenith"]
    def simweighter_corsika(frame):
        MCtype_corsika = np.array([])
        MCenergy_corsika = np.array([])
        CorsikaWeightMap = {k: [] for k in weight_keys}
        PolyplopiaPrimary = {k: [] for k in ["type", "energy", "zenith"]}
        if frame.Stop == frame.Physics and "FilterMask" in frame:
            MCtype_corsika = np.append(MCtype_corsika, frame["PolyplopiaPrimary"].type)
            MCenergy_corsika = np.append(MCenergy_corsika, frame["PolyplopiaPrimary"].energy)

            for k in weight_keys:
                CorsikaWeightMap[k].append(frame["CorsikaWeightMap"][k])

            PolyplopiaPrimary["zenith"].append(frame["PolyplopiaPrimary"].dir.zenith)
            PolyplopiaPrimary["type"].append(frame["PolyplopiaPrimary"].type)
            PolyplopiaPrimary["energy"].append(frame["PolyplopiaPrimary"].energy)
        fobj = dict(CorsikaWeightMap=CorsikaWeightMap, PolyplopiaPrimary=PolyplopiaPrimary)
        wobj = simweights.CorsikaWeighter(fobj, nfiles=1)
        weights_GaisserH4a = wobj.get_weights(simweights.GaisserH4a())[0]
        frame['SimWeight_GaisserH4a']=dataclasses.I3Double(weights_GaisserH4a)
        weights_Hoerandel5 = wobj.get_weights(simweights.Hoerandel5())[0]
        frame['SimWeight_Hoerandel5']=dataclasses.I3Double(weights_Hoerandel5)
        weights_Honda2004 = wobj.get_weights(simweights.Honda2004())[0]
        frame['SimWeight_Honda2004']=dataclasses.I3Double(weights_Honda2004)
    tray.AddModule(simweighter_corsika)

