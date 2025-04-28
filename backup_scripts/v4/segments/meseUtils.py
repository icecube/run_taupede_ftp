from icecube.icetray import traysegment, I3Module, I3Units
from icecube import dataclasses,dataio
from icecube import icetray
import collections
import numpy as np

###Functions needed to check boundary limit of particles.###

import matplotlib.path as mpltPath

def select(geometry):
        r"""Select IceCube DOMs.
        Select all DOMs with an OM type `IceCube` from the given
        detector geometry and sort the selected DOMs per string in
        ascending depth.
        Parameters
        ----------
        geometry : I3OMGeoMap or dict(OMKey, tuple(I3OMGeo, ...))
            Detector geometry
        Returns
        -------
        dict(int, list(tuple(OMKey, I3OMGeo)))
            Mapping of string numbers to sequences of IceCube DOMs
            arranged in ascending depth.
        """
        strings = collections.defaultdict(list)
        # print (type(geometry))
        for omkey, omgeo in geometry.items():
            if np.iterable(omgeo):
                omgeo = omgeo[0]

            if omgeo.omtype == dataclasses.I3OMGeo.IceCube:
                strings[omkey.string].append((omkey, omgeo))

        for doms in strings.values():
            doms.sort(
                key=lambda omgeo: omgeo[1].position.z, reverse=True)


        return strings
    
def boundaries(geometry):
#         Side and top boundaries
#         Find the veto's side and top boundaries.
#         Parameters
#         ----------
#         geometry : I3OMGeoMap or dict(OMKey, tuple(I3OMGeo, ...))
#             IC79 or IC86 detector geometry
#         Returns

#         -------
#         sides : set(int)
#             Sequence of string numbers of the outermost strings
#         top : float
#             Depth in detector coordinates of the first DOM on the
#             deepest non-DeepCore string minus the thickness given
#             by `top_layer`
        
        top_layer=90.*icetray.I3Units.m,
        dust_layer=(-220.*icetray.I3Units.m,-100.*icetray.I3Units.m)

        strings = select(geometry)
        top = min(strings[s][0][1].position.z for s in strings if s <= 78)
   
        neighbors = collections.defaultdict(int)
        dmax = 160.*icetray.I3Units.m

        for string in strings:
            pos = strings[string][0][1].position

            for other in strings:
                if other != string:
                    opos = strings[other][0][1].position

                    # The defined maximum inter-string spacing of 160m between
                    # neighboring DOM assures the "missing" strings of the full
                    # hexagon are treated correctly.
                    if np.hypot(pos.x - opos.x, pos.y - opos.y) < dmax:
                        neighbors[string] += 1

        # The outermost strings have less than six neighbors.
        sides = set(string for string in neighbors if neighbors[string] < 6)
        boundary_x=[]
        boundary_y=[]
        # print('Boundary strings',sides)
        for side_string in sides:
            pos=strings[side_string][0][1].position
            boundary_x.append(pos.x)
            boundary_y.append(pos.y)

        # return sides, top - top_layer[0]
        return boundary_x,boundary_y

def get_surface_det(gcdFile=None):
    
    from icecube import MuonGun
    gcdFile=gcdFile
    bound_2D=[]
    MuonGunGCD='/cvmfs/icecube.opensciencegrid.org/data/GCD/GeoCalibDetectorStatus_AVG_55697-57531_PASS2_SPE_withScaledNoise.i3.gz'
    surface_det = MuonGun.ExtrudedPolygon.from_file(MuonGunGCD, padding=0)##Build Polygon from I3Geometry (convex hull)
    f = dataio.I3File(MuonGunGCD)
    omgeo = f.pop_frame(icetray.I3Frame.Geometry)['I3Geometry'].omgeo
    surface_det_x,surface_det_y=boundaries(omgeo)#getting this from omgeo gives concave hull instead of convex hull
    x=[(surface_det_x[i],surface_det_y[i])for i in range(len(surface_det_x))]###getting only x and y
    bound_2D=mpltPath.Path(x)#Projection of detector on x,y plane
    
    return bound_2D, surface_det

def boundary_check(particle1,gcdFile=None):
    ####checks if particle is inside the detector###
    gcdFile=gcdFile
    bound_2D,surface_det = get_surface_det(gcdFile=gcdFile)
    inlimit = False  
    if ((particle1.pos.z <=max(surface_det.z)) and (particle1.pos.z>=min(surface_det.z))):
        if bound_2D.contains_points([(particle1.pos.x, particle1.pos.y)]):
            inlimit=True

    return inlimit
###begin traysegments###

@traysegment
def RemoveUnwantedKeys(tray, name):
    """
    Remove  unwanted keys accumulated through the various processing levels.
    """
    from icecube import icetray, dataclasses
    
    UnwantedKeys=['BeaconLaunches',
    'CVMultiplicity', 'CVStatistics', 'CalibratedSLC',
    'CascadeContainmentTagging_L2', 'CascadeContainmentTagging_Singles_L2',
    'CascadeDipoleFit_L2', 'CascadeDipoleFit_L2Params', 'CascadeFillRatio_L2',
    'CascadeImprovedLineFit_L2', 'CascadeImprovedLineFit_L2Params', 'CascadeLast_IC_Singles_L2',
    'CascadeLast_IC_Singles_L2Params', 'CascadeLast_L2', 'CascadeLast_L2Params',
    'CascadeLineFitSplit1_L2', 'CascadeLineFitSplit1_L2Params', 'CascadeLineFitSplit2_L2', 'CascadeLineFitSplit2_L2Params',
    'CascadeLineFit_L2', 'CascadeLineFit_L2Params', 'CascadeLlhVertexFitSplit1_L2', 'CascadeLlhVertexFitSplit1_L2Params',
    'CascadeLlhVertexFitSplit2_L2', 'CascadeLlhVertexFitSplit2_L2Params',
    'CascadeLlhVertexFit_IC_Singles_L2', 'CascadeLlhVertexFit_IC_Singles_L2Params',
    'CascadeLlhVertexFit_L2', 'CascadeLlhVertexFit_L2Params',
    'CascadeSplitPulses_L21', 'CascadeSplitPulses_L22',
    'CascadeToISplit1_L2', 'CascadeToISplit1_L2Params', 'CascadeToISplit2_L2', 'CascadeToISplit2_L2Params',
    'CorsikaMoonMJD', 'CorsikaSunMJD',
    'CscdL2_Topo1stPulse_HLC0', 'CscdL2_Topo1stPulse_HLCSplitCount', 'CscdL2_Topo_HLCFirstPulses',
    'EHEATWDPortiaPulseSRT', 'EHEATWDPulseSeriesSRT', 'EHEBestPortiaPulseSRT',
    'EHEDSTShieldParameters_ImpLF', 'EHEDSTShieldParameters_SPE12', 'EHEFADCPortiaPulseSRT',
    'EHEFADCPulseSeriesSRT', 'EHEOpheliaBTWSRT', 'EHEOpheliaParticleBTWSRT', 'EHEOpheliaParticleSRT',
    'EHEOpheliaParticleSRT_ImpLF', 'EHEOpheliaSRT', 'EHEOpheliaSRT_ImpLF', 'EHEPortiaEventSummary',
    'EHEPortiaEventSummarySRT',
    'Estres_CausalQTot', 'Estres_Homogenized_QTot',
    'FilterMask_NullSplit0',
    'FiniteRecoCuts', 'FiniteRecoFit', 'FiniteRecoLlh', 'Homogenized_QTot', 'HuberFit',
    'LargestOMKey', 'LineFit', 'LineFitEHE', 'LineFitEHEParams', 'LineFitParams',
    'MPEFit', 'MPEFitCharacteristics', 'MPEFitCramerRaoParams', 'MPEFitFitParams', 'MPEFitMuEX',
    'OfflinePulsesHLC', 'OfflinePulsesSLC',
    'OnlineL2_BayesianFit', 'OnlineL2_BayesianFitFitParams',
    'OnlineL2_BestFit', 'OnlineL2_BestFitFitParams', 'OnlineL2_BestFit_CramerRaoParams',
    'OnlineL2_BestFit_CramerRao_cr_azimuth', 'OnlineL2_BestFit_CramerRao_cr_zenith',
    'OnlineL2_BestFit_DirectHitsA', 'OnlineL2_BestFit_DirectHitsB', 'OnlineL2_BestFit_DirectHitsC',
    'OnlineL2_BestFit_DirectHitsD', 'OnlineL2_BestFit_DirectHitsE', 'OnlineL2_BestFit_MuEx', 'OnlineL2_BestFit_MuEx_r',
    'OnlineL2_BestFit_Name', 'OnlineL2_CleanedMuonPulses', 'OnlineL2_HitMultiplicityValues', 'OnlineL2_HitMultiplicityValuesIC',
    'OnlineL2_HitStatisticsValues', 'OnlineL2_HitStatisticsValuesIC', 'OnlineL2_MPEFit', 'OnlineL2_MPEFitFitParams',
    'OnlineL2_SPE2itFit', 'OnlineL2_SPE2itFitFitParams', 'OnlineL2_SplineMPE', 'OnlineL2_SplineMPEFitParams',
    'OnlineL2_SplineMPE_Characteristics', 'OnlineL2_SplineMPE_CharacteristicsIC', 'OnlineL2_SplineMPE_CharacteristicsNoRCut',
    'OnlineL2_SplineMPE_CharacteristicsNoRCutIC', 'OnlineL2_SplineMPE_CramerRaoParams',
    'OnlineL2_SplineMPE_CramerRao_cr_azimuth', 'OnlineL2_SplineMPE_CramerRao_cr_zenith',
    'OnlineL2_SplineMPE_DirectHitsA', 'OnlineL2_SplineMPE_DirectHitsB',
    'OnlineL2_SplineMPE_DirectHitsC', 'OnlineL2_SplineMPE_DirectHitsD', 'OnlineL2_SplineMPE_DirectHitsE',
    'OnlineL2_SplineMPE_DirectHitsICA', 'OnlineL2_SplineMPE_DirectHitsICB', 'OnlineL2_SplineMPE_DirectHitsICC',
    'OnlineL2_SplineMPE_DirectHitsICD', 'OnlineL2_SplineMPE_DirectHitsICE', 'OnlineL2_SplineMPE_MuE',
    'OnlineL2_SplineMPE_MuEx', 'OnlineL2_SplineMPE_MuEx_r', 'OnlineL2_SplineMPE_TruncatedEnergy_AllBINS_MuEres',
    'OnlineL2_SplineMPE_TruncatedEnergy_AllBINS_Muon', 'OnlineL2_SplineMPE_TruncatedEnergy_AllBINS_Neutrino',
    'OnlineL2_SplineMPE_TruncatedEnergy_AllDOMS_MuEres', 'OnlineL2_SplineMPE_TruncatedEnergy_AllDOMS_Muon',
    'OnlineL2_SplineMPE_TruncatedEnergy_AllDOMS_Neutrino', 'OnlineL2_SplineMPE_TruncatedEnergy_BINS_MuEres',
    'OnlineL2_SplineMPE_TruncatedEnergy_BINS_Muon', 'OnlineL2_SplineMPE_TruncatedEnergy_BINS_Neutrino',
    'OnlineL2_SplineMPE_TruncatedEnergy_DOMS_MuEres', 'OnlineL2_SplineMPE_TruncatedEnergy_DOMS_Muon',
    'OnlineL2_SplineMPE_TruncatedEnergy_DOMS_Neutrino', 'OnlineL2_SplineMPE_TruncatedEnergy_ORIG_Muon',
    'OnlineL2_SplineMPE_TruncatedEnergy_ORIG_Neutrino', 'OnlineL2_SplineMPE_TruncatedEnergy_ORIG_dEdX',
    'OnlineL2_SplitGeo1_BayesianFit', 'OnlineL2_SplitGeo1_BayesianFitFitParams',
    'OnlineL2_SplitGeo1_Linefit', 'OnlineL2_SplitGeo1_SPE2itFit', 'OnlineL2_SplitGeo1_SPE2itFitFitParams',
    'OnlineL2_SplitGeo2_BayesianFit', 'OnlineL2_SplitGeo2_BayesianFitFitParams', 'OnlineL2_SplitGeo2_Linefit',
    'OnlineL2_SplitGeo2_SPE2itFit', 'OnlineL2_SplitGeo2_SPE2itFitFitParams', 'OnlineL2_SplitTime1_BayesianFit',
    'OnlineL2_SplitTime1_BayesianFitFitParams', 'OnlineL2_SplitTime1_Linefit',
    'OnlineL2_SplitTime1_SPE2itFit', 'OnlineL2_SplitTime1_SPE2itFitFitParams',
    'OnlineL2_SplitTime2_BayesianFit', 'OnlineL2_SplitTime2_BayesianFitFitParams',
    'OnlineL2_SplitTime2_Linefit', 'OnlineL2_SplitTime2_SPE2itFit', 'OnlineL2_SplitTime2_SPE2itFitFitParams',
    'PassedAnyFilter', 'PassedConventional', 'PassedKeepSuperDSTOnly',
    'PoleCascadeFilter_CscdLlh', 'PoleEHESummaryPulseInfo', 'PoleMuonLinefit', 'PoleMuonLlhFit', 'PoleMuonLlhFitFitParams',
    'RTTWOfflinePulses_FR_WIMP', 'SPEFit12EHE', 'SPEFit12EHEFitParams',
    'SPEFit2', 'SPEFit2Characteristics', 'SPEFit2CramerRaoParams', 'SPEFit2FitParams', 'SPEFit2MuEX_FSS',
    'SPEFitSingle', 'SPEFitSingleEHE', 'SPEFitSingleEHEFitParams', 'SPEFitSingleFitParams', 'SRTInIcePulses',
    'SRTInIcePulses_IC_Singles_L2CleanedKeys', 'SRTInIcePulses_WODCCleanedKeys',
    'SplineMPE_Estres', 'SplineMPE_EstresFitParams', 
    'TWOfflinePulsesHLC', 'TWOfflinePulses_FR_WIMP', 'TWOfflinePulses_FR_WIMPTimeRange',
    'UncleanedInIcePulsesTimeRange', 'WIMPrecoTopoSplitSplitCount',
    'splittedDOMMap','splittedDOMMapSRT',
    'IceTopRawData', 'JEBEventInfo', 'OfflineIceTopHLCPulseInfo', 
    'OfflineIceTopHLCTankPulses', 'OfflineIceTopHLCVEMPulses', 'OfflineIceTopSLCVEMPulses', 'TankPulseMergerExcludedStations',
    'CleanTriggerHierarchy_IT']
    
    tray.AddModule('Delete', 'initial_clean', Keys=UnwantedKeys)


@traysegment
def MCInfo(tray, name,gcdFile=None):
    """
    Get info about neutrino type and interaction type and save info in frame
    """
    from icecube import dataclasses, MuonGun
    import numpy as np
    def check_type(frame):
        if "I3MCWeightDict" in frame:
            if frame["I3MCWeightDict"]["PrimaryNeutrinoType"]==12 or frame["I3MCWeightDict"]["PrimaryNeutrinoType"]==-12:
                frame["IsCascade_true"]= icetray.I3Bool(True)
                frame["IsTrack_true"]= icetray.I3Bool(False)

            if (frame["I3MCWeightDict"]["PrimaryNeutrinoType"]==14 or frame["I3MCWeightDict"]["PrimaryNeutrinoType"]==-14) and frame["I3MCWeightDict"]["InteractionType"]==2.0:           
                frame["IsCascade_true"]= icetray.I3Bool(True)
                frame["IsTrack_true"]= icetray.I3Bool(False)
            
            elif (frame["I3MCWeightDict"]["PrimaryNeutrinoType"]==14 or frame["I3MCWeightDict"]["PrimaryNeutrinoType"]==-14) and frame["I3MCWeightDict"]["InteractionType"]==1.0:
                    frame["IsCascade_true"]= icetray.I3Bool(False)
                    frame["IsTrack_true"]= icetray.I3Bool(True)

            if (frame["I3MCWeightDict"]["PrimaryNeutrinoType"]==16 or frame["I3MCWeightDict"]["PrimaryNeutrinoType"]==-16) and frame["I3MCWeightDict"]["InteractionType"]==2.0: 
                frame["IsCascade_true"]= icetray.I3Bool(True)
                frame["IsTrack_true"]= icetray.I3Bool(False)
            
            elif  (frame["I3MCWeightDict"]["PrimaryNeutrinoType"]==16 or frame["I3MCWeightDict"]["PrimaryNeutrinoType"]==-16) and frame["I3MCWeightDict"]["InteractionType"]==1.0:           
                frame["IsCascade_true"]= icetray.I3Bool(False)
                frame["IsTrack_true"]= icetray.I3Bool(True)
            
            elif  (frame["I3MCWeightDict"]["PrimaryNeutrinoType"]==16 or frame["I3MCWeightDict"]["PrimaryNeutrinoType"]==-16) and frame["I3MCWeightDict"]["InteractionType"]==3.0:###these are all glashow events                       
                frame["IsCascade_true"]= icetray.I3Bool(False)
                frame["IsTrack_true"]= icetray.I3Bool(True)
            
            frame["MCPrimaryType"]=dataclasses.I3Double(frame["I3MCWeightDict"]["PrimaryNeutrinoType"])
            frame["NEvents"]=dataclasses.I3Double(frame["I3MCWeightDict"]["NEvents"])
            frame["PrimaryEnergyMC"]=dataclasses.I3Double(frame["I3MCWeightDict"]["PrimaryNeutrinoEnergy"])
            frame["PrimaryZenithMC"]=dataclasses.I3Double(frame["I3MCWeightDict"]["PrimaryNeutrinoZenith"])
            frame["PrimaryAzimuthMC"]=dataclasses.I3Double(frame["I3MCWeightDict"]["PrimaryNeutrinoAzimuth"])
            frame["OneWeight"]=dataclasses.I3Double(frame["I3MCWeightDict"]["OneWeight"])
    tray.Add(check_type)

    def collectStats(frame):
        if 'I3MCTree_preMuonProp' in frame:
                mctree = frame['I3MCTree_preMuonProp']
                bound_2D,surface_det = get_surface_det(gcdFile=gcdFile)
                e_dep_total=0
                neutrino = None
                for p in mctree:
                    depth = mctree.depth(p)
                    if (depth == 0):
                        if neutrino is None and len(mctree.get_daughters(p))>0:
                            neutrino = p
                            Truth_Energies=p.energy
                            intersection=surface_det.intersection(p.pos, p.dir)#points of intersection
                            z_inter=p.pos.z-intersection.first*np.cos(p.dir.zenith)
                            depth=1948.07-z_inter
                            frame["PrimaryDepthMC"]=dataclasses.I3Double(depth)
                            break
                
    tray.AddModule(collectStats)
    def add_primary(frame):
        if "I3MCTree_preMuonProp" in frame:###switch to premuon prop
            def sanitize(particle):
                if particle is None:
                    return dataclasses.I3Particle()
                else:
                    return particle

            mcTree = frame["I3MCTree_preMuonProp"]
            primary = None
            neutrino = None
            for p in mcTree:
                if mcTree.depth(p) != 0: continue

                if p.is_neutrino:
                    if neutrino is None or p.energy > neutrino.energy:
                        neutrino = p

                if primary is None or p.energy > primary.energy:
                    primary = p
            if "MCPrimary" not in frame:
                frame["MCPrimary"] = sanitize(primary)
    tray.Add(add_primary)

@traysegment
def Intersections(tray, name, TimePadding = 60.*I3Units.m/dataclasses.I3Constants.c,TrackName="MMCTrackList",
                  OutputTimeWindow="ContainedTimeRange"):
    """
    We are checking if the events are inside the detector volume or not. This is just for sanity checks. 
    The analysis does not use this explicitly
    """
    from icecube import VHESelfVeto,simclasses
    import math
    def _FindDetectorVolumeIntersections(frame, recoParticle, geometry):
        intersectionPoints = VHESelfVeto.IntersectionsWithInstrumentedVolume(geometry, recoParticle)
        intersectionTimes = []
        IntersectionPoint=dataclasses.I3Position()
        for intersectionPoint in intersectionPoints:
            vecX = intersectionPoint.x - recoParticle.pos.x
            vecY = intersectionPoint.y - recoParticle.pos.y
            vecZ = intersectionPoint.z - recoParticle.pos.z

            prod = vecX*recoParticle.dir.x + vecY*recoParticle.dir.y + vecZ*recoParticle.dir.z
            dist = math.sqrt(vecX**2 + vecY**2 + vecZ**2)
            if prod < 0.: dist *= -1.

            if abs(prod-dist) > 1e-3*icetray.I3Units.m:
                    raise RuntimeError("intersection points are not on track")

            intersectionTimes.append(dist/dataclasses.I3Constants.c + recoParticle.time)

        min_index = 0
        for i in range(len(intersectionTimes)):
            if i==0:
                min_time = intersectionTimes[i]
            else:
                if intersectionTimes[i]<min_time:
                    min_index = i
                    min_time = intersectionTimes[i]

        if len(intersectionTimes)==0:
            IntersectionPoint = dataclasses.I3Position()
        else:
            IntersectionPoint= dataclasses.I3Position(intersectionPoints[min_index])
        sortedTimes = sorted(intersectionTimes)
        return sortedTimes,IntersectionPoint

    def FindDetectorVolumeIntersections(frame, TrackName="", OutputTimeWindow=None, TimePadding=0.):
            if TrackName in frame:
                if OutputTimeWindow is not None:
                        twName = OutputTimeWindow
                else:
                        twName = TrackName + "TimeRange"
                tWindow_List=dataclasses.I3TimeWindowSeries()
                IntersectionPoint_List=dataclasses.I3VectorI3Position()
                for track in frame[TrackName]:
                    theTrack = track
                    geometry=frame['I3Geometry']

                    times,IntersectionPoint = _FindDetectorVolumeIntersections(frame, theTrack.GetI3Particle(), geometry)
                    try:
                        if len(times) == 0:
                                #raise RuntimeError("track does not intersect the detector volume")
                                tWindow_List.append(dataclasses.I3TimeWindow())
                                IntersectionPoint_List.append(IntersectionPoint)
                        elif len(times) == 1:
                                raise RuntimeError("tracks with only one intersection are not supported")
                        else:
                                tWindow = dataclasses.I3TimeWindow(times[0]-TimePadding, times[-1]+TimePadding)
                                tWindow_List.append(dataclasses.I3TimeWindow(tWindow))
                                IntersectionPoint_List.append(IntersectionPoint)
                    except:
                        continue
                                # frame[twName] = tWindow
                frame["IntersectionPoints"]=IntersectionPoint_List
                frame[twName]=tWindow_List   
    tray.AddModule(FindDetectorVolumeIntersections,TimePadding = 60.*I3Units.m/dataclasses.I3Constants.c,
        TrackName="MMCTrackList",
        OutputTimeWindow="ContainedTimeRange")

def getRecoPulses(frame,name):
        pulses = frame[name]
        if pulses.__class__ == dataclasses.I3RecoPulseSeriesMapMask:
                pulses = pulses.apply(frame)
        return pulses

@traysegment
def getCCWD(tray,name,Pulses="TWTSInIcePulses_NoDC",Track="TrackFit",gcdFile=None):
    from icecube import  phys_services
    
    
    def ComputeChargeWeightedDist(frame, Pulses, Track):
        if(not frame.Stop==icetray.I3Frame.Physics):
                    return
        if(not frame.Has(Pulses)):
                    return
        if(not frame.Has(Track)):
                    return
        pulses=getRecoPulses(frame,Pulses)
        track=frame[Track]
        if(track.__class__==dataclasses.I3String):
                    Track=track.value
                    if(not frame.Has(Track)):
                            return
                    track=frame[Track]
        geo=frame.Get('I3Geometry')
        omgeo=geo.omgeo

        Qtot=0
        AvgDistQ=0
        for dom in pulses:
                    DomPosition=omgeo[dom[0]].position
                    Dist=phys_services.I3Calculator.closest_approach_distance(track,DomPosition)
                    Qdom=0
                    for pulse in dom[1]:
                            Qdom+=pulse.charge
                    Qtot+=Qdom
                    AvgDistQ+=Dist*Qdom
        if(Qtot==0):
                    AvgDistQ=NaN
        else:
                    AvgDistQ/=Qtot
        if not frame.Has(Track+'_AvgDistQ'):
                    frame.Put(Track+"_AvgDistQ",dataclasses.I3Double(AvgDistQ))
        if not frame.Has(Pulses+"_Qtot"):
                    frame.Put(Pulses+"_Qtot",dataclasses.I3Double(Qtot))
    
    tray.AddModule(ComputeChargeWeightedDist,name,Pulses=Pulses,Track=Track)

@traysegment
def get_track_length(tray,name,gcdFile=None):
    """
    Track length for cascade/track discrimination
    """
    gcdFile=gcdFile
    def track_length(frame):
        if "Millipede_SplineMPE_TWTS" in frame:  
            losses = frame["Millipede_SplineMPE_TWTS"]
            first_loss_found=False
            for i in range(0,len(losses)):
                if first_loss_found==False and losses[i].energy>1 and boundary_check(losses[i],gcdFile):
                   first_loss = losses[i]
                   first_loss_found=True
                if losses[i].energy>1 and boundary_check(losses[i],gcdFile):
                    last_loss = losses[i]
            if first_loss_found==False:
                frame["TrackLength"] = dataclasses.I3Double(0.)
            else:
                dist = (first_loss.pos-last_loss.pos).magnitude
                frame["TrackLength"] = dataclasses.I3Double(dist)    
        else:
                frame["TrackLength"] = dataclasses.I3Double(0.)
    tray.Add(track_length)

@traysegment
def do_L4_filter(tray,name):
    
    def L4_filter(frame):
        ###########
        ##  HESE
        ###########
        if (frame["IsHESE_ck"].value==True):##On 18 Feb 2022 changed Hese condition to just Hese_ck
            frame["IsHese"]=icetray.I3Bool(True)
            print (frame["I3EventHeader"].run_id, frame["I3EventHeader"].event_id, " HESE")
        else:
            frame["IsHese"]=icetray.I3Bool(False)

        ##  Downgoing Track Veto
        ###########
        qtot = frame['HomogenizedQTot_TWTS'].value
        veto_charge = frame['VetoTrackTWTSInIceL5VetoCharge'].value
        if (qtot<1000. and veto_charge>=0.5):
            frame["IsVetoTrack_New"]=icetray.I3Bool(True)
        elif (qtot>=1000. and veto_charge>=2):
            frame["IsVetoTrack_New"]=icetray.I3Bool(True)
        else:
            frame["IsVetoTrack_New"]=icetray.I3Bool(False)

        ###########
        ##  Starting events
        ###########

        if frame["IsVetoTrack_New"].value==False:
            frame["IsStartingEvent_L4"]=icetray.I3Bool(True)
            print (frame["I3EventHeader"].run_id, frame["I3EventHeader"].event_id, " Starting Event",frame['HomogenizedQTot'].value,frame['HomogenizedQTot_toposplit'].value,frame['HomogenizedQTot_TWTS'].value)
        else:
            frame["IsStartingEvent_L4"]=icetray.I3Bool(False)


        ###########
        ##  Upgoing muons
        ###########
        track_veto_charge=frame['VetoTrackTWTSInIceL5VetoCharge'].value
        if ((( frame['UpgoingTrackTWTSInIceL5VetoCharge'].value > 10 and frame['UpgoingTrackTWTSInIceL5VetoCharge'].value > \
            track_veto_charge and frame['UpgoingTrackTWTSInIceL5VetoChannels'].value > 3 ) \
        or ( frame['UpgoingTrackMilliTWTSInIceVetoCharge'].value > 10 and frame['UpgoingTrackMilliTWTSInIceVetoCharge'].value > \
            track_veto_charge and frame['UpgoingTrackMilliTWTSInIceVetoChannels'].value > 3 ))\
        and (frame["TrackFit"].dir.zenith>1.5)):
            frame["IsUpgoingMuon_L4"]=icetray.I3Bool(True)
            print ("Upgoing Muon event!", frame["I3EventHeader"].run_id, frame["I3EventHeader"].event_id)
        else:
            frame["IsUpgoingMuon_L4"]=icetray.I3Bool(False)

        if frame["IsHese"].value==True:
            return True
        else:
            if (frame['IsVetoTrack_New'].value==False or frame["IsUpgoingMuon_L4"].value==True):
                return True
            else:
                return False
    tray.Add(L4_filter)

@traysegment
def cascade_track_discrimination(tray,name):
    
    def discriminate(frame):
        if frame["IsUpgoingMuon_L4"].value==True:
            frame["IsCascade_recoL4"]=icetray.I3Bool(False)
            frame["IsTrack_recoL4"]=icetray.I3Bool(True)

        if ((frame["L5MonopodFit4_AvgDistQ"].value-30)*0.9<frame["TrackFit_AvgDistQ"].value):
            ###events falling below line cascade distance = m*track distance + c (m=1/0.9, c=30) are cascades
            is_track_reco=False
            is_cascade_reco=True
        else:
            is_track_reco=True
            is_cascade_reco=False

        milli_outgoing_q = frame['StartingTrackHLCMilliTWTSInIceVetoCharge'].value
        mono_outgoing_q = frame["StartingTrackHLCTWTSInIceL5VetoCharge"].value

        if frame["IsHese"].value==True:#HESE event
            if not frame.Has("IsTrack_recoL4"):###if upgoing_muon selection has already determined this as track, avoid this event
                ###selecting hese track events
                if (is_track_reco==True and milli_outgoing_q>50.): ##distance condition and outgoing charge>50
                    frame["IsCascade_recoL4"]=icetray.I3Bool(False)
                    frame["IsTrack_recoL4"]=icetray.I3Bool(True)
                ##track length condition (500 instead of 550) and outgoing charge>30
                elif frame["TrackLength"].value>500 and  mono_outgoing_q>30.:
                    frame["IsCascade_recoL4"]=icetray.I3Bool(False)
                    frame["IsTrack_recoL4"]=icetray.I3Bool(True)
                ##millipede outgoing charge>200, monopod outgoing charge >60, millipede outgoing charge>monopod outgoing charge
                elif  (milli_outgoing_q>200 and mono_outgoing_q>60 and milli_outgoing_q>mono_outgoing_q):
                    frame["IsCascade_recoL4"]=icetray.I3Bool(False)
                    frame["IsTrack_recoL4"]=icetray.I3Bool(True)
                ##selecting cascade
                ####monopod outgoing charge=0
                elif mono_outgoing_q==0:
                    frame["IsCascade_recoL4"]=icetray.I3Bool(True)
                    frame["IsTrack_recoL4"]=icetray.I3Bool(False)
                ##everything else is cascade
                else:
                    frame["IsCascade_recoL4"]=icetray.I3Bool(True)
                    frame["IsTrack_recoL4"]=icetray.I3Bool(False) 

        ####non hese starting events
        if (frame["IsStartingEvent_L4"].value==True and frame["IsHese"].value==False and frame["IsUpgoingMuon_L4"].value==False):
            ##selecting track
            ##satisfy cascade/track distance 
            if is_track_reco==True:
                frame["IsCascade_recoL4"]=icetray.I3Bool(False)
                frame["IsTrack_recoL4"]=icetray.I3Bool(True)
            ##monopod outgoing charge >1.5
            elif mono_outgoing_q>1.5:
                frame["IsCascade_recoL4"]=icetray.I3Bool(False)
                frame["IsTrack_recoL4"]=icetray.I3Bool(True)
            ###selecting cascade
            ###monopod outgoing charge =0
            elif mono_outgoing_q==0:
                frame["IsCascade_recoL4"]=icetray.I3Bool(True)
                frame["IsTrack_recoL4"]=icetray.I3Bool(False)
            ##everything else are cascades
            else:
                frame["IsCascade_recoL4"]=icetray.I3Bool(True)
                frame["IsTrack_recoL4"]=icetray.I3Bool(False)

    tray.Add(discriminate)
    
@traysegment
def coincident_event_cut1(tray,name):
    
    def remove_coinc(frame):
        """
        If there are two coincident events, we will get huge Q_weighted_distance
        """
        if frame["L5MonopodFit4_AvgDistQ"].value>150. and frame["TrackFit_AvgDistQ"].value>110:#if CWD is too large, might not be causally connected
            if  frame['HomogenizedQTot'].value<6000:
                return False
            else:
                return True
        else:
            return True
    tray.Add(remove_coinc)
    
@traysegment
def coincident_event_cut2(tray,name):
    import numpy as np
    def coincident_track_cut(frame):
        """
        If two track fits have large angle difference between them, possibly coincident events. 
        Here, events with angle diff>30 deg are rejected
        """
        def vector(zenith, azimuth):
            return [np.sin(zenith)*np.cos(azimuth), np.sin(zenith)*np.sin(azimuth), np.cos(zenith)]

        def unit_vector(vector):
            """ Returns the unit vector of the vector.  """
            return vector / np.linalg.norm(vector)

        def angle_between(v1, v2):
            v1_u = unit_vector(v1)
            v2_u = unit_vector(v2)
            return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

        if frame.Has("SPEFit4_offline") and frame.Has("TrackFit"):
            vec1 = vector(frame["SPEFit4_offline"].dir.zenith,frame["SPEFit4_offline"].dir.azimuth)
            vec2 = vector(frame["TrackFit"].dir.zenith,frame["TrackFit"].dir.azimuth)
            angle_diff = angle_between(vec1,vec2)

            if frame["IsUpgoingMuon_L4"].value==True and angle_diff>0.5:
                print ("track event removed: "+str(frame["I3EventHeader"].run_id)+","+ str(frame["I3EventHeader"].event_id)+","+ str(np.cos(frame["TrackFit"].dir.zenith)))
                return False

            else:
                return True

    tray.Add(coincident_track_cut)

@traysegment
def coincident_event_cut3(tray,name):
    
    def coincident_rlogLcut(frame):
     
        if 'SPEFit4_offlineFitParams' in frame and 'SPEFit4_rlogL' not in frame :#E_reco Energ
            frame['SPEFit4_rlogL']=dataclasses.I3Double(frame['SPEFit4_offlineFitParams'].rlogl)
        if 'SPEFit4_offlineFitParams' in frame and frame['SPEFit4_offlineFitParams'].rlogl >=8.5:
            return False
        if frame["IsHESE_ck"].value==False:
            if (frame['SPEFit4_offlineFitParams'].rlogl>=8 and frame['IsUpgoingMuon_L4'].value==True):#Tighter cut on misreconstructed upgoing muons
                return False

    tray.Add(coincident_rlogLcut)
    
@traysegment
def corsika_bundle_checks(tray,name):
    """
    If it is single muons, we choose muongun over corsika. So this part is used to reject single muon events.
    """
    from icecube.simprod import segments
    from icecube import  phys_services
    from icecube import MuonGun
    def I3MCTpmp_2_I3MCT(frame):
        del frame["MMCTrackList"]
        del frame["I3MCTree"]
    tray.Add(I3MCTpmp_2_I3MCT,Streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics])
    randomService = phys_services.I3SPRNGRandomService(
        seed = 10000,
        nstreams = 200000000,
        streamnum = 100014318)
    tray.Add(segments.PropagateMuons, 'PropagateMuons',
       RandomService=randomService,
       SaveState=True,
       InputMCTreeName="I3MCTree_preMuonProp",
       OutputMCTreeName="I3MCTree")
    def muonstodet(frame, surface):
        detmu = MuonGun.muons_at_surface(frame, surface)
        frame['Muons_at_Surface'] = dataclasses.I3Double(len(detmu))
        for i in range(len(detmu)):
          frame['EnteringMuon_'+str(i)] = detmu[i]
    def ncut(frame):
        if frame.Has('EnteringMuon_1'):
            return True
        else:
            return False
    surface_det = MuonGun.Cylinder(1400*I3Units.m, 700*I3Units.m, dataclasses.I3Position(0.0, 0.0, 0.0))##used in estes
    tray.context['I3RandomService'] = phys_services.I3GSLRandomService(seed = 1)
    tray.Add(muonstodet, surface=surface_det)
    tray.AddModule(ncut, 'ncut')