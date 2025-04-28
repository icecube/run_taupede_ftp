#!/mnt/ceph1-npx/user/tvaneede/software/py_venvs/event_generator_py3-v4.3.0/bin/python

# icecube imports
from icecube import dataio, icetray, dataclasses
from icecube import phys_services, photonics_service, millipede, VHESelfVeto
from icecube.photonics_service import I3PhotoSplineService
from icecube.dataclasses import I3Double, I3Particle, I3Direction, I3Position, I3VectorI3Particle, I3Constants, I3VectorOMKey
from icecube.dataclasses import I3RecoPulse, I3RecoPulseSeriesMap, I3RecoPulseSeriesMapMask, I3TimeWindow, I3TimeWindowSeriesMap
from icecube.icetray import I3Units, I3Frame, I3ConditionalModule, traysegment
from I3Tray import I3Tray
from icecube.millipede import MonopodFit, MuMillipedeFit, TaupedeFit, HighEnergyExclusions, MillipedeFitParams
from icecube.sim_services.label_events import MCLabeler, MuonLabels, CascadeLabels
from icecube.sim_services.label_events import ClassificationConverter

from segments.VHESelfVeto import SelfVetoWrapper
from segments.MCPreProc import mcpreproc
from segments.Level2RecoWrapper import Level2ReconstructionWrapper
from segments.Level3RecoWrapper import Level3ReconstructionWrapper
from segments.MonopodWrapper import MonopodWrapper
from segments.TaupedeWrapper import TaupedeWrapper
from segments.MillipedeWrapper import MillipedeWrapper
from segments.TrueObservables import calculatetrueobservables
from segments.RecoObservables import calculaterecoobservables
from segments.AddOutgoingParticles import mctreeinfo
from segments.MCinfo_NL import mcinfo
from segments.MCInfoWrapper import MCInfoWrapper
from segments.PassingFraction import penetrating_depth, PassingFraction, add_primary
from segments.FinalEventClassification import checkfinaltopology
# python system imports
import sys, os, datetime
from glob import glob
from os.path import join
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import numpy as np
import pickle
from optparse import OptionParser

parser = OptionParser()
usage = """%prog [options]"""
parser.set_usage(usage)

parser.add_option("-i","--infiles",type=str,help="Input Files",dest="infiles",action='store')
parser.add_option("-o","--outfile",type=str,help="Output File",dest="outfile",action='store')
parser.add_option("-g","--gcdfile",type=str,help="GCD File",dest="GCDfile",action='store',default='/data/user/nlad/Ternary_Classifier/Ternary_Classifier/GeoCalibDetectorStatus_2020.Run134142.Pass2_V0.i3.gz')
parser.add_option("--flavor", choices=('NuE', 'NuMu', 'NuTau','MuonGun','data') ,dest='flavor')
parser.add_option("--icemodel", choices=('Spice_3.2.1','Bfr', 'ftp-v3'),default='ftp-v3', help="Specify which ice model to use",dest='icemodel')
parser.add_option("--year",type=str ,help="Year of operation.",default='IC86_2011',dest='year')
parser.add_option("--innerboundary", type=float, default=550., help="Inner detector boundary to determine contained energy depositions.",dest='innerboundary')
parser.add_option("--outerboundary", type=float, default=650., help="Outer detector boundary to determine contained energy depositions.",dest='outerboundary')

(opts,args) = parser.parse_args()
outfile = opts.outfile
if not opts.GCDfile=='None':
    gcdFile=opts.GCDfile
    infiles=[gcdFile, opts.infiles]
else:
    print('No GCDFile')
    infiles=opts.infiles

# initialize timers
starttimes, stoptimes = {}, {}
timekeys = ['level2','level3','monopod', 'taupede', 'millipede','Reco_observables']
for timekey in timekeys:
    starttimes[timekey] = []
    stoptimes[timekey] = []


################################################################
########################## ICETRAY #############################
################################################################
starttime = datetime.datetime.now()
tray = I3Tray()

tray.AddModule('I3Reader', 'reader', FilenameList=infiles)

def timer(frame,tag,key):
    if tag == 'start':
        starttimes[key].append(datetime.datetime.now())
    elif tag == 'stop':
        stoptimes[key].append(datetime.datetime.now())


def filter_nullsplit(frame):
    if frame["I3EventHeader"].sub_event_stream=='NullSplit':
        return False
    else:
        eventid = frame['I3EventHeader'].event_id
        #print('Interaction tyoe is',frame['I3MCWeightDict']['InteractionType'])
        print("*******Currently processing frame %s*******" %eventid)


tray.Add(filter_nullsplit)
 
pulses = 'SplitInIcePulses'
photons_per_bin = 5
shower_spacing = 5

# Icemodel stuff
ice_model = opts.icemodel

if ice_model == 'Spice_3.2.1':
    base = os.path.expandvars('/cvmfs/icecube.opensciencegrid.org/data/photon-tables/splines/cascade_single_spice_3.2.1_flat_z20_a5.%s.fits')
    base_eff = os.path.expandvars('/cvmfs/icecube.opensciencegrid.org/data/photon-tables/splines/cascade_effectivedistance_spice_3.2.1_z20.%s.fits')
    tiltdir = os.path.expandvars('/cvmfs/icecube.opensciencegrid.org/users/vbasu/meta-projects/combo3/icetray/photonics-service/resources/tilt/')
    pxs = I3PhotoSplineService(base % "abs", base % "prob", effectivedistancetable=base_eff % "eff", timingSigma=0, tiltTableDir=tiltdir)
    #uncomment this following line if effectove distance is NOT to be used, by default, one should always use it for this script!
    #pxs = I3PhotoSplineService(base % "abs", base % "prob", timingSigma=0, tiltTableDir=tiltdir)
elif ice_model == 'Bfr':  
    base = os.path.expandvars('/cvmfs/icecube.opensciencegrid.org/data/photon-tables/splines/cascade_single_spice_bfr-v2_flat_z20_a5.%s.fits')
    base_eff = os.path.expandvars('/cvmfs/icecube.opensciencegrid.org/data/photon-tables/splines/cascade_effectivedistance_spice_bfr-v2_z20.%s.fits')
    tiltdir = os.path.expandvars('/cvmfs/icecube.opensciencegrid.org/users/vbasu/meta-projects/combo3/icetray/photonics-service/resources/tilt/')
    pxs = I3PhotoSplineService(base % "abs", base % "prob", effectivedistancetable=base_eff % "eff", timingSigma=0, tiltTableDir=tiltdir)
    #uncomment this following line if effectove distance is NOT to be used, by default, one should always use it for this script!
    #pxs = I3PhotoSplineService(base % "abs", base % "prob", timingSigma=0, tiltTableDir=tiltdir)
elif ice_model == "ftp-v3":
    base_dir = os.path.expandvars('$I3_DATA/photon-tables/splines')
    base = os.path.join(base_dir, 'cascade_single_spice_ftp-v1_flat_z20_a5.%s.fits')
    base_eff = os.path.join(base_dir, 'cascade_effectivedistance_spice_ftp-v1_z20.%s.fits')
    tiltdir = os.path.expandvars('$I3_SRC/photonics-service/resources/tilt/')
    pxs = I3PhotoSplineService(base % "abs", base % "prob",
                            effectivedistancetable=base_eff % "eff",
                            timingSigma=0, tiltTableDir=tiltdir,
                            effectivedistancetableprob=base_eff % "prob",
                            effectivedistancetabletmod=base_eff % "tmod")



tray.Add('Delete', keys=['BrightDOMs', 'DeepCoreDOMs', 'SaturatedDOMs'])
excludedDOMs = tray.Add(HighEnergyExclusions,
    Pulses = pulses,
    BadDomsList = 'BadDomsList',
    CalibrationErrata = 'CalibrationErrata',
    ExcludeBrightDOMs = 'BrightDOMs',
    ExcludeDeepCore = False,
    ExcludeSaturatedDOMs = 'SaturatedDOMs',
    SaturationWindows = 'SaturationTimes') 

millipede_params = {'Pulses': pulses, 'PartialExclusion' : False , 'CascadePhotonicsService' : pxs, 'ExcludedDOMs': excludedDOMs}



################################################################
####################### HELPER CLASSES #########################
################################################################
if opts.GCDfile=='None':
    infile = dataio.I3File(infiles)
    frame=infile.pop_frame(icetray.I3Frame.Geometry)
    
    geometry = frame['I3Geometry'].omgeo
    
else:
    gcdfile = dataio.I3File(opts.GCDfile)
    frame = gcdfile.pop_frame()
    while 'I3Geometry' not in frame:
        frame = gcdfile.pop_frame()
    geometry = frame['I3Geometry'].omgeo

if opts.year == 'IC79_2010':
    strings = [2, 3, 4, 5, 6, 13, 21, 30, 40, 50, 59, 67, 74, 73, 72, 78, 77, 76, 75, 68, 60, 51, 41, 32, 23, 15, 8]
else:
    strings = [1, 2, 3, 4, 5, 6, 13, 21, 30, 40, 50, 59, 67, 74, 73, 72, 78, 77, 76, 75, 68, 60, 51, 41, 31, 22, 14, 7]

outerbounds = {}
cx, cy = [], []
for string in strings:
    omkey = icetray.OMKey(string, 1)
    if geometry.has_key(omkey):
        x, y = geometry[omkey].position.x, geometry[omkey].position.y
        outerbounds[string] = (x, y)
        cx.append(x)
        cy.append(y)
cx, cy = np.asarray(cx), np.asarray(cy)
order = np.argsort(np.arctan2(cx, cy))
outeredge_x = cx[order]
outeredge_y = cy[order]

# ################################################################
#      #################### SelfVeto ######################
# ################################################################

if opts.flavor != 'data':
    tray.Add('TauGeneratesMuon', If=lambda f : 'I3MCTree' in f)

## both these variables are not in the new files
# tray.Add(lambda frame : 'VHESelfVeto' in frame and not frame['VHESelfVeto'].value) 
# tray.Add(lambda frame : 'CausalQTot' in frame and frame['CausalQTot'].value >= 6000)


################################################################
     #################### MC PREPROC ######################
################################################################
if opts.flavor != 'data':
    tray.AddModule(mctreeinfo,'Add vertices')
    tray.Add(MCInfoWrapper,'Add MC info')

    tray.Add(mcinfo, 'mcpreproc_',
                        outeredge_x=outeredge_x, outeredge_y=outeredge_y,
                        innerboundary=opts.innerboundary, outerboundary=opts.outerboundary,
                        dataset=opts.flavor,
                        ethreshold=1e3,
                        PhotonsPerBin=5,
                        ShowerSpacing=5)

################################################################
       ############## LEVEL 2 RECONSTRUCTION #############
################################################################
tray.Add(timer, tag='start', key='level2')
tray.Add(Level2ReconstructionWrapper, 'level2reco', Pulses = pulses)
tray.Add(timer, tag='stop', key='level2')

################################################################
       ############## LEVEL 3 RECONSTRUCTION #############
################################################################

tray.Add(timer, tag='start', key='level3')
tray.Add(Level3ReconstructionWrapper, 'CombinedCascadeSeed_L3',pulses=pulses) 
tray.Add(timer, tag='stop', key='level3')

################################################################
       ############## MONOPOD RECONSTRUCTION ############
################################################################

tray.Add(timer, tag='start', key='monopod')
tray.Add(MonopodWrapper, 'HESEMonopodFit',
     Seed = 'CombinedCascadeSeed_L3',
     Iterations = 4,
     PhotonsPerBin = photons_per_bin ,
     **millipede_params) 
tray.Add(timer, tag='stop', key='monopod')

################################################################
    ############## TAUPEDE RECONSTRUCTION ################
################################################################
tray.Add(timer, tag='start', key='taupede')
tray.Add(TaupedeWrapper, 'HESETaupedeFit',
     Seed = 'HESEMonopodFit',
     Iterations = 4,
     PhotonsPerBin = photons_per_bin ,
     innerboundary=opts.innerboundary,
     outerboundary=opts.outerboundary,
     outeredge_x=outeredge_x,
     outeredge_y=outeredge_y,
     **millipede_params)
tray.Add(timer, tag='stop', key='taupede')

################################################################
      ########## MILLIPEDE ENERGY RECONSTRUCTIONS ###########
################################################################

# tray.Add(timer, tag='start', key='millipede')
# tray.Add(MillipedeWrapper, 'HESEMillipedeFit',
#     Seeds = ['HESEMonopodFit', 'HESETaupedeFit', 'SPEFit16'],
#     PhotonsPerBin = 5,
#     ShowerSpacing = 5,
#     innerboundary=opts.innerboundary,
#     outerboundary=opts.outerboundary,
#     outeredge_x=outeredge_x,
#     outeredge_y=outeredge_y,
#     **millipede_params)
# tray.Add(timer, tag='stop', key='millipede')

################################################################
    ############## TRUE OBSERVABLES ##############
################################################################
# if opts.flavor != 'data':
#     tray.Add(calculatetrueobservables,
#         'calc_true_observables',
#         innerboundary=opts.innerboundary,
#         outeredge_x=outeredge_x, 
#         outeredge_y=outeredge_y)

#     tray.AddModule(add_primary)
#     tray.AddModule(penetrating_depth, gcd=gcdFile, depth_name_suffix='')


################################################################
   ########### RECONSTRUCTED OBSERVABLES ###########
################################################################
    
# tray.Add(timer, tag='start', key='Reco_observables')
# tray.AddModule(calculaterecoobservables,
# 	   'calc_reco_observables',
# 	   innerboundary=opts.innerboundary,
# 	   outeredge_x=outeredge_x,
# 	   outeredge_y=outeredge_y)
# tray.Add(timer, tag='stop', key='Reco_observables')   

# tray.Add(checkfinaltopology)

deletekeys =['CalibratedSLC', 'FilterMask_NullSplit0','ClusterCleaningExcludedTanks','I3MCTree_preMuonProp_RNGState','SimTrimmer','IceTopPulses',\
           'IceTopRawData','OfflineIceTopHLCTankPulses','OfflineIceTopSLCVEMPulses','InIceDSTPulses','HESE_SPEFitSingle',\
           'MPEFitCramerRaoParams','CascadeContainmentTagging_L2','OnlineL2_SplineMPE','CascadeDipoleFit_L2Params','CascadeLast_IC_Singles_L2Params',\
       'CascadeImprovedLineFit_L2','CascadeImprovedLineFit_L2Params','CascadeLast_IC_Singles_L2','HESE_VHESelfVetoVertexPos','CascadeLast_L2Params',\
       'CascadeLineFitSplit1_L2','OnlineL2_SplineMPE_CramerRao_cr_azimuth','CascadeLineFitSplit1_L2Params','OnlineL2_SplineMPE_DirectHitsICA',\
       'CascadeLineFitSplit2_L2','CascadeLineFitSplit2_L2Params','OnlineL2_SplineMPE_CharacteristicsNoRCutIC','CascadeLineFit_L2',\
'CascadeLlhVertexFitSplit2_L2','TankPulseMergerExcludedTanks','EHEDSTShieldParameters_SPE12','Estres_CausalQTot','CascadeLlhVertexFitSplit2_L2Params',\
'MPEFitMuEX','CascadeLlhVertexFit_IC_Singles_L2','CascadeLlhVertexFit_L2','CascadeLlhVertexFit_L2Params','CascadeLlhVertexFit_IC_Singles_L2Params',\
'EHEOpheliaParticleSRT','CascadeToISplit1_L2Params','CascadeLlhVertexFitSplit1_L2Params','CascadeToISplit2_L2','HESE_MuonImprovedLineFit',\
 'CascadeToISplit2_L2Params','OnlineL2_SplineMPE_CharacteristicsIC','CorsikaMoonMJD','OnlineL2_SplitTime1_SPE2itFitFitParams','PassedKeepSuperDSTOnly',\
 'CscdL2_Topo1stPulse_HLC0','CscdL2_Topo1stPulse_HLCSplitCount','OnlineL2_SplitGeo1_BayesianFitFitParams','EHEATWDPortiaPulseSRT','EHEATWDPulseSeriesSRT',\
'EHEBestPortiaPulseSRT','EHEDSTShieldParameters_ImpLF','EHEFADCPulseSeriesSRT','EHEOpheliaParticleBTWSRT','EHEOpheliaParticleSRT_ImpLF',\
'EHEPortiaEventSummary','HESE_CascadeLlhVertexFitParams','EHEOpheliaBTWSRT','HESE_HomogenizedQTot','HESE_CausalQTot','OnlineL2_SplineMPE_CramerRaoParams',\
'OnlineL2_SplineMPE_TruncatedEnergy_ORIG_Neutrino','OnlineL2_BestFit_Name','CascadeDipoleFit_L2','SaturationWindows','HESE_VHESelfVetoVertexTime',\
'HESE_llhratio','OnlineL2_SplineMPE_TruncatedEnergy_DOMS_Neutrino','Estres_Homogenized_QTot','OnlineL2_BestFit_CramerRaoParams','PoleEHEOphelia_ImpLF',\
 'OnlineL2_SplineMPE_MuEx_r','I3DST','OnlineL2_SplineMPE_TruncatedEnergy_DOMS_Muon','CscdL2_Topo_HLCFirstPulses','OnlineL2_SplitTime1_Linefit',\
 'CascadeSplitPulses_L22','HESE_SPEFit2','splittedDOMMap','OfflineIceTopHLCVEMPulses','IceTop_SLC_InTime','SplineMPE_Estres','OnlineL2_HitStatisticsValues',\
'OnlineL2_HitStatisticsValuesIC','OnlineL2_SplineMPE_DirectHitsA','LineFitEHE','MPEFit','MPEFitCharacteristics','OnlineL2_BayesianFit',\
 'OnlineL2_SplitGeo2_BayesianFit','OnlineL2_SplitTime1_BayesianFitFitParams','SPEFitSingleEHEFitParams','OnlineL2_BestFitFitParams',\
 'OnlineL2_BestFit_CramerRao_cr_azimuth','OnlineL2_HitMultiplicityValues','OnlineL2_BestFit_CramerRao_cr_zenith','OnlineL2_HitMultiplicityValuesIC',\
 'CleanIceTopRawData','OnlineL2_BestFit_DirectHitsA','OnlineL2_BestFit_DirectHitsB','OnlineL2_BestFit_DirectHitsC',\
'OnlineL2_BestFit_DirectHitsD','CascadeLineFit_L2Params','OnlineL2_SplitGeo1_SPE2itFitFitParams','OnlineL2_BestFit_DirectHitsE','OnlineL2_BestFit_MuEx',\
 'OnlineL2_BestFit_MuEx_r','CascadeLlhVertexFitSplit1_L2','OnlineL2_CleanedMuonPulses','EHEPortiaEventSummarySRT',
 'OnlineL2_SplineMPE_TruncatedEnergy_AllDOMS_Muon','CascadeContainmentTagging_Singles_L2','OnlineL2_MPEFit','OnlineL2_MPEFitFitParams',\
 'OnlineL2_SplineMPEFitParams','PoleEHESummaryPulseInfo','EHEOpheliaSRT','OnlineL2_BestFit','OnlineL2_SPE2itFit','PoleMuonLlhFit',\
 'OnlineL2_SplineMPE_Characteristics','OnlineL2_SplineMPE_TruncatedEnergy_BINS_Muon','OnlineL2_SplineMPE_CharacteristicsNoRCut',\
 'OnlineL2_SplineMPE_TruncatedEnergy_BINS_Neutrino','OnlineL2_SplineMPE_CramerRao_cr_zenith','CorsikaSunMJD','HESE_CascadeLlhVertexFit',
 'OnlineL2_SplineMPE_DirectHitsB','OnlineL2_SplineMPE_DirectHitsC','OnlineL2_SplitGeo2_BayesianFitFitParams',
 'OnlineL2_SplineMPE_DirectHitsD','CascadeFillRatio_L2','OnlineL2_SplineMPE_DirectHitsE','OnlineL2_SplineMPE_DirectHitsICB','OnlineL2_SplineMPE_DirectHitsICC',\
 'OnlineL2_SplineMPE_DirectHitsICD','OnlineL2_SplineMPE_DirectHitsICE','OnlineL2_BayesianFitFitParams','OnlineL2_SplineMPE_MuEx',\
 'BeaconLaunches','OnlineL2_SplineMPE_TruncatedEnergy_ORIG_dEdX','SplineMPE_EstresFitParams','splittedDOMMapSRT','OnlineL2_SplineMPE_TruncatedEnergy_AllBINS_MuEres',\
'OnlineL2_SplitGeo1_BayesianFit','OnlineL2_SplitGeo1_Linefit','OnlineL2_SplitGeo1_SPE2itFit','SPEFitSingleFitParams','OnlineL2_SplitGeo2_SPE2itFit',\
'OnlineL2_SplitGeo2_SPE2itFitFitParams','CascadeLast_L2','OnlineL2_SplitTime1_SPE2itFit','OnlineL2_SplitTime2_BayesianFit',\
'OnlineL2_SplitTime2_SPE2itFit','OnlineL2_SplitTime2_BayesianFitFitParams','OnlineL2_SplineMPE_TruncatedEnergy_ORIG_Muon','OnlineL2_SplitTime2_Linefit',\
 'OnlineL2_SplitTime2_SPE2itFitFitParams','HESE_MuonImprovedLineFitParams','PoleCascadeFilter_CscdLlh','PoleEHEOpheliaParticle_ImpLF','PoleMuonLinefit',\
 'PoleMuonLlhFitFitParams','SPEFit12EHE','SPEFit12EHEFitParams','SPEFit2Characteristics','MPEFitFitParams','SPEFit2CramerRaoParams',
'SPEFit2CramerRao_ITParams','OnlineL2_SplitGeo2_Linefit','SPEFit2MuEX_FSS','SPEFitSingle','SPEFitSingleEHE','OnlineL2_SplitTime1_BayesianFit',\
 'OnlineL2_SplineMPE_MuE',]

tray.Add('Delete', keys=deletekeys)                                                                   

tray.AddModule('I3Writer', 'writer',DropOrphanStreams=[icetray.I3Frame.DAQ],
               Streams=[icetray.I3Frame.Geometry, icetray.I3Frame.Calibration,
                        icetray.I3Frame.DetectorStatus, icetray.I3Frame.DAQ, icetray.I3Frame.Physics, icetray.I3Frame.Stream('M')],
                #Streams=[icetray.I3Frame.DAQ, icetray.I3Frame.Physics,icetray.I3Frame.Simulation,
              #icetray.I3Frame.Stream('M')],
                   Filename=opts.outfile)

tray.AddModule('TrashCan', 'yeswecan')
tray.Execute(12)
tray.Finish()
duration = datetime.datetime.now() - starttime
print("\t\tFinished I3Tray..")
print("")
print("This took:",duration)
print("")
print("Timing information for each modules is as follows:")
print("")
for timekey in timekeys:

    if len(starttimes[timekey]) == 0:
        continue
    tstart, tstop = np.asarray(starttimes[timekey]), np.asarray(stoptimes[timekey])
    if len(tstart) != len(tstop):
        durations = tstop - tstart[:len(tstop)]
    else:
        durations = tstop - tstart

    print ("\t{} took {}".format(timekey,durations.sum()))