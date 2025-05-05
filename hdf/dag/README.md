# Versions

Here we can create dag jobs for converting the icetray reconstructed files to hdf5 files, while adding some extra information. 

## Workflow

The main scripts 
- jobs.py: neutrino simulations
- jobs_muon.py: muongun simulations 
search for the available files and create the dag files that call the wrappers
- snowstorm.sh: FullPID_NNMFit_hdf.py
- muongun.sh: MuonGun_NNMFit_hdf.py

This creates hdf5 files for every run and files within a subfolder. The following script
- merge.py
merges these into one dataframe

## Versions

### v0.0

Only 1 file for muons and tau neutrinos in Double Cascade channel. Checking if everything works, added the vertex positions.

### v0.1

All files in a folder for muons, to see if I get any events

### v1.0

All files for DoubleCascade channel, neutrinos only. 

### v1.1 

Mistake fixed with missing SnowStormDict