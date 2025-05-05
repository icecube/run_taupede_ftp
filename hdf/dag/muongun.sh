#!/bin/bash

set -e  # Exit on error

# # my evt gen env
# /cvmfs/icecube.opensciencegrid.org/py3-v4.3.0/RHEL_7_x86_64/metaprojects/icetray/v1.12.0/env-shell.sh
# source /data/user/tvaneede/software/py_venvs/event_generator_py3-v4.3.0/bin/activate
# /cvmfs/icecube.opensciencegrid.org/py3-v4.3.0/RHEL_7_x86_64/metaprojects/icetray/v1.12.0/env-shell.sh /mnt/ceph1-npx/user/tvaneede/software/py_venvs/event_generator_py3-v4.3.0/bin/python /data/user/tvaneede/GlobalFit/create_hdf_files/MuonGun_NNMFit_hdf.py $@

# neha env
# echo "hoi"
# /cvmfs/icecube.opensciencegrid.org/py3-v4.2.1/RHEL_7_x86_64/metaprojects/icetray/v1.4.1/env-shell.sh
# echo "hier?"
/cvmfs/icecube.opensciencegrid.org/py3-v4.2.1/RHEL_7_x86_64/metaprojects/icetray/v1.4.1/env-shell.sh /cvmfs/icecube.opensciencegrid.org/py3-v4.2.1/RHEL_7_x86_64/bin/python /data/user/tvaneede/GlobalFit/create_hdf_files/MuonGun_NNMFit_hdf.py $@