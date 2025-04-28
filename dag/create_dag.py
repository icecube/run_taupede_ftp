import sys, os, glob

# set the inputs
reco_version = "v3"
icemodel = "ftp-v3"

# # SnowStorm_NuTau_highE
# simulation_dataset   = "22635" 
# simulation_subfolder = "0000000-0000999"
# simulation_flavor    = "NuTau"

# # SnowStorm_NuTau_midE
# simulation_dataset   = "22634" 
# simulation_subfolder = "0000000-0000999"
# simulation_flavor    = "NuTau"

# SnowStorm_NuE_highE
simulation_dataset   = "22612" 
simulation_subfolder = "0000000-0000999"
simulation_flavor    = "NuE"

# SnowStorm_NuE_midE
simulation_dataset   = "22613" 
simulation_subfolder = "0000000-0000999"
simulation_flavor    = "NuE"

# fixed paths
work_path = "/data/user/tvaneede/GlobalFit/run_taupede_ftp"
reco_input_path = f"/data/sim/IceCube/2023/filtered/level2/neutrino-generator/{simulation_dataset}/{simulation_subfolder}"

reco_out_path = f"{work_path}/reco_files/{reco_version}/{simulation_dataset}/{simulation_subfolder}"

# fixed dag paths
dag_base_path = "/scratch/tvaneede/reco/run_taupede_ftp"
dag_name = f"reco_dag_{reco_version}_{simulation_dataset}_{simulation_subfolder}"

dag_path      = f"{dag_base_path}/{dag_name}"
log_dir       = f"{dag_path}/logs"
backup_path   = f"{work_path}/backup_scripts/{reco_version}"

# creating folders and copying scripts
print("creating", dag_path)
os.system(f"mkdir -p {dag_path}")
os.system(f"mkdir -p {log_dir}")
os.system(f"mkdir -p {reco_out_path}")
os.system(f"mkdir -p {backup_path}")
os.system(f"cp reco.sub {dag_path}")

# backup scripts
os.system(f"cp {work_path}/dag/reco.sub {backup_path}")
os.system(f"cp {work_path}/dag/wrapper.sh {backup_path}")
os.system(f"cp {work_path}/HESE_Taupede.py {backup_path}")


# create the dag job
outfile = open(f"{dag_path}/submit.dag", 'w')

infiles_list = glob.glob(f"{reco_input_path}/Level2_{simulation_flavor}_*.i3.zst")
print(f"found {len(infiles_list)} files")

import re
infiles_list = sorted(infiles_list, key=lambda x: int(re.search(r'\.(\d{6})\.i3\.zst$', x).group(1)))

i = 0
for INFILES in infiles_list:

    filename = os.path.basename(INFILES)
    JOBID = filename.split("_")[2] # gives the run number
    OUTFILE = f"{reco_out_path}/Reco_{simulation_flavor}_{JOBID}_out.i3.bz2"

    # print(INFILES)
    # print(filename, JOBID)
    # print(OUTFILE)

    outfile.write(f"JOB {JOBID} reco.sub\n")
    outfile.write(f'VARS {JOBID} LOGDIR="{log_dir}"\n')
    outfile.write(f'VARS {JOBID} JOBID="{JOBID}"\n')
    outfile.write(f'VARS {JOBID} INFILES="{INFILES}"\n')
    outfile.write(f'VARS {JOBID} OUTFILE="{OUTFILE}"\n')
    outfile.write(f'VARS {JOBID} FLAVOR="{simulation_flavor}"\n')
    outfile.write(f'VARS {JOBID} ICEMODEL="{icemodel}"\n')

    i+=1
    if i == 50: break
