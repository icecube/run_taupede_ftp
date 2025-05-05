
import os

runs = {
    "Muon"  : ["21315", "21316", "21317" ]
}

base_path =  "/data/user/tvaneede/GlobalFit/create_hdf_files"
this_path = f"{base_path}/dag"

def get_job_list( flavor = "Muon", reco_type = "RecowithBfr", simulation_type = "Baseline" ):

    jobs = {}

    directory = f"/data/ana/Diffuse/GlobalFit_Flavor/taupede/MuonGun/{reco_type}"

    for run in runs[flavor]:
        jobs[run] = {}

        directory_run = f"{directory}/{run}"

        folders = [f for f in os.listdir(directory_run) if os.path.isdir(os.path.join(directory_run, f))]
        for folder in folders:
            jobs[run][folder] = os.path.join(directory_run, folder)

        print("run", run, "folders", len(folders))

    return jobs

###
### Topology DoubleCascades
###
flavor = "Muon"
topology = "DoubleCascades"
reco_type = "RecowithBfr"
simulation_type = "Baseline"
dag_version = "v0.0"

# fixed dag paths
dag_path = f"/scratch/tvaneede/create_hdf_files/{dag_version}/muongun_{reco_type}_{simulation_type}"
log_dir       = f"{dag_path}/logs"

# creating folders and copying scripts
print("creating", dag_path)
os.system(f"mkdir -p {dag_path}")
os.system(f"mkdir -p {log_dir}")
os.system(f"cp muongun.sub {dag_path}")

outfile = open(f"{dag_path}/submit.dag", 'w')

for topology in ["DoubleCascades", "Cascades", "Tracks"]:


    jobs = get_job_list( flavor = flavor, reco_type = reco_type, simulation_type = simulation_type )
    for run in jobs:
        for folder in jobs[run]:
            # cmd = f"python {base_path}/MuonGun_NNMFit_hdf.py --Table {reco_type} --Topology {topology} --Dataset {run} --subfolder {folder}"
            # cmd = f"{this_path}/muongun.sh --Table {reco_type} --Topology {topology} --Dataset {run} --subfolder {folder}"
            # print(cmd)
            # os.system(cmd)

            jobid = f"hdf_{topology}_{flavor}_{run}_{folder}"
            outfile.write(f"JOB {jobid} muongun.sub\n")
            outfile.write(f'VARS {jobid} LOGDIR="{log_dir}"\n')
            outfile.write(f'VARS {jobid} JOBID="{jobid}"\n')
            outfile.write(f'VARS {jobid} Table="{reco_type}"\n')
            outfile.write(f'VARS {jobid} Topology="{topology}"\n')
            outfile.write(f'VARS {jobid} Dataset="{run}"\n')
            outfile.write(f'VARS {jobid} subfolder="{folder}"\n')

            break
        break
    break