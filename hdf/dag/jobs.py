
import os

runs = {
    "NuTau" : ["22635", "22634"],
    "NuMu"  : ["22645", "22644"],
    "NuE"   : ["22612", "22613"],
}

base_path =  "/data/user/tvaneede/GlobalFit/run_taupede_ftp/hdf"
this_path = f"{base_path}/dag"

def get_job_list( flavor = "NuTau", reco_version = "v5" ):

    jobs = {}

    directory = f"/data/user/tvaneede/GlobalFit/run_taupede_ftp/reco_files/{reco_version}"

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
dag_version = "v0_bdt"
reco_version = "v5"

outpath = f"/data/user/tvaneede/GlobalFit/run_taupede_ftp/hdf/output/{reco_version}_bdt"

# fixed dag paths
dag_path = f"/scratch/tvaneede//run_taupede_ftp_hdf/{dag_version}/snowstorm_{reco_version}"
log_dir       = f"{dag_path}/logs"

# creating folders and copying scripts
print("creating", dag_path)
os.system(f"mkdir -p {dag_path}")
os.system(f"mkdir -p {log_dir}")
os.system(f"mkdir -p {outpath}")
os.system(f"cp snowstorm.sub {dag_path}")

dagfile = open(f"{dag_path}/submit.dag", 'w')

for flavor in runs:

    print(f"Flavor: {flavor}")

    jobs = get_job_list( flavor = flavor, reco_version = reco_version )
    for run in jobs:
        for folder in jobs[run]:

            inpath = jobs[run][folder]
            outfile = f"{outpath}/{flavor}_{run}_{folder}.hdf5"

            cmd = f"{this_path}/snowstorm.sh --Inpath {inpath} --Outfile {outfile}"
            print(cmd)
            os.system(cmd)

            # jobid = f"hdf_{flavor}_{run}_{folder}"
            # dagfile.write(f"JOB {jobid} snowstorm.sub\n")
            # dagfile.write(f'VARS {jobid} JOBID="{jobid}"\n')
            # dagfile.write(f'VARS {jobid} LOGDIR="{log_dir}"\n')
            # dagfile.write(f'VARS {jobid} Inpath="{inpath}"\n')
            # dagfile.write(f'VARS {jobid} Outfile="{outfile}"\n')

            break
        break
    break
