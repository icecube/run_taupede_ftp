import os, sys

# for topology in ["DoubleCascades"]:
reco_type = "RecowithBfr"
simulation_type = "Baseline"

reco_dir_nu   = f"/data/ana/Diffuse/GlobalFit_Flavor/taupede/SnowStorm/{reco_type}/{simulation_type}"
reco_dir_muon = f"/data/ana/Diffuse/GlobalFit_Flavor/taupede/MuonGun//{reco_type}/"

runs = {
    "NuTau" : ["22049", "22050", "22085", "22086"],
    "NuMu"  : ["22043", "22044", "22079", "22080"],
    "NuE"   : ["22046", "22047", "22082", "22083"],
    "Muon"  : ["21315", "21316", "21317" ]
}

for flavor in runs:

    print(10*"-", flavor)

    for run_number in runs[flavor]:

        total_files = 0

        if flavor == "Muon":
            reco_run_dir = f"{reco_dir_muon}/{run_number}" 
        else:
            reco_run_dir = f"{reco_dir_nu}/{run_number}" 

        # print(reco_run_dir)

        folders = [f for f in os.listdir(reco_run_dir) if os.path.isdir(os.path.join(reco_run_dir, f))]

        for folder in folders:
            folder_path = f"{reco_run_dir}/{folder}"
            total_files += len([f for f in os.listdir(folder_path) if os.path.isfile(os.path.join(folder_path, f))])

        print("run_number", run_number, "files", total_files)
