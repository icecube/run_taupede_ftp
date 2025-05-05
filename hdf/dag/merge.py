import os, sys

# for topology in ["DoubleCascades"]:
reco_type = "RecowithBfr"
simulation_type = "Baseline"


reco_dir_nu   = f"/data/ana/Diffuse/GlobalFit_Flavor/taupede/SnowStorm/{reco_type}/{simulation_type}"
reco_dir_muon = f"/data/ana/Diffuse/GlobalFit_Flavor/taupede/MuonGun//{reco_type}/"

hdf_dir_nu    = f"/data/user/tvaneede/datasets/taupede/SnowStorm/NoDeepCore/hdf_files/{reco_type}/{simulation_type}"
hdf_dir_muon  = f"/data/user/tvaneede/datasets/taupede/MuonGun/NoDeepCore/hdf_files/{reco_type}/"
hdf_dir_merg  = f"/data/user/tvaneede/datasets/taupede/merged"

runs = {
    "NuTau" : ["22049", "22050", "22085", "22086"],
    "NuMu"  : ["22043", "22044", "22079", "22080"],
    "NuE"   : ["22046", "22047", "22082", "22083"],
    # "Muon"  : ["21315", "21316", "21317" ]
}

def merge_folders( flavor, run_number, topology ):
    
    if flavor == "Muon":
        reco_run_dir = f"{reco_dir_muon}/{run_number}" 
        hdf_run_dir  = f"{hdf_dir_muon}/{run_number}" 
    else:
        reco_run_dir = f"{reco_dir_nu}/{run_number}" 
        hdf_run_dir  = f"{hdf_dir_nu}/{run_number}" 

    merged_outfile = f"{hdf_run_dir}/{run_number}_{topology}.hdf5" 

    hdf_files_to_merge = []
    cmd = f"hdfwriter-merge -o {merged_outfile}"

    # cross check which folders are in a run
    folders = [f for f in os.listdir(reco_run_dir) if os.path.isdir(os.path.join(reco_run_dir, f))]

    print(20*"-")
    print(flavor, run_number, folders)
    print(merged_outfile)

    # check if all folders have a hdf5 file
    for folder in folders:
        hdf_file_test = f"{hdf_run_dir}/{run_number}_{folder}_{topology}.hdf5"
        if os.path.exists(hdf_file_test): 
            # print( hdf_file_test, "exists" )
            hdf_files_to_merge.append( hdf_file_test )
            cmd+= f" {hdf_file_test}"
        else:
            sys.exit( hdf_file_test + " not exist" )

    # print(cmd)
    os.system(cmd)

# def merge_runs( topology ):

#     merged_outfile = f"{hdf_dir_merg}/{topology}.hdf5" 

#     cmd = f"hdfwriter-merge -o {merged_outfile}"

#     for flavor in runs:
#         for run_number in runs[flavor]:
#             hdf_file_test = f"{hdf_dir_nu}/{run_number}/{run_number}_{topology}.hdf5"
#             if os.path.exists(hdf_file_test): 
#                 print( hdf_file_test, "exists" )
#                 cmd+= f" {hdf_file_test}"
#             else:
#                 sys.exit( hdf_file_test + " not exist" )

    # print(cmd)


# merge_folders( "NuTau", "22050", topology )

topology = "DoubleCascades"

for flavor in runs:
    for run_number in runs[flavor]:
        merge_folders( flavor, run_number, topology )

# merge_runs("DoubleCascades")