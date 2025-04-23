# Create DAG for running EventGenerator Reconstruction

The following scripts:
- create_dag.py: python script that loops over the input files and creates a `dag` folder. Make sure to run this script on npx, because it will create this folder in /scratch. From that folder on the submit node, you can run `condor_submit_dag submit.day`.
- reco.sub: will be copied to your dag folder in /scratch. Contains the specifics for the jobs on the cluster.
- wrapper.sh: a simple .sh wrapper around the script that you want to run

