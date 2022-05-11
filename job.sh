#!/bin/bash -l
# Batch script to run a serial job on Legion with the upgraded
# software stack under SGE.
# 1. Force bash as the executing shell.
#$ -S /bin/bash
# 2. Request ten minutes of wallclock time (format hours:minutes:seconds).
#$ -l h_rt=24:00:0
# 3. Request 1 gigabyte of RAM (must be an integer)
#$ -l mem=4G
# 4. Request 15 gigabyte of TMPDIR space (default is 10 GB)
#$ -l tmpfs=15G
# 5. Set the name of the job.
#$ -N 11052022_models_feature_extraction
# .Set up array so we can loop over structures
#$ -t 1-5:50
# 6. Set the working directory to somewhere in your scratch space.  This is
# a necessary step with the upgraded software stack as compute nodes cannot
# write to $HOME.
# Replace "<your_UCL_id>" with your UCL user ID :)
#$ -wd /home/regmili/Scratch/WetLabModels
# 7. Your work *must* be done in $TMPDIR
#cd $TMPDIR #
# 8. Run the application.
#load modules
conda deactivate
cd /home/regmili/Scratch/WetLabModels/
mkdir -p all_pdbs/
mkdir -p anarci/
mkdir -p renumbered_pdb/
mkdir -p features_logs/
mkdir -p model_outputs/
mkdir -p Results_by_pdb/distances/
mkdir -p Results_by_pdb/energies/
mkdir -p Results_by_pdb/sequences/
mkdir -p many_models_logs/
mkdir -p outputs_for_model/
cd /home/regmili/Scratch/WetLabModels/fasta_files/
ls | grep ".fasta">/home/regmili/Scratch/WetLabModels/list_to_model.txt
module unload python3/recommended
module load hmmer
input1="/home/regmili/Scratch/WetLabModels/list_to_model.txt"
today=$(date +"%m-%d-%Y")
for (( i=$SGE_TASK_ID; i<$SGE_TASK_ID+100; i++ ))
 do
   conda activate TCRpMHC_env
   cd /home/regmili/Scratch/WetLabModels/fasta_files/
   number=$i
   line="`sed -n ${number}p $input1`"
   id=${line%.fasta}
   echo $id
   structure_name="${id}_TCR-pMHC.pdb"
   /home/regmili/bin/modeller9.23/bin/modpy.sh python /home/regmili/Scratch/PDB/tcrpmhc_models/tcrpmhc_models/__main__.py -t "$line" -n "$id"  -p /home/regmili/Scratch/WetLabModels/model_outputs/
   conda deactivate
   conda activate pyrosetta
   cd /home/regmili/Scratch/TCR-structure-binding-prediction
   git checkout WetLabModels
   cp "/home/regmili/Scratch/WetLabModels/model_outputs/${structure_name}" "/home/regmili/Scratch/WetLabModels/all_pdbs/${structure_name}"
   python ./processing_and_prediction/01_dist_atchley.py ${structure_name} WetLabModels/ > /home/regmili/Scratch/WetLabModels/features_logs/${id}_log.txt
   python ./processing_and_prediction/02_energies.py ${structure_name} WetLabModels/ >> /home/regmili/Scratch/WetLabModels/features_logs/${id}_log.txt
   python ./processing_and_prediction/03_extract_template_info.py ${structure_name} WetLabModels/ >> /home/regmili/Scratch/WetLabModels/features_logs/${id}_log.txt
 done
# # Make sure you have given enough time for the copy to complete!
# #$ -m e #send emails at beginning, end, suspension or abortion of job
# # Make sure you have given enough time for the copy to complete!
# #$ -m e #send emails at beginning, end, suspension or abortion of job