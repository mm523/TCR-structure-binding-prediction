# Extract-info-from-PDB

These are the scripts that allow to start off with a TCR sequence, model it and extract interface information (https://www.frontiersin.org/articles/10.3389/fphys.2021.730908/full). The script relies on TCRpMHCmodels, which I obtained in command line form from Kamilla Jensen (main author).

The main pipeline scripts are in the **processing_and_prediction** folder and are the following:

* **01_dist_atchley.py** renumbers the PDB to IMGT standards, finds pairwise distances and extracts the sequence information
* **02_energies.py** takes the renumbered structure and calculates the parwise energies
* **03_extract_template_info.py** takes template information from TCRpMHCmodels and saves it into a .csv file.
* **04_prediction_set_as_argument.py** makes the prediction given the input model and saves the output to a .txt file. It takes as argument to the function the subset of features you want to use. 
* **train_model_diff_sets.py** trains the models. It takes as argument to the function the subset of features you want to use.

Note: *at* and *gt* determine if all structures are used in the script, or just the ones with good templates (gt). 
The other files in this folder are (mostly) associated functions.

The prepare_fasta_files folder contains older scripts that I used to create the datasets. 

The parse_results_files go into each log file created for a sequence and obtain the prediction for that sequence.

other contains some other useful scripts. 
