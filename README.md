# FinalDegreeProject
This repository includes the supplementary data used for the Final Degree Project.


## Pipeline folder.
This folder includes the files needed to run the pipeline until obtaining the MMGBSA input that has to be calculated using another computer cluster.
This Pipeline is able to calcualted the intermediate steps needed to get all the descriptors for each Molecular Dynamic.

* **Execution_Analysis_Results.py**: Python Script that is able to execute all the Analysis for the fragments specified in the arguments.

* **en_Create_Pipeline_Files**: Example file that includes all the arguments for the fragments we want to run the pipeline.

* **Create_Pipeline_Files.py**: Program that is able to create the input files for each fragment we want to run the pipeline using the arguments specificed in the en_Craete_Pipeline_Files file.

* **R1_pipeline_fdMD_OneTrajREM.py**: First program that calculates the main descriptors of each trajectory.

* **R2_pipeline_fdMD_ReactiveTraj.py**: Second program executed by R1_pipeline_fdMD_OneTrajREM.py, that using the data obtained, discards the non-reactive trajectories.

* **R3_pipeline_fdMD_MMPBSA_DYN_OneLig.py**: Third program executed by R2_pipeline_fdMD_ReactiveTraj.py that is able to remove these non-reactive trajectories and create the input files to calculate the MMGBSA (Energy of Binding) of the Reactive Trajectories.

## MD_Analysis folder.

* **en_R1_fdMD_OneTrajREM_jaen**: Example file that runs the R1_fdMD_OneTrajREM.py program and has the arguments needed.

* **R1_fdMD_OneTrajREM.py**: First program run to calcualte the descriptors of each Molecular Dynamic.

* **en_R2_fdMD_ReactiveTraj**: Example file that has the arguments needed to run R2_fdMD_ReactiveTraj.py program and execute it.

* **R2_fdMD_ReactiveTraj.py**: Second program executed to discard the non-reactive trajectories.

* **en_R3_fdMD_MMPBSA_DYN_OneLig**: Example file that has the arguments needed to run R3_fdMD_MMPBSA_DYN_OneLig.py and execute it.

* **R3_fdMD_MMPBSA_DYN_OneLig.py**: Third program needed to create the MMGBSA (Energy of Binding) input for the reactive trajectories.

* **R4_fdMD_MMGBSA_Analize.py**: Fourth program that is able to process all the descriptors and create the average of all of them using the last 10% of the MD.

* **en_Generate_Summary_Excel**: Example file that is able to execute the generation of a Summary Excel file with the four Molecular Dynamics of each fragment.

* **Generate_Summary_Excel.py**: Program that generates the Summary Excel file with the four Molecular Dyanmics of each fragment.

* **Generate_Average_Excel_file.py**: Python program that is able to transform the Summary Excel file into an Average Excel file to know which would be the most suitable Binding site.

**Note that:** After executing R3_fdMD_MMPBSA_DYN_OneLig.py file we need to transfer the files to another computer cluster and execute the MMGBSA for each reactive trajectory alone.
