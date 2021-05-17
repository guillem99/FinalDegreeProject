#!/usr/bin/env python3

# Import libraries neeeded.
import numpy as np
import os, sys, stat
from datetime import datetime, date, time
from math import sqrt
import  matplotlib
matplotlib.use('Agg')
import  matplotlib.pyplot as plt
from argparse import ArgumentParser
from os import listdir
from os.path import isfile, join


####################################################################
#####  You need to have R1_*, R2_* and R3_* in this directory  #####
#####                                                          #####
####################################################################

def main() :
  # Call the funciton to print the program information.
  prog_info ()

  # Define a counter to 0 to count the cpu we are suing (if we are in a qsub system).  
  count = 0

  # Read Arguments
  args = cmdlineparse()

  ## General Arguments
  # Read molecule used.
  molecule = args.molecule

  # Split (if we want to iterate for more than moleclue).
  molecules = molecule.split(",")
  
  # Read file R1
  R1 = args.R1
  
  # Read file R2
  R2 = args.R2
  
  # Read file R3
  R3 = args.R3
  
  # Parameters that are different for each different fragment.
  dyn_1 = args.dyn_1
  dyn_1_ = dyn_1.split(",")
  dyn_2 = args.dyn_2
  dyn_2_ = dyn_2.split(",")
  dyn_3 = args.dyn_3
  dyn_3_ = dyn_3.split(",")
  dyn_4 = args.dyn_4
  dyn_4_ = dyn_4.split(",")  
  dir_traj_dyn_1 = args.dir_traj_dyn_1
  dir_traj_dyn_1_ = dir_traj_dyn_1.split(",")
  dir_traj_dyn_2 = args.dir_traj_dyn_2
  dir_traj_dyn_2_ = dir_traj_dyn_2.split(",")
  dir_traj_dyn_3 = args.dir_traj_dyn_3
  dir_traj_dyn_3_ = dir_traj_dyn_3.split(",")
  dir_traj_dyn_4 = args.dir_traj_dyn_4    
  dir_traj_dyn_4_ = dir_traj_dyn_4.split(",")
  top_name = args.top_name
  tops_name = top_name.split(",")
  Atom_Name = args.Atom_Name
  Atom_Names = Atom_Name.split(",")
  prefix_traj = args.prefix_traj  
  prefix_trajs = prefix_traj.split(",")
  fragment = args.fragment
  fragments = fragment.split(",")
  prefix_filpdb = args.prefix_filrep
  prefix_filreps = prefix_filpdb.split(",")
  name_atom = args.name_atomrep
  name_atomreps = name_atom.split(",")
  fil_prot = args.protein_pdb
  
  # Parameters that are the same for each different fragment.
  SuperposeProt = args.SuperposeProt
  ref_pdb = args.ref_pdb
  ini_traj = args.ini_traj
  end_traj = args.end_traj
  ini_read = args.ini_read
  end_read = args.end_read
  inc_read = args.inc_read
  full_traj = args.full_traj
  last_pdb = args.last_pdb
  clean = args.clean
  pocket = args.use_pocket
  xray = args.XRay
  dis_min = args.dis_lig_prot_min
  snaps_byone = args.snaps_by_one
  num_ns_anal = args.num_ns_anal
  percent_anal = args.percent_anal
  num_res_prot = int(args.num_res_prot)
  queue_system = args.queue_system
  dis_plots = args.dis_plots
  dis_plots = dis_plots.split(",") 
  info_from_top = args.info_from_top
  num_plots = args.num_plots

  ## Check if R1, R2 and R3 file exists in the directory.
  # Get the path of our directory. 
  mypath = os.getcwd()

  # List the files that are in our directory.
  mypath_files = [f for f in listdir(mypath) if isfile(join(mypath, f))]
  
  # Check if R1 exists in the list of file of our directory.
  if (R1 not in mypath_files):
    # Print that R1 does not exist.
    print("There is no such file called ", R1)

    # Close.
    sys.exit()
  
  # Check if R2 exists in the list of file of our directory.
  if (R2 not in mypath_files):
    # Print that R2 does not exist.
    print("There is no such file called ", R2)

    # Close.
    sys.exit()
  
  # Check if R3 exists in the list of file of our directory.
  if (R3 not in mypath_files):
    # Print that R1 does not exist.
    print("There is no such file called ", R3)

    # Close.
    sys.exit()
  

  # Iterate all the molecules to copy files R1, R2 and R3 and create the file en_R1..., en_R2..., etc.
  for ii in range(0, len(molecules)):
    # Catch the parameters for molecule ii.
    molecule = molecules[ii]
    top_name = tops_name[ii]
    Atom_Name = Atom_Names[ii]
    prefix_traj = prefix_trajs[ii]
    fragment = fragments[ii]
    prefix_filpdb = prefix_filreps[ii]
    name_atom = name_atomreps[ii]
    dyn_1 = dyn_1_[ii]
    dyn_2 = dyn_2_[ii]
    dyn_3 = dyn_3_[ii]
    dyn_4 = dyn_4_[ii]
    dir_traj_1 = dir_traj_dyn_1_[ii]
    dir_traj_2 = dir_traj_dyn_2_[ii] 
    dir_traj_3 = dir_traj_dyn_3_[ii]
    dir_traj_4 = dir_traj_dyn_4_[ii]
    
    # Print the molecule name we are iterating. 
    print("For Molecule named = " + molecule)

    ## Copy files R1, R2 and R3 in each directory of the molecule. 
    # Change directory to the initial directory.
    os.chdir(mypath)
    
    # Copy files R1, R2 and R3 to dyn_1
    os.system("cp " + R1 + " " + dyn_1)
    os.system("cp " + R2 + " " + dyn_1)
    os.system("cp " + R3 + " " + dyn_1)
   
    # Print that we have copied R1, R2 and R3 in dyn_1.
    print("-- Copied R1, R2 and R3 in " + dyn_1)
    
   # Copy files R1, R2 and R3 to dyn_2.
    os.system("cp " + R1 + " " + dyn_2)
    os.system("cp " + R2 + " " + dyn_2)
    os.system("cp " + R3 + " " + dyn_2)

    # Print that we have copied R1, R2 and R3 in dyn_2.
    print("-- Copied R1, R2 and R3 in " + dyn_2)

    # Copy files R1, R2 and R3 to dyn_3.
    os.system("cp " + R1 + " " + dyn_3)
    os.system("cp " + R2 + " " + dyn_3)
    os.system("cp " + R3 + " " + dyn_3)

    # Print that we have copied R1, R2 and R3 in dyn_3.
    print("-- Copied R1, R2 and R3 in " + dyn_3)

    # Copy files R1, R2 and R3 to dyn_4.
    os.system("cp " + R1 + " " + dyn_4)
    os.system("cp " + R2 + " " + dyn_4)
    os.system("cp " + R3 + " " + dyn_4)

    # Print that we have copied R1, R2 and R3 in dyn_4.
    print("-- Copied R1, R2 and R3 in " + dyn_4)

    # Once copied the R1, R2 and R3 files, we can create the en_R1..., en_R2... and en_R3... files.
    # Create en_R1... in dyn_1.
    # Change directory to the initial one.
    os.chdir(dyn_1)

    ## Create en_R1 for dyn_1
    # Define file input.
    file_inp  = "en_" + R1[:-3]

    # Define file output.
    file_out  = "Info_" + R1[:-3]

    # Open the inpput file.
    file_     = open (file_inp,"w")

    # Add line by line.
    line = "#!/bin/bash -f " + "\n"
    file_.write(line)
    line = "# "+ "\n"
    file_.write(line)
    line = "export DIRI=\"" + dyn_1 + "\"" + "\n"
    file_.write(line)
    line = "cd $DIRI" + "\n" + "\n" + "\n"
    file_.write(line)

    # Avoid these lines depending on the server we are working.
    '''
    line = "source /home/jaime/.bashrc" + "\n"
    file_.write(line)
    line = "source  /home/prog/amber18/miniconda/bin/activate py3" + "\n"
    file_.write(line)
    line = "conda info --envs" + "\n" + "\n"
    file_.write(line)
    line = "export  AMBERHOME=/home/prog/amber18" + "\n"
    file_.write(line)
    line = "source $AMBERHOME/amber.sh" + "\n" + "\n" + "\n"
    file_.write(line)
    '''

    # Get the name of the en_R2 file.
    en_R2 = "en_" + R2[:-3]
    
    # Add the line that includes all the arguments needed for executing R1 file.
    line = "./" + R1 + " -info_from_top " + str(info_from_top) + " -top_name " + str(top_name) + " -Atom_Name " + str(Atom_Name) + " -SuperposeProt " + str(SuperposeProt) + " -ref_pdb " + str(ref_pdb) + " -dir_traj " + str(dir_traj_1) + " -prefix_traj " + str(prefix_traj) + " -ini_read "+ str(ini_read) + " -end_read "+ str(end_read) + " -ini_traj " + str(ini_traj) + " -end_traj " + str(end_traj) + " -inc_read " + str(inc_read) + " -full_traj " + str(full_traj) + " -last_pdb " + str(last_pdb) + " -clean " + str(clean) + " -en_R2 " + str(en_R2) + " -queue_system " + str(queue_system) 
    
    # If the queue system is qsub, add the cpu we want to use.
    if (queue_system == "qsub"):
      # Read the cpus we want to use.
      cpu = args.cpu
      
      # Split the cpus available.
      cpus = cpu.split(",")

      # If the count is greater than the total of cpus (start again).
      if count >= len(cpus):
        # Define count 0 to repeat the cpus.
        count = 0

      # Get the cpu we want in that case.
      cpu = cpus[count]

      # Add one value to the cppu (to use the next one in the next step).
      count += 1

      # Add the final argument in the line added to the en_R1... input file.
      line += " -cpu " + cpu

    # Add the final string in the input line for en_R1... file.
    line += " > " + file_out + "\n"

    # Write the line in the input file.
    file_.write(line)

    # Define the name of the running file.
    run_cpptraj = file_inp

    # Change the mode of the input file to be able to execute it.
    os.system('chmod u+x {} '.format(run_cpptraj) )
 
    # Print which file we have created.
    print("-- Created en_R1 in " + dyn_1)

    # Create en_R1 in dyn_2.
    # Change to dyn_2 directory. 
    os.chdir(dyn_2)
    
    # Define the input file.
    file_inp  = "en_" + R1[:-3]

    # Define the output file.
    file_out  = "Info_" + R1[:-3]

    # Open the input file.
    file_     = open (file_inp,"w")
    
    # Add line by line.
    line = "#!/bin/bash -f " + "\n"
    file_.write(line)
    line = "# "+ "\n"
    file_.write(line)

    '''    
    line = "export DIRI=\"" + dyn_2 + "\"" + "\n"
    file_.write(line)
    line = "cd $DIRI" + "\n" + "\n" + "\n"
    file_.write(line)

    line = "source /home/jaime/.bashrc" + "\n"
    file_.write(line)
    line = "source  /home/prog/amber18/miniconda/bin/activate py3" + "\n"
    file_.write(line)
    line = "conda info --envs" + "\n" + "\n"
    file_.write(line)
    line = "export  AMBERHOME=/home/prog/amber18" + "\n"
    file_.write(line)
    line = "source $AMBERHOME/amber.sh" + "\n" + "\n" + "\n"
    file_.write(line)
    '''   
    # Get the name of the en_R2 file.
    en_R2 = "en_" + R2[:-3]

    # Add the line with all the arguments needed to execute the R1 file.
    line = "./" + R1 + " -info_from_top " + str(info_from_top) + " -top_name " + str(top_name) + " -Atom_Name " + str(Atom_Name) + " -SuperposeProt " + str(SuperposeProt) + " -ref_pdb " + str(ref_pdb) + " -dir_traj " + str(dir_traj_2) + " -prefix_traj " + str(prefix_traj) + " -ini_read "+ str(ini_read) + " -end_read "+ str(end_read) + " -ini_traj " + str(ini_traj) + " -end_traj " + str(end_traj) + " -inc_read " + str(inc_read) + " -full_traj " + str(full_traj) + " -last_pdb " + str(last_pdb) + " -clean " + str(clean) + " -en_R2 " + str(en_R2) + " -queue_system " + str(queue_system) 

    # If the queue system is qsub, add the cpu we want to use.
    if (queue_system == "qsub"):
      # Read the cpus we want to use.
      cpu = args.cpu
      
      # Split the cpus available.
      cpus = cpu.split(",")

      # If the count is greater than the total of cpus (start again).
      if count >= len(cpus):
        # Define count 0 to repeat the cpus.
        count = 0

      # Get the cpu we want in that case.
      cpu = cpus[count]

      # Add one value to the cppu (to use the next one in the next step).
      count += 1

      # Add the final argument in the line added to the en_R1... input file.
      line += " -cpu " + cpu

    # Add the final string in the input line for en_R1... file.
    line += " > " + file_out + "\n"

    # Write the line in the input file.
    file_.write(line)

    # Get the file name we want to run.
    run_cpptraj = file_inp

    # Change the mode of the input file to be able to execute it.
    os.system('chmod u+x {} '.format(run_cpptraj) )

    # Print the name of the en_R1 file.
    print("-- Created en_R1 in " + dyn_2)

    # Creamos en_R1 en dyn_3
    os.chdir(dyn_3)

    ## Create en_R1 for dyn_2
    # Define the input file name.
    file_inp  = "en_" + R1[:-3]

    # Define the output file name.
    file_out  = "Info_" + R1[:-3]

    # Open the input file.
    file_     = open (file_inp,"w")
    
    # Write line by line to the input file.
    line = "#!/bin/bash -f " + "\n"
    file_.write(line)
    line = "# "+ "\n"
    file_.write(line)
    line = "export DIRI=\"" + dyn_3 + "\"" + "\n"
    file_.write(line)
    line = "cd $DIRI" + "\n" + "\n" + "\n"
    file_.write(line)

    # Avoid these lines depending on which server we are workinf with.
    '''
    line = "source /home/jaime/.bashrc" + "\n"
    file_.write(line)
    line = "source  /home/prog/amber18/miniconda/bin/activate py3" + "\n"
    file_.write(line)
    line = "conda info --envs" + "\n" + "\n"
    file_.write(line)
    line = "export  AMBERHOME=/home/prog/amber18" + "\n"
    file_.write(line)
    line = "source $AMBERHOME/amber.sh" + "\n" + "\n" + "\n"
    file_.write(line)
    '''

    # Get the name of the en_R2 file.
    en_R2 = "en_" + R2[:-3]

    # Add the line that includes all the arguments needed for executing R1 file.
    line = "./" + R1 + " -info_from_top " + str(info_from_top) + " -top_name " + str(top_name) + " -Atom_Name " + str(Atom_Name) + " -SuperposeProt " + str(SuperposeProt) + " -ref_pdb " + str(ref_pdb) + " -dir_traj " + str(dir_traj_3) + " -prefix_traj " + str(prefix_traj) + " -ini_read "+ str(ini_read) + " -end_read "+ str(end_read) + " -ini_traj " + str(ini_traj) + " -end_traj " + str(end_traj) + " -inc_read " + str(inc_read) + " -full_traj " + str(full_traj) + " -last_pdb " + str(last_pdb) + " -clean " + str(clean) + " -en_R2 " + str(en_R2) + " -queue_system " + str(queue_system)

    # If the queue system is qsub, add the cpu we want to use.
    if (queue_system == "qsub"):
      # Read the cpus we want to use.
      cpu = args.cpu
      
      # Split the cpus available.
      cpus = cpu.split(",")

      # If the count is greater than the total of cpus (start again).
      if count >= len(cpus):
        # Define count 0 to repeat the cpus.
        count = 0

      # Get the cpu we want in that case.
      cpu = cpus[count]

      # Add one value to the cppu (to use the next one in the next step).
      count += 1

      # Add the final argument in the line added to the en_R1... input file.
      line += " -cpu " + cpu

    # Add the final string in the input line for en_R1... file.
    line += " > " + file_out + "\n"

    # Write the line in the input file.
    file_.write(line)

    # Define the name of the program we want to run. 
    run_cpptraj = file_inp

    # Change the mode of the input file to be able to execute it. 
   os.system('chmod u+x {} '.format(run_cpptraj) )

    # Print the name of the file we have created.
    print("-- Created en_R1 in " + dyn_3)
   
    # Creamos en_R1 en dyn_4
    # Change directory to the dyn_4.
    os.chdir(dyn_4)

    # Define the input file.
    file_inp  = "en_" + R1[:-3]

    # Define the output file.
    file_out  = "Info_" + R1[:-3]

    # Open the input file.
    file_     = open (file_inp,"w")
    
    # Write line by line.
    line = "#!/bin/bash -f " + "\n"
    file_.write(line)
    line = "# "+ "\n"
    file_.write(line)
    line = "export DIRI=\"" + dyn_4 + "\"" + "\n"
    file_.write(line)
    line = "cd $DIRI" + "\n" + "\n" + "\n"
    file_.write(line)

    # Avoid these lines depending on the server we are working.
    '''
    line = "source /home/jaime/.bashrc" + "\n"
    file_.write(line)
    line = "source  /home/prog/amber18/miniconda/bin/activate py3" + "\n"
    file_.write(line)
    line = "conda info --envs" + "\n" + "\n"
    file_.write(line)
    line = "export  AMBERHOME=/home/prog/amber18" + "\n"
    file_.write(line)
    line = "source $AMBERHOME/amber.sh" + "\n" + "\n" + "\n"
    file_.write(line)
    '''    

    # Get the name of the en_R2 file.
    en_R2 = "en_" + R2[:-3]

    # Define the line that allow to run the R1 program.
    line = "./" + R1 + " -info_from_top " + str(info_from_top) + " -top_name " + str(top_name) + " -Atom_Name " + str(Atom_Name) + " -SuperposeProt " + str(SuperposeProt) + " -ref_pdb " + str(ref_pdb) + " -dir_traj " + str(dir_traj_4) + " -prefix_traj " + str(prefix_traj) + " -ini_read "+ str(ini_read) + " -end_read "+ str(end_read) + " -ini_traj " + str(ini_traj) + " -end_traj " + str(end_traj) + " -inc_read " + str(inc_read) + " -full_traj " + str(full_traj) + " -last_pdb " + str(last_pdb) + " -clean " + str(clean) + " -en_R2 " + str(en_R2) + " -queue_system " + str(queue_system)

    # If the queue system is qsub, add the cpu we want to use.
    if (queue_system == "qsub"):
      # Read the cpus we want to use.
      cpu = args.cpu
      
      # Split the cpus available.
      cpus = cpu.split(",")

      # If the count is greater than the total of cpus (start again).
      if count >= len(cpus):
        # Define count 0 to repeat the cpus.
        count = 0

      # Get the cpu we want in that case.
      cpu = cpus[count]

      # Add one value to the cppu (to use the next one in the next step).
      count += 1

      # Add the final argument in the line added to the en_R1... input file.
      line += " -cpu " + cpu

    # Add the final string in the input line for en_R1... file.
    line += " > " + file_out + "\n"

    # Write the line in the input file.
    file_.write(line)

    # Define the file we want to run.
    run_cpptraj = file_inp

    # Change the mode of the input file to be able to execute it.
    os.system('chmod u+x {} '.format(run_cpptraj) )

    # Print the file we have created.
    print("-- Created en_R1 in " + dyn_4)

    # Change directory to our initial directory.
    os.chdir(mypath)
  
    ## Create the executable file R2_* (en_R2_) in dyn_1
    # Change directory to the initial one.
    os.chdir(dyn_1)

    # Define th einput file.
    file_inp  = "en_" + R2[:-3]

    # Define the output file.
    file_out  = "Info_" + R2[:-3]

    # Open the input file.
    file_     = open (file_inp,"w")

    # Write line by line into the input file.
    line = "#!/bin/bash -f " + "\n"
    file_.write(line)
    line = "# "+ "\n"
    file_.write(line)
    line = "export DIRI=\"" + dyn_1 + "\"" + "\n"
    file_.write(line)
    line = "cd $DIRI" + "\n" + "\n" + "\n"
    file_.write(line)
  
    # Avoid these lines depending on the server we are using.
    '''
    line = "source /home/jaime/.bashrc" + "\n"
    file_.write(line)
    line = "source  /home/prog/amber18/miniconda/bin/activate py3" + "\n"
    file_.write(line)
    line = "conda info --envs" + "\n" + "\n"
    file_.write(line)
    line = "export  AMBERHOME=/home/prog/amber18" + "\n"
    file_.write(line)
    line = "source $AMBERHOME/amber.sh" + "\n" + "\n" + "\n"
    file_.write(line)
    '''   
    
    # Get the name of the en_R3 file. 
    en_R3  = "en_" + R3[:-3]

    # Define the line that allow us to run the R2 program.
    line = "./" + R2 + " -end_traj " + str(end_traj) + " -top_name " + str(top_name) + " -pocket " + str(pocket) + " -num_res_prot " + str(num_res_prot) + " -xray " + str(xray) + " -fil_prot " + str(fil_prot) + " -dis_min " + str(dis_min) + " -dir_traj " + str(dir_traj_1) + " -prefix_filpdb " + str(prefix_filpdb) + " -name_atom " + str(name_atom) + " -num_plots " + str(num_plots) + " -dis_plots " + str(dis_plots[0]) + "," + str(dis_plots[1]) + " -snaps_byone " + str(snaps_byone) + " -num_ns_anal " + str(num_ns_anal) + " -percent_anal " + str(percent_anal) + " -en_R3 " + str(en_R3) + " -queue_system " + str(queue_system) + " -fragment " + str(fragment)

    # If the queue system is qsub, add the cpu we want to use.
    if (queue_system == "qsub"):
      # Read the cpus we want to use.
      cpu = args.cpu
      
      # Split the cpus available.
      cpus = cpu.split(",")

      # If the count is greater than the total of cpus (start again).
      if count >= len(cpus):
        # Define count 0 to repeat the cpus.
        count = 0

      # Get the cpu we want in that case.
      cpu = cpus[count]

      # Add one value to the cppu (to use the next one in the next step).
      count += 1

      # Add the final argument in the line added to the en_R1... input file.
      line += " -cpu " + cpu

    # Add the final string in the input line for en_R1... file.
    line += " > " + file_out + "\n"

    # Write the line in the input file.
    file_.write(line)

    # Define the file we want to run.
    run_cpptraj = file_inp

    # Change the mode of the input file to be able to execute it.
    os.system('chmod u+x {} '.format(run_cpptraj) )

    # Print the line we haev created.
    print("-- Created en_R2 in " + dyn_1)


    ## Create the executable file R2_* (en_R2_) in dyn_2
    # Change directory to the dyn_2.
    os.chdir(dyn_2)

    # Define th einput file.
    file_inp  = "en_" + R2[:-3]

    # Define the output file.
    file_out  = "Info_" + R2[:-3]

    # open the input file.
    file_     = open (file_inp,"w")

    # Write line by line in the input line.
    line = "#!/bin/bash -f " + "\n"
    file_.write(line)
    line = "# "+ "\n"
    file_.write(line)
    line = "export DIRI=\"" + dyn_2 + "\"" + "\n"
    file_.write(line)
    line = "cd $DIRI" + "\n" + "\n" + "\n"
    file_.write(line)

    # Avoid these lines depending on the server you are working with.
    '''
    line = "source /home/jaime/.bashrc" + "\n"
    file_.write(line)
    line = "source  /home/prog/amber18/miniconda/bin/activate py3" + "\n"
    file_.write(line)
    line = "conda info --envs" + "\n" + "\n"
    file_.write(line)
    line = "export  AMBERHOME=/home/prog/amber18" + "\n"
    file_.write(line)
    line = "source $AMBERHOME/amber.sh" + "\n" + "\n" + "\n"
    file_.write(line)
    '''    

    # Get the name of the en_R3 file.
    en_R3  = "en_" + R3[:-3]

    # Add the line that allow us to run R2.
    line = "./" + R2 + " -end_traj " + str(end_traj) + " -top_name " + str(top_name) + " -pocket " + str(pocket) + " -num_res_prot " + str(num_res_prot) + " -xray " + str(xray) + " -fil_prot " + str(fil_prot) + " -dis_min " + str(dis_min) + " -dir_traj " + str(dir_traj_2) + " -prefix_filpdb " + str(prefix_filpdb) + " -name_atom " + str(name_atom) + " -num_plots " + str(num_plots) + " -dis_plots " + str(dis_plots[0]) + "," + str(dis_plots[1]) + " -snaps_byone " + str(snaps_byone) + " -num_ns_anal " + str(num_ns_anal) + " -percent_anal " + str(percent_anal) + " -en_R3 " + str(en_R3) + " -queue_system " + str(queue_system) + " -fragment " + str(fragment) 

    # If the queue system is qsub, add the cpu we want to use.
    if (queue_system == "qsub"):
      # Read the cpus we want to use.
      cpu = args.cpu
      
      # Split the cpus available.
      cpus = cpu.split(",")

      # If the count is greater than the total of cpus (start again).
      if count >= len(cpus):
        # Define count 0 to repeat the cpus.
        count = 0

      # Get the cpu we want in that case.
      cpu = cpus[count]

      # Add one value to the cppu (to use the next one in the next step).
      count += 1

      # Add the final argument in the line added to the en_R1... input file.
      line += " -cpu " + cpu

    # Add the final string in the input line for en_R1... file.
    line += " > " + file_out + "\n"

    # Write the line in the input file.
    file_.write(line)

    # Define the line to run the input file.
    run_cpptraj = file_inp

    # Change the mode of the input file to be able to execute it.
    os.system('chmod u+x {} '.format(run_cpptraj) )

    # Print which file we have created.
    print("-- Created en_R2 in " + dyn_2)

    ## Create the executable file R2_* (en_R2_) in dyn_3
    # Change the directory to dyn_3
    os.chdir(dyn_3)

    # Define the input file.
    file_inp  = "en_" + R2[:-3]

    # Define the output file.
    file_out  = "Info_" + R2[:-3]

    # Open the input file.
    file_     = open (file_inp,"w")

    # Write line by line.
    line = "#!/bin/bash -f " + "\n"
    file_.write(line)
    line = "# "+ "\n"
    file_.write(line)
    line = "export DIRI=\"" + dyn_3 + "\"" + "\n"
    file_.write(line)
    line = "cd $DIRI" + "\n" + "\n" + "\n"
    file_.write(line)

    # Avoid these lines depending on the server we are working.
    '''
    line = "source /home/jaime/.bashrc" + "\n"
    file_.write(line)
    line = "source  /home/prog/amber18/miniconda/bin/activate py3" + "\n"
    file_.write(line)
    line = "conda info --envs" + "\n" + "\n"
    file_.write(line)
    line = "export  AMBERHOME=/home/prog/amber18" + "\n"
    file_.write(line)
    line = "source $AMBERHOME/amber.sh" + "\n" + "\n" + "\n"
    file_.write(line)
    '''    
    
    # Get the name of the en_R3 file.
    en_R3  = "en_" + R3[:-3]

    # Define the line that allow us to run R2 file.
    line = "./" + R2 + " -end_traj " + str(end_traj) + " -top_name " + str(top_name) + " -pocket " + str(pocket) + " -num_res_prot " + str(num_res_prot) + " -xray " + str(xray) + " -fil_prot " + str(fil_prot) + " -dis_min " + str(dis_min) + " -dir_traj " + str(dir_traj_3) + " -prefix_filpdb " + str(prefix_filpdb) + " -name_atom " + str(name_atom) + " -num_plots " + str(num_plots) + " -dis_plots " + str(dis_plots[0]) + "," + str(dis_plots[1]) + " -snaps_byone " + str(snaps_byone) + " -num_ns_anal " + str(num_ns_anal) + " -percent_anal " + str(percent_anal) + " -en_R3 " + str(en_R3) + " -queue_system " + str(queue_system) + " -fragment " + str(fragment) 

    # If the queue system is qsub, add the cpu we want to use.
    if (queue_system == "qsub"):
      # Read the cpus we want to use.
      cpu = args.cpu
      
      # Split the cpus available.
      cpus = cpu.split(",")

      # If the count is greater than the total of cpus (start again).
      if count >= len(cpus):
        # Define count 0 to repeat the cpus.
        count = 0

      # Get the cpu we want in that case.
      cpu = cpus[count]

      # Add one value to the cppu (to use the next one in the next step).
      count += 1

      # Add the final argument in the line added to the en_R1... input file.
      line += " -cpu " + cpu

    # Add the final string in the input line for en_R1... file.
    line += " > " + file_out + "\n"

    # Write the line in the input file.
    file_.write(line)

    # Define the variable taht allow us to execute the inpput file.
    run_cpptraj = file_inp

    # Change the mode of the input file to be able to execute it.
    os.system('chmod u+x {} '.format(run_cpptraj) )

    # Print the file we have created.
    print("-- Created en_R2 in " + dyn_3)

    ## Create the executable file R2_* (en_R2_) in dyn_4
    # Change to the directory dyn_4.
    os.chdir(dyn_4)

    # DEfine the input file.
    file_inp  = "en_" + R2[:-3]

    # Define the output file.
    file_out  = "Info_" + R2[:-3]

    # Open the input file.
    file_     = open (file_inp,"w")

    # Write line by line.
    line = "#!/bin/bash -f " + "\n"
    file_.write(line)
    line = "# "+ "\n"
    file_.write(line)
    line = "export DIRI=\"" + dyn_4 + "\"" + "\n"
    file_.write(line)
    line = "cd $DIRI" + "\n" + "\n" + "\n"
    file_.write(line)

    .#Avoid these lines depeneding on the server we are.
    '''
    line = "source /home/jaime/.bashrc" + "\n"
    file_.write(line)
    line = "source  /home/prog/amber18/miniconda/bin/activate py3" + "\n"
    file_.write(line)
    line = "conda info --envs" + "\n" + "\n"
    file_.write(line)
    line = "export  AMBERHOME=/home/prog/amber18" + "\n"
    file_.write(line)
    line = "source $AMBERHOME/amber.sh" + "\n" + "\n" + "\n"
    file_.write(line)
    '''    

    # Get the name of the en_R3 file.
    en_R3  = "en_" + R3[:-3]

    # Add the line that allow us to run R2 file.
    line = "./" + R2 + " -end_traj " + str(end_traj) + " -top_name " + str(top_name) + " -pocket " + str(pocket) + " -num_res_prot " + str(num_res_prot) + " -xray " + str(xray) + " -fil_prot " + str(fil_prot) + " -dis_min " + str(dis_min) + " -dir_traj " + str(dir_traj_4) + " -prefix_filpdb " + str(prefix_filpdb) + " -name_atom " + str(name_atom) + " -num_plots " + str(num_plots) + " -dis_plots " + str(dis_plots[0]) + "," + str(dis_plots[1]) + " -snaps_byone " + str(snaps_byone) + " -num_ns_anal " + str(num_ns_anal) + " -percent_anal " + str(percent_anal) + " -en_R3 " + str(en_R3) + " -queue_system " + str(queue_system) + " -fragment " + str(fragment)

    # If the queue system is qsub, add the cpu we want to use.
    if (queue_system == "qsub"):
      # Read the cpus we want to use.
      cpu = args.cpu
      
      # Split the cpus available.
      cpus = cpu.split(",")

      # If the count is greater than the total of cpus (start again).
      if count >= len(cpus):
        # Define count 0 to repeat the cpus.
        count = 0

      # Get the cpu we want in that case.
      cpu = cpus[count]

      # Add one value to the cppu (to use the next one in the next step).
      count += 1

      # Add the final argument in the line added to the en_R1... input file.
      line += " -cpu " + cpu

    # Add the final string in the input line for en_R1... file.
    line += " > " + file_out + "\n"

    # Write the line in the input file.
    file_.write(line)

    # Define the variable that allow us to run the input file.
    run_cpptraj = file_inp

    # Change the mode of the input file to be able to execute it.
    os.system('chmod u+x {} '.format(run_cpptraj) )

    # Print the file we have created.
    print("-- Created en_R2 in " + dyn_4)


    # Create en_R3... in dyn_1
    # Change directory to dyn_1.
    os.chdir(dyn_1)

    # Define the input file.
    file_inp  = "en_" + R3[:-3]

    # Define the output file.
    file_out  = "Info_" + R3[:-3]

    # Open the input file.
    file_     = open (file_inp,"w")

    # Write line by line.
    line = "#!/bin/bash -f " + "\n"
    file_.write(line)
    line = "# "+ "\n"
    file_.write(line)
    line = "export DIRI=\"" + dyn_1 + "\"" + "\n"
    file_.write(line)
    line = "cd $DIRI" + "\n" + "\n" + "\n"
    file_.write(line)
    
    # Avoid these lines depending on the server we are working.
    '''
    line = "source /home/jaime/.bashrc" + "\n"
    file_.write(line)
    line = "source  /home/prog/amber18/miniconda/bin/activate py3" + "\n"
    file_.write(line)
    line = "conda info --envs" + "\n" + "\n"
    file_.write(line)
    line = "export  AMBERHOME=/home/prog/amber18" + "\n"
    file_.write(line)
    line = "source $AMBERHOME/amber.sh" + "\n" + "\n" + "\n"
    file_.write(line)
    '''
    # Add the line that allows the execution of R3.
    line = "./" + R3 + " > " + file_out + "\n"
 
    # Write line.
    file_.write(line)

    # Define the name of the input file to run.
    run_cpptraj = file_inp

    # Change the mode of the input file to be able to execute it.
    os.system('chmod u+x {} '.format(run_cpptraj) )

    # Print the file we have created.
    print("-- Created en_R3 in " + dyn_1)


    # Create en_R3... in dyn_2
    # Change directory to dyn_2.
    os.chdir(dyn_2)

    # Define the input file.
    file_inp  = "en_" + R3[:-3]

    # Define the output file.
    file_out  = "Info_" + R3[:-3]

    # Open the input file.
    file_     = open (file_inp,"w")

    # Write line by line.
    line = "#!/bin/bash -f " + "\n"
    file_.write(line)
    line = "# "+ "\n"
    file_.write(line)
    line = "export DIRI=\"" + dyn_2 + "\"" + "\n"
    file_.write(line)
    line = "cd $DIRI" + "\n" + "\n" + "\n"
    file_.write(line)

    # Avoid these lines depending on the server we are working.
    '''    
    line = "source /home/jaime/.bashrc" + "\n"
    file_.write(line)
    line = "source  /home/prog/amber18/miniconda/bin/activate py3" + "\n"
    file_.write(line)
    line = "conda info --envs" + "\n" + "\n"
    file_.write(line)
    line = "export  AMBERHOME=/home/prog/amber18" + "\n"
    file_.write(line)
    line = "source $AMBERHOME/amber.sh" + "\n" + "\n" + "\n"
    file_.write(line)
    '''

    # Add the line that allows the execution of R3.
    line = "./" + R3 + " > " + file_out + "\n"

    # Write line.
    file_.write(line)

    # Define the name of the input file to run.
    run_cpptraj = file_inp

    # Change the mode of the input file to be able to execute it.
    os.system('chmod u+x {} '.format(run_cpptraj) )

    # Print the file we have created.
    print("-- Created en_R3 in " + dyn_2)


    # Create en_R3... in dyn_3
    # Change directory to dyn_3.
    os.chdir(dyn_3)

    # Define the input file.
    file_inp  = "en_" + R3[:-3]

    # Define the output file.
    file_out  = "Info_" + R3[:-3]

    # Open the input file.
    file_     = open (file_inp,"w")

    # Write line by line.
    line = "#!/bin/bash -f " + "\n"
    file_.write(line)
    line = "# "+ "\n"
    file_.write(line)
    line = "export DIRI=\"" + dyn_3 + "\"" + "\n"
    file_.write(line)
    line = "cd $DIRI" + "\n" + "\n" + "\n"
    file_.write(line)

    # Avoid these lines depending on the server we are working.
    '''
    line = "source /home/jaime/.bashrc" + "\n"
    file_.write(line)
    line = "source  /home/prog/amber18/miniconda/bin/activate py3" + "\n"
    file_.write(line)
    line = "conda info --envs" + "\n" + "\n"
    file_.write(line)
    line = "export  AMBERHOME=/home/prog/amber18" + "\n"
    file_.write(line)
    line = "source $AMBERHOME/amber.sh" + "\n" + "\n" + "\n"
    file_.write(line)
    '''

    # Add the line that allows the execution of R3.
    line = "./" + R3 + " > " + file_out + "\n"

    # Write line.
    file_.write(line)

    # Define the name of the input file to run.
    run_cpptraj = file_inp

    # Change the mode of the input file to be able to execute it.
    os.system('chmod u+x {} '.format(run_cpptraj) )

    # Print the file we have created.
    print("-- Created en_R3 in " + dyn_3)


    # Create en_R3... in dyn_4
    # Change directory to dyn_4.
    os.chdir(dyn_4)

    # Define the input file.
    file_inp  = "en_" + R3[:-3]

    # Define the output file.
    file_out  = "Info_" + R3[:-3]

    # Open the input file.
    file_     = open (file_inp,"w")

    # Write line by line.
    line = "#!/bin/bash -f " + "\n"
    file_.write(line)
    line = "# "+ "\n"
    file_.write(line)
    line = "export DIRI=\"" + dyn_4 + "\"" + "\n"
    file_.write(line)
    line = "cd $DIRI" + "\n" + "\n" + "\n"
    file_.write(line)

    # Avoid these lines depending on the server we are working.
    '''
    line = "source /home/jaime/.bashrc" + "\n"
    file_.write(line)
    line = "source  /home/prog/amber18/miniconda/bin/activate py3" + "\n"
    file_.write(line)
    line = "conda info --envs" + "\n" + "\n"
    file_.write(line)
    line = "export  AMBERHOME=/home/prog/amber18" + "\n"
    file_.write(line)
    line = "source $AMBERHOME/amber.sh" + "\n" + "\n" + "\n"
    file_.write(line)
    '''

    # Add the line that allows the execution of R3.
    line = "./" + R3 + " > " + file_out + "\n"

    # Write line.
    file_.write(line)

    # Define the name of the input file to run.
    run_cpptraj = file_inp

    # Change the mode of the input file to be able to execute it.
    os.system('chmod u+x {} '.format(run_cpptraj) )

    # Print the file we have created.
    print("-- Created en_R3 in " + dyn_4)

  # Change the directory to the initial one.
  os.chdir(mypath)
  
 
  # Return the main function.
  return 



def cmdlineparse():
  # Arguments used in the Command line to get parameters. 
  parser = ArgumentParser(description="command line arguments")
  
  parser.add_argument("-molecule"    , dest="molecule" , required=True, help="Molecule name we are using as a ligand in that run" )
  parser.add_argument("-R1"    , dest="R1" , required=True, help="Version of the R1_ file we want to execute" )
  parser.add_argument("-R2"    , dest="R2" , required=True, help="Version of the R2_ file we want to execute" )
  parser.add_argument("-R3"    , dest="R3" , required=True, help="Version of the R3_ file we want to execute" )
  parser.add_argument("-dyn_1"    , dest="dyn_1" , required=True, help="Directory where is located the data of dyn_1" )
  parser.add_argument("-dyn_2"    , dest="dyn_2" , required=True, help="Directory where is located the data of dyn_2" )
  parser.add_argument("-dyn_3"    , dest="dyn_3" , required=True, help="Directory where is located the data of dyn_3" )
  parser.add_argument("-dyn_4"    , dest="dyn_4" , required=True, help="Directory where is located the data of dyn_4" )
  parser.add_argument("-queue_system"    , dest="queue_system" , required=True, help="Qhich system we are using" )
  parser.add_argument("-cpu"    ,    dest="cpu"          , required=False, help = " CPu we want to use to execute the files ")
  parser.add_argument("-dir_traj_dyn_1" ,    dest="dir_traj_dyn_1"      , required=True, help=" DIR witn the trajectories  ")
  parser.add_argument("-dir_traj_dyn_2" ,    dest="dir_traj_dyn_2"      , required=True, help=" DIR witn the trajectories  ")
  parser.add_argument("-dir_traj_dyn_3" ,    dest="dir_traj_dyn_3"      , required=True, help=" DIR witn the trajectories  ")
  parser.add_argument("-dir_traj_dyn_4" ,    dest="dir_traj_dyn_4"      , required=True, help=" DIR witn the trajectories  ")
   
  # Arguments for R1 
  parser.add_argument("-top_name"     , dest="top_name"     , required=True, help=" Name of the Topological file ")
  parser.add_argument("-info_from_top", dest="info_from_top", required=True, help=" TRUE : data from TOP    file ")
  parser.add_argument("-num_res_prot" , dest="num_res_prot" , required=False, help=" Number of Protein Residues ")
  parser.add_argument("-Atom_Name"    , dest="Atom_Name"    , required=True,  help="  Atom to analyse (C99) ")
  parser.add_argument("-SuperposeProt", dest="SuperposeProt", required=True,  help="first or reference or xray  ")
  parser.add_argument("-ref_pdb" ,      dest="ref_pdb"      , required=True,  help=" file with the reference = Protein+Ligand ")
  parser.add_argument("-prefix_traj" , dest="prefix_traj"   , required=True, help=" Prefix of the trajectories ")
  parser.add_argument("-ini_traj"  ,   dest="ini_traj"      , required=True, help=" INI for reading trajectories   ")
  parser.add_argument("-end_traj"  ,   dest="end_traj"      , required=True, help=" IFI for reading trajectories   ")
  parser.add_argument("-ini_read" ,    dest="ini_read"      , required=True, help=" INI for reading trajectories ")
  parser.add_argument("-inc_read" ,    dest="inc_read"      , required=True, help=" INC for reading trajectories ")
  parser.add_argument("-end_read" ,    dest="end_read"      , required=True, help=" IFI  ")

  parser.add_argument("-full_traj"      , dest="full_traj"       , required=True , help=" Reading the full traj ( yes)  ")
  parser.add_argument("-num_snaps_traj" , dest="num_snaps_traj"  , required=False, help=" Only these snapshots  ")
  parser.add_argument("-last_pdb"       , dest="last_pdb"        , required=True , help=" Generate a PDB with the las structure ")
  parser.add_argument("-num_snaps_total", dest="num_snaps_total" , required=False, help=" Total number of snapshots ")
  parser.add_argument("-clean"    ,    dest="clean"          , required=False, help=" Remove all intermediate files ")



  # Arguments of R2
  parser.add_argument("-pocket"    , dest="use_pocket" , required=True, help=" Use pocket to define distances (False) " )
  parser.add_argument("-fil_pocket", dest="file_pocket", required=False,help=" File containing the aa of the pocket ")

  parser.add_argument("-xray"      , dest="XRay"       , required=True, help=" Use XRay structure as reference (False) " )
  parser.add_argument("-fil_prot"  , dest="protein_pdb", required=True, help=" Reference Protein structure ")

  parser.add_argument("-dis_min"   , dest="dis_lig_prot_min", required=True, help=" Minimum distance Ligand-Prot ")

  parser.add_argument("-prefix_filpdb", dest="prefix_filrep", required=True, help=" Posfix added to pdb files ")
  parser.add_argument("-name_atom"    , dest="name_atomrep" , required=True, help=" Atom to analyse (C99) ")

  parser.add_argument("-num_plots"    , dest="num_plots"  , required=True, help=" Num Plots to be done (2)  ")
  parser.add_argument("-dis_plots"    , dest="dis_plots"  , required=True, help=" Distances to be Plot  [5,25] ")

  parser.add_argument("-snaps_byone"  , dest="snaps_by_one"   , required=True, help=" Snaps by nanosecond ")
  parser.add_argument("-num_ns_anal"  , dest="num_ns_anal"    , required=True, help=" Num ns to analyse   ")
  parser.add_argument("-percent_anal" , dest="percent_anal"   , required=True, help=" PerCent out of the limits ")
  parser.add_argument("-fragment" , dest="fragment" , required=False, help=" Name of the fragment used ")


  args=parser.parse_args()
  
  return args

# Fucntion to priunt the initial information of the program.
def prog_info () :
  t = datetime.now()
  format_date = "%d-%m-%Y %H:%M:%S"
  ti= t.strftime(format_date)
  print ( "... ... ... ... ... ... ... " )
  print ( "...  Pipeline Execution ... " )
  print ( "...     ( 2020 )        ... " )
  print ( "... ... ... ... ... ... ... " )
  print ( " ", ti)
  print ( "... ... ... ... ... ... " )
  print ( " " )

if __name__ == '__main__':

  # Call the main function to do execute all the process.
  main()



