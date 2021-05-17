#!/usr/bin/env python3
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

  prog_info ()
  
  count = 1
  
  # Read Arguments
  args = cmdlineparse()

  # General Arguments
  molecule = args.molecule
  molecules = molecule.split(",")
  
  R1 = args.R1
  R2 = args.R2
  R3 = args.R3

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
  
  cpus = args.cpus
  cpus = cpus.split(",")
  count_cpu = 0 
  owner_execution = args.owner_execution
  last_program_output = args.last_program_output
  which_script = args.which_script
   
  
  
  # Check if files exists
  mypath = os.getcwd()
  mypath_files = [f for f in listdir(mypath) if isfile(join(mypath, f))]
  if (R1 not in mypath_files):
    print("There is no such file called ", R1)
    sys.exit()
  if (R2 not in mypath_files):
    print("There is no such file called ", R2)
    sys.exit()
  if (R3 not in mypath_files):
    print("There is no such file called ", R3)
    sys.exit()
   

  ii = 0
  molecule = molecules[ii]
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
  top_name = tops_name[ii]
  Atom_Name = Atom_Names[ii]
  prefix_traj = prefix_trajs[ii]
  info_from_top = args.info_from_top
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
  fragment = fragments[ii]
  prefix_filpdb = prefix_filreps[ii]
  name_atom = name_atomreps[ii]    
  pocket = args.use_pocket
  xray = args.XRay
  fils_prot = args.protein_pdb
  fils_prot = fils_prot.split(",")
  fil_prot = fils_prot[ii] 
  dis_min = args.dis_lig_prot_min
  num_plots = args.num_plots
  dis_plots = args.dis_plots
  dis_plots = dis_plots.split(",")
  snaps_byone = args.snaps_by_one
  num_ns_anal = args.num_ns_anal
  percent_anal = args.percent_anal
  num_res_prot = int(args.num_res_prot)
    
  ## Put all dyns together
  dyns = list()
  dir_trajs = list()
  for i in range(0, len(molecules)):
    dyn_1 = dyn_1_[i]
    dyn_2 = dyn_2_[i]
    dyn_3 = dyn_3_[i]
    dyn_4 = dyn_4_[i]
    dyns.append(dyn_1)
    dyns.append(dyn_2)
    dyns.append(dyn_3)
    dyns.append(dyn_4)
    dir_traj_dyn_1 = dir_traj_dyn_1_[i]
    dir_traj_dyn_2 = dir_traj_dyn_2_[i]
    dir_traj_dyn_3 = dir_traj_dyn_3_[i]
    dir_traj_dyn_4 = dir_traj_dyn_4_[i]
    dir_trajs.append(dir_traj_dyn_1)
    dir_trajs.append(dir_traj_dyn_2)
    dir_trajs.append(dir_traj_dyn_3)
    dir_trajs.append(dir_traj_dyn_4)

  ## Get the CPU's that have some running program
  os.chdir(mypath)
  os.system("qstat > data.txt")
  jobs_running = dict()
  with open("data.txt", "r") as f:
    for line in f:
      t = line.split()   
      if (len(t) == 1 and t[0].startswith("en_R")):
        line = f.readline()
        new_line = line.split()        
        if (new_line[1] == owner_execution and new_line[3] == "R"):
          jobs_running[new_line[2]] = t[0]
  
  #os.system("rm data.txt")
  cpus_used = list(jobs_running)
  
  for c in cpus_used:
      if c in cpus:
          cpus.remove(c)

  ## cpus_available are the ones we can execute
  cpus_available = cpus  
  
  
  
  # Si es la primera vez que lo ejecutamos:
  # el last_program_output estara vacio y lanzamos a ejecutar los primeros programas
  ii = 0
  dyns_removed = list()
  
  if (last_program_output == "None"):
    for c in cpus_available:
      if (len(dyns_removed)%4 == 0) and len(dyns_removed) != 0:
          ii += 1
      molecule = molecules[ii]
      top_name = tops_name[ii]
      Atom_Name = Atom_Names[ii]
      prefix_traj = prefix_trajs[ii]
      fragment = fragments[ii]
      prefix_filpdb = prefix_filreps[ii]
      name_atom = name_atomreps[ii]
      fil_prot = fils_prot[ii]
  
      dyns_ = dyns
      
      for dyn in dyns_:
        if dyn == "None":
          print(str(count) + ". I have executed en_R1 file in this directory: " + dyn)
          count += 1
          continue
        # which dir_traj we are using??
        which_dyn = dyn[-5:]
        which_mol = prefix_traj[-3:]
        dir_traj = next((s for s in dir_trajs if (which_dyn in s and which_mol in s)), None)

        
        # Move the files R1, R2 and R3 to each directory
        os.system("cp " + R1 + " " + dyn)
        os.system("cp " + R2 + " " + dyn)
        os.system("cp " + R3 + " " + dyn)
        os.chdir(dyn)  
        
        ## Create en_R1 (R1_...)        
        file_inp  = "en_" + R1[:-3] 
        file_out  = "Info_" + R1[:-3]
        file_     = open (file_inp,"w")

        line = "#!/bin/bash -f " + "\n"
        file_.write(line)
        line = "# "+ "\n"
        file_.write(line)
        line = "export DIRI=\"" + dyn + "\"" + "\n"
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

        en_R2 = "en_" + R2[:-3]
        
        line = "./" + R1 + " -info_from_top " + str(info_from_top) + " -top_name " + str(top_name) + " -Atom_Name " + str(Atom_Name) + " -SuperposeProt " + str(SuperposeProt) + " -ref_pdb " + str(ref_pdb) + " -dir_traj " + str(dir_traj) + " -prefix_traj " + str(prefix_traj) + " -ini_read "+ str(ini_read) + " -end_read "+ str(end_read) + " -ini_traj " + str(ini_traj) + " -end_traj " + str(end_traj) + " -inc_read " + str(inc_read) + " -full_traj " + str(full_traj) + " -last_pdb " + str(last_pdb) + " -clean " + str(clean) + " -en_R2 " + str(en_R2) + " -cpu " + str(c) + " > " + file_out + "\n"
        file_.write(line)

        
        run_cpptraj = file_inp
        os.system('chmod u+x {} '.format(run_cpptraj) )
        
        print(str(count) + ". I have created R1 and en_R1 file in this directory: " + dyn)
        print()
        count += 1

        ## Create the executable file R2_* (en_R2_) 
        file_inp  = "en_" + R2[:-3]
        file_out  = "Info_" + R2[:-3]
        file_     = open (file_inp,"w")

        line = "#!/bin/bash -f " + "\n"
        file_.write(line)
        line = "# "+ "\n"
        file_.write(line)
        line = "export DIRI=\"" + dyn + "\"" + "\n"
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

        en_R3  = "en_" + R3[:-3]
        
        line = "./" + R2 + " -top_name " + str(top_name) + " -pocket " + str(pocket) + " -num_res_prot " + str(num_res_prot) + " -xray " + str(xray) + " -fil_prot " + str(fil_prot) + " -dis_min " + str(dis_min) + " -dir_traj " + str(dir_traj) + " -prefix_filpdb " + str(prefix_filpdb) + " -name_atom " + str(name_atom) + " -num_plots " + str(num_plots) + " -dis_plots " + str(dis_plots[0]) + "," + str(dis_plots[1]) + " -snaps_byone " + str(snaps_byone) + " -num_ns_anal " + str(num_ns_anal) + " -percent_anal " + str(percent_anal) + " -en_R3 " + str(en_R3) + " -cpu " + str(c) + " -fragment " + str(fragment) + " > " + file_out + "\n"
        file_.write(line)


        run_cpptraj = file_inp
        os.system('chmod u+x {} '.format(run_cpptraj) )


        print(str(count) + ". I have created R2 and en_R2 file in this directory: " + dyn)
        print()
        count += 1
        
        # Create the executable file R3
        file_inp  = "en_" + R3[:-3]
        file_out  = "Info_" + R3[:-3]
        file_     = open (file_inp,"w")

        line = "#!/bin/bash -f " + "\n"
        file_.write(line)
        line = "# "+ "\n"
        file_.write(line)
        line = "export DIRI=\"" + dyn + "\"" + "\n"
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


        line = "./" + R3 + " > " + file_out + "\n"
        file_.write(line)


        run_cpptraj = file_inp
        os.system('chmod u+x {} '.format(run_cpptraj) )
        # Don't run the R3 file
        #os.system(run_cpptraj)
        
        print(str(count) + ". I have created R3 and en_R3 file in this directory: " + dyn)
        print()
        count += 1
        
        file_inp  = "en_" + R1[:-3]
        # Execute file R1
        os.system("qsub -q " + c + " " + file_inp)
        print("qsub -q " + c + " " + file_inp)
        print(str(count) + ". I have executed en_R1 file in this directory: " + dyn)
        print()
        count += 1
        dyns_removed.append(dyn)
        dyns_.remove(dyn)
        dir_trajs.remove(dir_traj)
        break
  else: # Si el programa ha sido ejecutado antes, miramos que ficheros han sido ejecutados y seguimos ejecutando
    dyns_executed = list()
    dyns_executed_correctly = list()
    dyns_to_execute = list()
    i = 0
    os.chdir(mypath)
    
    with open(last_program_output, "r") as f:      
      for line in f:
        t = line.split()
        if (len(t) == 10 and t[2] == "have" and t[3] == "executed" and t[5] == "file" and t[8] == "directory:"):
          if t[9] not in dyns_executed:
            dyns_executed.append(t[9])
            print(str(count) + ". I have executed en_R1 file in this directory: " + t[9])  

    # Has executed well?? If not, run it again    
    for dyn in dyns_executed:
      if dyn == "None":
        dyns_executed_correctly.append(dyn)
        continue
      os.chdir(dyn)
      try:
        with open("Info_" + R3[:-3], "r") as f:
          for line in f:
            t = line.split()
            if (len(t) > 0):
              if (t[0] == "R3" and t[1] == "executed" and t[2] == "correctly"):
                dyns_executed_correctly.append(dyn)
            
        print("Executed correctly = " + dyn)

      except:
        dyns_to_execute.append(dyn)
        print(dyn + " has not executed well!! Take a look at the output file.")
        
    print()    
    os.chdir(mypath)
    
    for dyn in dyns_executed_correctly:
      dyns_removed.append(dyn)
      dyns.remove(dyn)
    
    for c in cpus_available:
      if (len(dyns_removed)%4 == 0) and len(dyns_removed) != 0:
        ii += 1
      molecule = molecules[ii]
      top_name = tops_name[ii]
      Atom_Name = Atom_Names[ii]
      prefix_traj = prefix_trajs[ii]
      fragment = fragments[ii]
      prefix_filpdb = prefix_filreps[ii]
      name_atom = name_atomreps[ii]
      fil_prot = fils_prot[ii]
        
      dyns_ = dyns

      for dyn in dyns_:
          if dyn == "None":
            print(str(count) + ". I have executed en_R1 file in this directory: " + dyn)
            count += 1
            continue
          if which_script == "first":
            os.chdir(mypath)
            # Move the files R1, R2 and R3 to each directory
            os.system("cp " + R1 + " " + dyn)
            os.system("cp " + R2 + " " + dyn)
            os.system("cp " + R3 + " " + dyn)
            os.chdir(dyn)  
            
            # which dir_traj we are using??
            which_dyn = dyn[-5:]
            which_mol = prefix_traj[-3:]
            dir_traj = next((s for s in dir_trajs if (which_dyn in s and which_mol in s)), None)           
 
            ## Create en_R1 (R1_...)        
            file_inp  = "en_" + R1[:-3] 
            file_out  = "Info_" + R1[:-3]
            file_     = open (file_inp,"w")

            line = "#!/bin/bash -f " + "\n"
            file_.write(line)
            line = "# "+ "\n"
            file_.write(line)
            line = "export DIRI=\"" + dyn + "\"" + "\n"
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

            en_R2 = "en_" + R2[:-3]
            
            line = "./" + R1 + " -info_from_top " + str(info_from_top) + " -top_name " + str(top_name) + " -Atom_Name " + str(Atom_Name) + " -SuperposeProt " + str(SuperposeProt) + " -ref_pdb " + str(ref_pdb) + " -dir_traj " + str(dir_traj) + " -prefix_traj " + str(prefix_traj) + " -ini_read "+ str(ini_read) + " -end_read "+ str(end_read) + " -ini_traj " + str(ini_traj) + " -end_traj " + str(end_traj) + " -inc_read " + str(inc_read) + " -full_traj " + str(full_traj) + " -last_pdb " + str(last_pdb) + " -clean " + str(clean) + " -en_R2 " + str(en_R2) + " -cpu " + str(c) + " > " + file_out + "\n"
            file_.write(line)

            
            run_cpptraj = file_inp
            os.system('chmod u+x {} '.format(run_cpptraj) )
            
            print(str(count) + ". I have created R1 and en_R1 file in this directory: " + dyn)
            print()
            count += 1

            ## Create the executable file R2_* (en_R2_) 
            file_inp  = "en_" + R2[:-3]
            file_out  = "Info_" + R2[:-3]
            file_     = open (file_inp,"w")

            line = "#!/bin/bash -f " + "\n"
            file_.write(line)
            line = "# "+ "\n"
            file_.write(line)
            line = "export DIRI=\"" + dyn + "\"" + "\n"
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

            en_R3  = "en_" + R3[:-3]
            
            line = "./" + R2 + " -dir_traj " + str(dir_traj) + " -top_name  " + str(top_name) + " -pocket " + str(pocket) + " -xray " + str(xray) + " -num_res_prot " + str(num_res_prot) + " -fil_prot " + str(fil_prot) + " -dis_min " + str(dis_min) + " -prefix_filpdb " + str(prefix_filpdb) + " -name_atom " + str(name_atom) + " -num_plots " + str(num_plots) + " -dis_plots " + str(dis_plots[0]) + "," + str(dis_plots[1]) + " -snaps_byone " + str(snaps_byone) + " -num_ns_anal " + str(num_ns_anal) + " -percent_anal " + str(percent_anal) + " -en_R3 " + str(en_R3) + " -cpu " + str(c) + " -fragment " + str(fragment) + " > " + file_out + "\n"
            file_.write(line)


            run_cpptraj = file_inp
            os.system('chmod u+x {} '.format(run_cpptraj) )


            print(str(count) + ". I have created R2 and en_R2 file in this directory: " + dyn)
            print()
            count += 1
            
            # Create the executable file R3
            file_inp  = "en_" + R3[:-3]
            file_out  = "Info_" + R3[:-3]
            file_     = open (file_inp,"w")

            line = "#!/bin/bash -f " + "\n"
            file_.write(line)
            line = "# "+ "\n"
            file_.write(line)
            line = "export DIRI=\"" + dyn + "\"" + "\n"
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


            line = "./" + R3 + " > " + file_out + "\n"
            file_.write(line)


            run_cpptraj = file_inp
            os.system('chmod u+x {} '.format(run_cpptraj) )
            # Don't run the R3 file
            #os.system(run_cpptraj)
            
            print(str(count) + ". I have created R3 and en_R3 file in this directory: " + dyn)
            print()
            count += 1
            
            file_inp  = "en_" + R1[:-3]
            # Execute file R1
            os.system("qsub -q " + c + " " + file_inp)
            print("qsub -q " + c + " " + file_inp)
            print(str(count) + ". I have executed en_R1 file in this directory: " + dyn)
            print()
            count += 1
            dyns_removed.append(dyn)
            dyns_.remove(dyn)
            dir_trajs.remove(dir_traj)
            break
          elif which_script == "second":
            os.chdir(mypath)
            # Move the files R2 and R3 to each directory
            os.system("cp " + R2 + " " + dyn)
            os.system("cp " + R3 + " " + dyn)
            os.chdir(dyn)  
            
            # which dir_traj we are using??
            which_dyn = dyn[-5:]
            which_mol = prefix_traj[-3:]
            dir_traj = next((s for s in dir_trajs if (which_dyn in s and which_mol in s)), None)

 
            ## Create the executable file R2_* (en_R2_) 
            file_inp  = "en_" + R2[:-3]
            file_out  = "Info_" + R2[:-3]
            file_     = open (file_inp,"w")

            line = "#!/bin/bash -f " + "\n"
            file_.write(line)
            line = "# "+ "\n"
            file_.write(line)
            line = "export DIRI=\"" + dyn + "\"" + "\n"
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

            en_R3  = "en_" + R3[:-3]
            
            line = "./" + R2 + " -dir_traj " + str(dir_traj) +" -top_name " + str(top_name) + " -pocket " + str(pocket) + " -xray " + str(xray) + " -num_res_prot " + str(num_res_prot) + " -fil_prot " + str(fil_prot) + " -dis_min " + str(dis_min) + " -prefix_filpdb " + str(prefix_filpdb) + " -name_atom " + str(name_atom) + " -num_plots " + str(num_plots) + " -dis_plots " + str(dis_plots[0]) + "," + str(dis_plots[1]) + " -snaps_byone " + str(snaps_byone) + " -num_ns_anal " + str(num_ns_anal) + " -percent_anal " + str(percent_anal) + " -en_R3 " + str(en_R3) + " -cpu " + str(c) + " -fragment " + str(fragment) + " > " + file_out + "\n"
            file_.write(line)


            run_cpptraj = file_inp
            os.system('chmod u+x {} '.format(run_cpptraj) )


            print(str(count) + ". I have created R2 and en_R2 file in this directory: " + dyn)
            print()
            count += 1
            
            # Create the executable file R3
            file_inp  = "en_" + R3[:-3]
            file_out  = "Info_" + R3[:-3]
            file_     = open (file_inp,"w")

            line = "#!/bin/bash -f " + "\n"
            file_.write(line)
            line = "# "+ "\n"
            file_.write(line)
            line = "export DIRI=\"" + dyn + "\"" + "\n"
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


            line = "./" + R3 + " > " + file_out + "\n"
            file_.write(line)


            run_cpptraj = file_inp
            os.system('chmod u+x {} '.format(run_cpptraj) )
            # Don't run the R3 file
            #os.system(run_cpptraj)
            
            print(str(count) + ". I have created R3 and en_R3 file in this directory: " + dyn)
            print()
            count += 1
            
            file_inp  = "en_" + R1[:-3]
            # Execute file R1
            os.system("qsub -q " + c + " " + file_inp)
            print("qsub -q " + c + " " + file_inp)
            print(str(count) + ". I have executed en_R1 file in this directory: " + dyn)
            print()
            count += 1
            dyns_removed.append(dyn)
            dyns_.remove(dyn)
            dir_trajs.remove(dir_traj)
            break
          elif which_script == "third":
            os.chdir(mypath)
            # Move the files R3 to each directory
            os.system("cp " + R3 + " " + dyn)
            os.chdir(dyn)  
            
            # which dir_traj we are using??
            which_dyn = dyn[-5:]
            which_mol = prefix_traj[-3:]
            dir_traj = next((s for s in dir_trajs if (which_dyn in s and which_mol in s)), None)
            
            # Create the executable file R3
            file_inp  = "en_" + R3[:-3]
            file_out  = "Info_" + R3[:-3]
            file_     = open (file_inp,"w")

            line = "#!/bin/bash -f " + "\n"
            file_.write(line)
            line = "# "+ "\n"
            file_.write(line)
            line = "export DIRI=\"" + dyn + "\"" + "\n"
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


            line = "./" + R3 + " > " + file_out + "\n"
            file_.write(line)


            run_cpptraj = file_inp
            os.system('chmod u+x {} '.format(run_cpptraj) )
            # Don't run the R3 file
            #os.system(run_cpptraj)
            
            print(str(count) + ". I have created R3 and en_R3 file in this directory: " + dyn)
            print()
            count += 1
            
            file_inp  = "en_" + R1[:-3]
            # Execute file R1
            os.system("qsub -q " + c + " " + file_inp)
            print("qsub -q " + c + " " + file_inp)
            print(str(count) + ". I have executed en_R1 file in this directory: " + dyn)
            print()
            count += 1
            dyns_removed.append(dyn)
            dyns_.remove(dyn)
            dir_trajs.remove(dir_traj)
            break
        

    

  os.chdir(mypath)
  
  



 




  return 



def cmdlineparse():

  parser = ArgumentParser(description="command line arguments")
  
  parser.add_argument("-molecule"    , dest="molecule" , required=True, help="Molecule name we are using as a ligand in that run" )
  parser.add_argument("-R1"    , dest="R1" , required=True, help="Version of the R1_ file we want to execute" )
  parser.add_argument("-R2"    , dest="R2" , required=True, help="Version of the R2_ file we want to execute" )
  parser.add_argument("-R3"    , dest="R3" , required=True, help="Version of the R3_ file we want to execute" )
  parser.add_argument("-dyn_1"    , dest="dyn_1" , required=True, help="Directory where is located the data of dyn_1" )
  parser.add_argument("-dyn_2"    , dest="dyn_2" , required=True, help="Directory where is located the data of dyn_2" )
  parser.add_argument("-dyn_3"    , dest="dyn_3" , required=True, help="Directory where is located the data of dyn_3" )
  parser.add_argument("-dyn_4"    , dest="dyn_4" , required=True, help="Directory where is located the data of dyn_4" )
  parser.add_argument("-cpus"    , dest="cpus" , required=True, help="List of Cpus we have available" )
  parser.add_argument("-owner_execution"    , dest="owner_execution" , required=True, help=" Who is the owner of the execution " )
  parser.add_argument("-last_program_output"    , dest="last_program_output" , required=True, help=" If the program has been executed before, we need to know the last file output of the program or you can leave as None" )
  parser.add_argument("-which_script"    , dest="which_script" , required=True, help=" Which script is the first to be launched (first, second, third)" )
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
  parser.add_argument("-end_read" ,    dest="end_read"      , required=True, help=" IFI for reading trajectories ")
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

  main()



