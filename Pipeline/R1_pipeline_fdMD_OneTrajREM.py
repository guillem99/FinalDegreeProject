#!/usr/bin/env python3
import numpy as np
import os, sys, stat
from datetime import datetime, date, time
from math import sqrt
import  matplotlib
matplotlib.use('Agg')
import  matplotlib.pyplot as plt
from argparse import ArgumentParser

""" __Input Data__  """
"""
        
"""
###########MAIN
def main() :

  global num_res_prot,num_lig,clean

  print ( " Using AMBER : ",amber_dir )
  mypath = os.getcwd()
  DIR      = os.getcwd()
  list_DIR = os.listdir(DIR)
  print ( " Initial Dir : ",DIR )
  print ( "  " )

  args = cmdlineparse()

  top_name      =      args.top_name  
  info_from_top =      args.info_from_top
  full_traj     =      args.full_traj
  last_pdb      =      args.last_pdb  
  clean         =      args.clean     

  if info_from_top == 'True' :
   info_from_top = True 
  else:
   info_from_top = False

  if full_traj == 'True':          
    full_traj = True
  else:
    full_traj = False

  if last_pdb  == 'True':          
    last_pdb  = True  
  else:
    last_pdb  = False

  if clean == 'True'   :          
    clean     = True  
  else:
    clean     = False 

  print ( " Using a system with : " )
  print ( " ... Topoplogy : {} ".format(top_name) )

  dir_traj    =      args.dir_traj
  top_file    =      dir_traj + "/" + top_name 

  if info_from_top :
    num_res_prot,num_lig,name_res_lig,pos_LIG,pos_lig_NoWat = GET_RES_from_top (top_file)  
  else :
    num_res_prot  = int ( args.num_res_prot )
    num_lig       = int ( args.num_lig      )
    name_res_lig  =       args.name_res_lig

  Atom_Name     =       args.Atom_Name

  SuperposeProt =       args.SuperposeProt

  Superpose = False
  if SuperposeProt   == "first" or SuperposeProt   == "reference" :
    Superpose = True
  elif SuperposeProt   == "xray" :
    Superpose = True

  if Superpose :
    ini_res_sup  = 1
    ifi_res_sup  = num_res_prot
    out_pdb_rmsd = False 
    if SuperposeProt == "reference" :
      ref_pdb        = args.ref_pdb
    elif SuperposeProt   == "xray" :
      ref_pdb        = args.ref_pdb
      lig_ref_pdb_sp = ref_pdb.split(".pdb")
      lig_ref_pdb    = lig_ref_pdb_sp[0] + "_LIG.pdb"
      out_pdb_rmsd   = True

      lig_ref_pdb_exist = lig_ref_pdb in list_DIR
      ref_pdb_exist     = ref_pdb     in list_DIR
      if (not lig_ref_pdb_exist) or ( not ref_pdb_exist ):
        print ( ' ... Using xray option but :' )
        if  not ref_pdb_exist :
          print ( " ... ... Reference PROTEIN is NOT present   : ",ref_pdb )
        if not lig_ref_pdb_exist :
          print ( " ... ... Reference LIGAND  is NOT present   : ",lig_ref_pdb )
        exit () 

    if args.ini_res_sup :
      ini_res_sup  = int ( ini_res_sup )
    if args.ifi_res_sup :
      ifi_res_sup  = int ( ifi_res_sup )

  prefix_traj =      args.prefix_traj
  ini_traj    = int ( args.ini_traj   )
  end_traj    = int ( args.end_traj   )
  num_traj    = end_traj - ini_traj + 1

  ini_read    = 1
  end_read    = "last"
  inc_read    = 1     
  if args.ini_read :
    ini_read  = int ( args.ini_read )
  if args.end_read :
    end_read  =     ( args.end_read )
    if end_read != "last" :
      end_read = int ( end_read ) 
      end_read = str ( end_read )
  if args.inc_read :
    inc_read  = int ( args.inc_read )

  using_full_traj   = True  
  num_snaps_traj    = -1
  if full_traj :
    num_snaps_traj  =  "last"
  else :
    if args.num_snaps_traj :
      num_snaps_traj  = int (args.num_snaps_traj)
      using_full_traj = False

  out_last_pdb    = False   
  num_snaps_total = -1 
  if last_pdb :
    out_last_pdb = True
    if args.num_snaps_total :
      num_snaps_total = int (args.num_snaps_total)

  if not info_from_top :
    print ( " ... Ligands : {:4d} Protein Atoms : {:6d} ".format(num_lig,num_res_prot) )
    print ( " ... Analysing the Atom/RES    : ",Atom_Name )

  if Superpose :
    print ( " ... Superposing Trajectory to : ",SuperposeProt )
    if SuperposeProt == "reference" :
      print ( " ... Using Protein as Reference from a COMPLEX   : ",ref_pdb )
    elif SuperposeProt == "xray" :
      print ( " ... Using Protein as Reference from a COMPLEX   : ",ref_pdb )
      print ( " ... Using Ligand  as Reference from a LIGAND    : ",lig_ref_pdb )
    print ( " ... superposing RES fom {:3d} to {:3d} ".format(ini_res_sup,ifi_res_sup) )

  print ( " ... Analysing Trajectories with Prefix   : ",prefix_traj )
  print ( " ... Trajectories files from {:3d} to {:3d} ".format(ini_traj,end_traj) )

  print ( " ... Using Snapshots from {:6d} to {} in {:2d} steps ".format(ini_read,end_read,inc_read ) )

  if out_last_pdb and num_snaps_total != -1 :
    print ( " ... Generating LAST pdb with the snapshot : {:6d} ".format(num_snaps_total) )

  if using_full_traj :
    print ( " ... Using the full lenght of the Total Trajectory " )
  else :
    print ( " ... Reading Only up to the : {:6s} snapshot".format(num_snaps_total) )

#
#...... RMSD to obtain a new .nc without WAT.....
#
  print ( " 1. cpptraj for trajectory : NoWat " )

  parm_wat      = "[Wat]"
  parm_NoWat    = "[NoWat]"
  parm_OneLigNW = "[OneLigNoWat]"
  parm_lig      = "[Lig]"

  file_cpptraj_inp  = "Ainptrj_nowat"
  file_cpptraj_out  = "Ainptrj_nowat.out"
  file_cpptraj      = open (file_cpptraj_inp,"w")

  line = " " + cpptraj + " >  " + file_cpptraj_out + " << EOF " + "\n"
  file_cpptraj.write(line)  
  line = " parm " + top_file + " " + parm_wat + "\n"
  file_cpptraj.write(line)  
  line = " " + "\n"
  file_cpptraj.write(line)  

  count_traj = 0
  for traj in range ( num_traj ) :
    indx_traj = ini_traj + count_traj
    name_traj = prefix_traj + "_" + str(indx_traj) + "_dyn.nc" 
    print ( "    ... Procesing Trajectory : ", name_traj )
    if amber_version > 16 :
      line = " trajin " + dir_traj + "/" + name_traj + " " + str(ini_read) + " " + end_read + " " + str (inc_read) + " parm " + parm_wat + "\n"
    else:
      line = " trajin " + dir_traj + "/" + name_traj + " " + str(ini_read) + " " + end_read + " " + str (inc_read) + " " + parm_wat + "\n"
    file_cpptraj.write(line)  
    count_traj += 1

  line = " " + "\n"
  file_cpptraj.write(line)  
  line = " autoimage " + "\n"
  file_cpptraj.write(line)  
  line = " strip " + amber_mask + " outprefix  NoWat " + "\n"
  top_NoWat = "NoWat." + top_name
  file_cpptraj.write(line)  
  line = " rms first out RMS_first.dat :" + str(ini_res_sup) + "-" + str(ifi_res_sup) + "@CA " + "\n"
  file_cpptraj.write(line)  

  traj_SupFirst_withoutWat = "RMSD_FIRST_NoWat.nc"
  line = " trajout " + traj_SupFirst_withoutWat + "\n"
  file_cpptraj.write(line)  

  line = " run " + "\n"
  file_cpptraj.write(line)  
  line = " " + "\n"
  file_cpptraj.write(line)  
  line = "EOF" + "\n"
  file_cpptraj.write(line)  

  file_cpptraj.close()  

  run_cpptraj = "./" + file_cpptraj_inp
  os.chmod  (run_cpptraj,0o0764)
  os.system (run_cpptraj) 

  if os.path.isfile(file_cpptraj_out) :
    file_cpptraj = open (file_cpptraj_out,"r")
    for line in file_cpptraj :
      if "frames and processed" in line :
        line_split      = line.split(" ")
        num_snaps_read = int ( line_split[5] ) 
        print ( " ... There are {:6d} Total Snapshots ".format(num_snaps_read) )
        if out_last_pdb :
          print ( " ... Generating LAST pdb with the snapshot : {:6d} ".format(num_snaps_read) )
        break
  else :
    print ( " --- cpptraj fails : STOP --- ) " ) 
    exit ()

  if clean :
    os.system  ('rm {} '.format(file_cpptraj_inp) )
    os.system  ('rm {} '.format(file_cpptraj_out) )

  if num_snaps_total < 0 :
    num_snaps_total  = num_snaps_read
  if num_snaps_traj != "last" :
    if num_snaps_traj  < 0 :
      num_snaps_traj   = num_snaps_read
  else:
    num_snaps_traj   = num_snaps_read

  """ 
    LIE : If the number of snapshots >  than max_snaps problems can appear
    In main ()  max_snaps     = 30000
  inc_read_LIE  = inc_read
  num_snaps_LIE = num_snaps_traj
  print ('1*',inc_read_LIE,num_snaps_LIE,max_snaps)
  if num_snaps_traj >= max_snaps :
    inc = 0 
    while num_snaps_LIE >= max_snaps and inc <= 100  :
      inc += 2 
      num_snaps_LIE = num_snaps_traj // 2
    if inc >= 100 :
      print ( ' ... There are Problems with LIE dimensions : STOP ! ' ) 
      exit ()
    inc_read_LIE = inc
    num_snaps_LIE = num_snaps_read // inc_read_LIE
    print ( " ... Using {:6d} Snapshots for LIE energy ".format(num_snaps_LIE) )
  """

  GET_lie_energy  (top_NoWat,traj_SupFirst_withoutWat) 

  print ( " 3. cpptraj for trajectory : Trajectory for Ligands alone " ) 
  
  count_ligand = num_res_prot + 1
  top_OneLigNoWat = "OneLig." + top_NoWat

  if num_lig != 1 :
    for ligand in range ( num_lig ) :
 
      if SuperposeProt   == "reference" or SuperposeProt   == "xray":
        name_nc_lig = "refer_lig_" + str(count_ligand) + ".nc"
      else :
        name_nc_lig =       "lig_" + str(count_ligand) + ".nc"

      file_cpptraj_inp  = "Binptrj_ligands"
      file_cpptraj_out  = "Binptrj_ligands.out"
      file_cpptraj      = open (file_cpptraj_inp,"w")

      line = " " + cpptraj + " >  " + file_cpptraj_out + " << EOF " + "\n"
      file_cpptraj.write(line)  
      line = " parm " + "NoWat." + top_name + " " + parm_NoWat + "\n"
      file_cpptraj.write(line)  
      line = " " + "\n"
      file_cpptraj.write(line)  
  
      if amber_version > 16 :
        line = " trajin " + traj_SupFirst_withoutWat + " " + "1" + " " + str(num_snaps_traj)  + " " + "1" + " parm " + parm_NoWat + "\n"
      else:
        line = " trajin " + traj_SupFirst_withoutWat + " " + "1" + " " + str(num_snaps_traj)  + " " + "1" + " " + parm_NoWat + "\n"
      file_cpptraj.write(line)  

      if ligand == 0 :
        line = " strip " + "!:1-" + str(num_res_prot) + "," + str(count_ligand) + " outprefix  OneLig " + "\n"
      else:
        line = " strip " + "!:1-" + str(num_res_prot) + "," + str(count_ligand) + "\n"
      file_cpptraj.write(line)  

      line = " trajout " + name_nc_lig + " " + parm_NoWat + "\n"
      file_cpptraj.write(line)  

      line = " run " + "\n"
      file_cpptraj.write(line)  
      line = " " + "\n"
      file_cpptraj.write(line)  
      line = "EOF" + "\n"
      file_cpptraj.write(line)  

      file_cpptraj.close()  

      run_cpptraj = "./" + file_cpptraj_inp
      os.chmod  (run_cpptraj,0o0764)
      os.system (run_cpptraj) 
      count_ligand += 1

    if clean :
      os.system  ('rm {} '.format(file_cpptraj_inp) )
      os.system  ('rm {} '.format(file_cpptraj_out) )

  else: 
    if SuperposeProt   == "reference" or SuperposeProt   == "xray":
      name_nc_lig = "refer_lig_" + str(count_ligand) + ".nc"
    else :
      name_nc_lig =       "lig_" + str(count_ligand) + ".nc"

    move_lig_nc    = " mv " + traj_SupFirst_withoutWat + " " + name_nc_lig
    copy_top_NoWat = " cp " + top_NoWat + " " + top_OneLigNoWat
    os.system (move_lig_nc   ) 
    os.system (copy_top_NoWat) 

  if SuperposeProt   == "reference" or SuperposeProt   == "xray":
    print ( " 4. cpptraj for trajectory : Protein Superposition for each Ligand alone  " ) 
  
    count_ligand = num_res_prot
    for ligand in range ( num_lig ) :
      count_ligand += 1
      name_nc_lig = "refer_lig_" + str(count_ligand) + ".nc"
      name_nc_out =       "lig_" + str(count_ligand) + ".nc"

      file_cpptraj_inp  = "Cinptrj_ligands_sup"
      file_cpptraj_out  = "Cinptrj_ligands_sup.out"
      file_cpptraj      = open (file_cpptraj_inp,"w")

      line = " " + cpptraj + " >  " + file_cpptraj_out + " << EOF " + "\n"
      file_cpptraj.write(line)  
      line = " parm " + "OneLig.NoWat." + top_name + " " + parm_OneLigNW + "\n"
      file_cpptraj.write(line)  
      line = " " + "\n"
      file_cpptraj.write(line)  
  
      if amber_version > 16 :
        line = " trajin " + name_nc_lig + " " + " parm " + parm_OneLigNW + "\n"
      else:
        line = " trajin " + name_nc_lig + " " + parm_OneLigNW + "\n"
      file_cpptraj.write(line)  
#     line = " reference " + "./" + ref_pdb + "\n"
      line = " reference " + ref_pdb + "\n"
      file_cpptraj.write(line)  
      line = " rms reference out RMS_reference.dat :" + str(ini_res_sup) + "-" + str(ifi_res_sup) + "@CA " + "\n"
      file_cpptraj.write(line)  
      line = " trajout " + name_nc_out + " " + parm_OneLigNW + "\n"
      file_cpptraj.write(line)  

      line = " run " + "\n"
      file_cpptraj.write(line)  
      line = " " + "\n"
      file_cpptraj.write(line)  
      line = "EOF" + "\n"
      file_cpptraj.write(line)  

      file_cpptraj.close()  

      run_cpptraj = "./" + file_cpptraj_inp
      os.system  ('chmod u+x {} '.format(run_cpptraj) )
      os.system (run_cpptraj) 
      
      os.remove ( name_nc_lig )

    if clean :
      os.system  ('rm {} '.format(file_cpptraj_inp) )
      os.system  ('rm {} '.format(file_cpptraj_out) )

  count_ligand = num_res_prot
  if Atom_Name != "NO" :
    print ( " 5. Generating files with the coordinates of the {} Atom".format(Atom_Name)  ) 
    for ligand in range ( num_lig ) :
      count_ligand += 1

      file_cpptraj_inp  = "Dinptrj_ligands_Atom"
      file_cpptraj_out  = "Dinptrj_ligands_Atom.out"
      file_cpptraj      = open (file_cpptraj_inp,"w")

      name_nc_lig     =  "lig_" + str(count_ligand) + ".nc"
      name_pdb_Atoms  =  "lig_" + str(count_ligand) + "_" + Atom_Name + ".pdb"

      line = " " + cpptraj + " >  " + file_cpptraj_out + " << EOF " + "\n"
      file_cpptraj.write(line)  
      line = " parm " + "OneLig.NoWat." + top_name + " " + parm_OneLigNW + "\n"
      file_cpptraj.write(line)  
      line = " " + "\n"
      file_cpptraj.write(line)  
 
      if amber_version > 16 :
        line = " trajin " + name_nc_lig + " " + " parm " + parm_OneLigNW + "\n"
      else:
        line = " trajin " + name_nc_lig + " " +  parm_OneLigNW + "\n"
      file_cpptraj.write(line)  

      line = " strip " + "!@" + Atom_Name + "\n"
      file_cpptraj.write(line)  
      line = " trajout " + name_pdb_Atoms + " pdb " + "\n"
      file_cpptraj.write(line)  

      line = " run " + "\n"
      file_cpptraj.write(line)  
      line = " " + "\n"
      file_cpptraj.write(line)  
      line = "EOF" + "\n"
      file_cpptraj.write(line)  

      file_cpptraj.close()  

      run_cpptraj = "./" + file_cpptraj_inp
      os.system  ('chmod u+x {} '.format(run_cpptraj) )
      os.system (run_cpptraj) 

      name_pdb_All    =  "lig_all_" + Atom_Name + ".pdb"

    if clean :
      os.system  ('rm {} '.format(file_cpptraj_inp) )
      os.system  ('rm {} '.format(file_cpptraj_out) )

  if out_last_pdb :
    print ( " 6. cpptraj for LAST pdb : For each Complex ( Protein + One_Ligand ) " )
  
    count_ligand = num_res_prot
    for ligand in range ( num_lig ) :
      count_ligand += 1
      name_nc_lig =       "lig_" + str(count_ligand) + ".nc"

      file_cpptraj_inp  = "Einptrj_ligands_last"
      file_cpptraj_out  = "Einptrj_ligands_last.out"
      file_cpptraj      = open (file_cpptraj_inp,"w")

      line = " " + cpptraj + " >  " + file_cpptraj_out + " << EOF " + "\n"
      file_cpptraj.write(line)  
      line = " parm " + "OneLig.NoWat." + top_name + " " + parm_OneLigNW + "\n"
      file_cpptraj.write(line)  
      line = " " + "\n"
      file_cpptraj.write(line)  
  
      if amber_version > 16 :
        line = " trajin " + name_nc_lig + " " + str(num_snaps_traj) + " " + str(num_snaps_traj) + " 1 " + " parm " + parm_OneLigNW + "\n"
      else:
        line = " trajin " + name_nc_lig + " " + str(num_snaps_traj) + " " + str(num_snaps_traj) + " 1 " + parm_OneLigNW + "\n"
      file_cpptraj.write(line)  
 
      name_pdb_last = "lig_" + str(count_ligand) + "_LAST.pdb" 
      line = " trajout " + name_pdb_last + " pdb " +  "\n"
      file_cpptraj.write(line)  

      line = " run " + "\n"
      file_cpptraj.write(line)  
      line = " " + "\n"
      file_cpptraj.write(line)  
      line = "EOF" + "\n"
      file_cpptraj.write(line)  

      file_cpptraj.close()  

      run_cpptraj = "./" + file_cpptraj_inp
      os.system  ('chmod u+x {} '.format(run_cpptraj) )
      os.system (run_cpptraj) 

    if clean :
      os.system  ('rm {} '.format(file_cpptraj_inp) )
      os.system  ('rm {} '.format(file_cpptraj_out) )

    print ( " 7. cpptraj SASA for LAST pdb : For each Complex ( Protein + One_Ligand ) " )
  
    count_ligand = num_res_prot
    num_res_tot  = num_res_prot + 1
    for ligand in range ( num_lig ) :
      count_ligand += 1
      name_pdb_last =     "lig_" + str(count_ligand) + "_LAST.pdb" 

      file_cpptraj_inp  = "Finptrj_ligands_sasa"
      file_cpptraj_out  = "Finptrj_ligands_sasa.out"
      file_cpptraj      = open (file_cpptraj_inp,"w")

      line = " " + cpptraj + " >  " + file_cpptraj_out + " << EOF " + "\n"
      file_cpptraj.write(line)  
      line = " parm " + "OneLig.NoWat." + top_name + " " + parm_OneLigNW + "\n"
      file_cpptraj.write(line)  
      line = " " + "\n"
      file_cpptraj.write(line)  
  
      if amber_version > 16 :
        line = " trajin " + name_pdb_last + " parm " + parm_OneLigNW + "\n"
      else:
        line = " trajin " + name_pdb_last            + parm_OneLigNW + "\n"
      file_cpptraj.write(line)  
 
      name_out_sasa = "lig_" + str(count_ligand) + "_SASA.out" 
      line = " molsurf  out  " + name_out_sasa + ' :1-' + str(num_res_tot) + "\n"
      file_cpptraj.write(line)  
      line = " run " + "\n"
      line = " molsurf  out  " + name_out_sasa + ' :1-' + str(num_res_prot) + "\n"
      file_cpptraj.write(line)  
      line = " run " + "\n"
      line = " molsurf  out  " + name_out_sasa + ' :' + str(num_res_tot) + "\n"
      file_cpptraj.write(line)  
      line = " run " + "\n"

      file_cpptraj.write(line)  
      line = " " + "\n"
      file_cpptraj.write(line)  
      line = "EOF" + "\n"
      file_cpptraj.write(line)  

      file_cpptraj.close()  

      run_cpptraj = "./" + file_cpptraj_inp
      os.system  ('chmod u+x {} '.format(run_cpptraj) )
      os.system (run_cpptraj) 

    if clean :
      os.system  ('rm {} '.format(file_cpptraj_inp) )
      os.system  ('rm {} '.format(file_cpptraj_out) )

    GET_dis_LIG_DCM (top_name,ini_read,end_read) 

    print ( " 9. cpptraj for ligand LAST Snapshot : one file with ALL LAST Ligands {}  ".format(name_res_lig) ) 

    file_cpptraj_inp  = "Ginptrj_all_ligands"
    file_cpptraj_out  = "Ginptrj_all_ligands.out"
    file_cpptraj      = open (file_cpptraj_inp,"w")

    line = " " + cpptraj + " >  " + file_cpptraj_out + " << EOF " + "\n"
    file_cpptraj.write(line)  
    line = " parm " + "OneLig.NoWat." + top_name + " " + parm_OneLigNW + "\n"
    file_cpptraj.write(line)  
    line = " " + "\n"
    file_cpptraj.write(line)  

    count_ligand = num_res_prot
    for ligand in range ( num_lig ) :
      count_ligand += 1
      name_nc_lig = "lig_" + str(count_ligand) + ".nc"
  
      if amber_version > 16 :
        line = " trajin " + name_nc_lig + " " + str(num_snaps_traj) + " " + str(num_snaps_traj) + " 1 " + " parm " + parm_OneLigNW + "\n"
      else:
        line = " trajin " + name_nc_lig + " " + str(num_snaps_traj) + " " + str(num_snaps_traj) + " 1 " + parm_OneLigNW + "\n"
      file_cpptraj.write(line)  

    line = " strip " + ":1-" + str(num_res_prot) + " outprefix  LigAlone " + "\n"
    top_LigAlone  = "LigAlone." + top_OneLigNoWat

    file_cpptraj.write(line)  
 
    name_pdb_last = "lig_All_" + name_res_lig + "_LAST.pdb" 
    line = " trajout " + name_pdb_last + " pdb " +  "\n"
    file_cpptraj.write(line)  

    line = " run " + "\n"
    file_cpptraj.write(line)  
    line = " " + "\n"
    file_cpptraj.write(line)  
    line = "EOF" + "\n"
    file_cpptraj.write(line)  

    file_cpptraj.close()  

    run_cpptraj = "./" + file_cpptraj_inp
    os.system  ('chmod u+x {} '.format(run_cpptraj) )
    os.system (run_cpptraj) 

    if clean :
      os.system  ('rm {} '.format(file_cpptraj_inp) )
      os.system  ('rm {} '.format(file_cpptraj_out) )

    print ( "10. cpptraj pdb for ALL Ligands {} alone the full trajectory ".format(name_res_lig) ) 

    file_cpptraj_inp  = "Hinptrj_full_ligands"
    file_cpptraj_out  = "Hinptrj_full_ligands.out"
    file_cpptraj      = open (file_cpptraj_inp,"w")

    line = " " + cpptraj + " >  " + file_cpptraj_out + " << EOF " + "\n"
    file_cpptraj.write(line)  
    line = " parm " + "OneLig.NoWat." + top_name + " " + parm_OneLigNW + "\n"
    file_cpptraj.write(line)  
    line = " " + "\n"
    file_cpptraj.write(line)  

    count_ligand = num_res_prot
    for ligand in range ( num_lig ) :
      count_ligand += 1
      name_nc_lig = "lig_" + str(count_ligand) + ".nc"
  
      if amber_version > 16 :
        line = " trajin " + name_nc_lig + " parm " + parm_OneLigNW + "\n"
      else:
        line = " trajin " + name_nc_lig + parm_OneLigNW + "\n"
      file_cpptraj.write(line)  

    line = " strip " + ":1-" + str(num_res_prot) + "\n"
    file_cpptraj.write(line)  
 
    name_pdb_last = "lig_All_" + name_res_lig + "_all.pdb" 
    line = " trajout " + name_pdb_last + " pdb " +  "\n"
    file_cpptraj.write(line)  

    line = " run " + "\n"
    file_cpptraj.write(line)  
    line = " " + "\n"
    file_cpptraj.write(line)  
    line = "EOF" + "\n"
    file_cpptraj.write(line)  

    file_cpptraj.close()  

    run_cpptraj = "./" + file_cpptraj_inp
    os.system  ('chmod u+x {} '.format(run_cpptraj) )
    os.system (run_cpptraj) 

    if clean :
      os.system  ('rm {} '.format(file_cpptraj_inp) )
      os.system  ('rm {} '.format(file_cpptraj_out) )

  if out_pdb_rmsd :
    print ( "11. RMSD_Reference for Ligand {} alone the full trajectory ".format(name_res_lig) ) 

    file_cpptraj_inp  = "Iinptrj_full_ligands_RMSD"
    file_cpptraj_out  = "Iinptrj_full_ligands_RMSD.out"
    file_cpptraj      = open (file_cpptraj_inp,"w")

    line = " " + cpptraj + " >  " + file_cpptraj_out + " << EOF " + "\n"
    file_cpptraj.write(line)  
    line = " parm " + top_LigAlone + " " + parm_lig + "\n"
    file_cpptraj.write(line)  
    line = " " + "\n"
    file_cpptraj.write(line)  

    name_pdb_last = "lig_All_" + name_res_lig + "_all.pdb" 
    if amber_version > 16 :
      line = " trajin " + name_pdb_last + " " + " parm " + parm_lig + "\n"
    else:
      line = " trajin " + name_pdb_last + " " + parm_lig + "\n"
    file_cpptraj.write(line)  

#   line = " reference " + "./" + lig_ref_pdb + "\n"
    line = " reference " + lig_ref_pdb + "\n"
    file_cpptraj.write(line)  

    line = " rms reference nofit out RMS_Lig_NoFit.dat " + ":1&!@H* " + "\n"
    file_cpptraj.write(line)  

    line = " run " + "\n"
    file_cpptraj.write(line)  
    line = " " + "\n"
    file_cpptraj.write(line)  
    line = "EOF" + "\n"
    file_cpptraj.write(line)  

    file_cpptraj.close()  

    run_cpptraj = "./" + file_cpptraj_inp
    os.chmod  (run_cpptraj,0o0764)
    os.system (run_cpptraj) 

    if clean :
      os.system  ('rm {} '.format(file_cpptraj_inp) )
      os.system  ('rm {} '.format(file_cpptraj_out) )

  # Change directory to the initial one.
  os.chdir(mypath)

  # Read final arguments
  en_R2 = args.en_R2 
  queue_system = args.queue_system
  os.chdir(mypath)
  
  # If the queue system is equal to qsub
  if (queue_system == "qsub"):
    # Read the argument to know which cpu we are going to use.
    cpu = args.cpu
    
    
    # Print what is going to be executed.
    print(queue_system + " -q " + cpu + " " + str(en_R2))

    # Execute next file R2.
    os.system(queue_system + " -q " + cpu + " " + str(en_R2))

 # If the queue system is not equal to qsub (it will be sbatch)
  else:
    # Print what is going to be executed.
    print(queue_system + " " + str(en_R2))

    # Execute next file R2
    os.system(queue_system + " " +  str(en_R2))

  # Recover the top_name file.
  os.system("cp " + top_file + " " + mypath)

  return  

def GET_lie_energy  (top_NoWat,traj_SupFirst_withoutWat) : 
  print ( " 2. LIE Interaction Energy ( vdW + Elec ) " )

  cc = '\n'
  parm_NoWat    = "[NoWat]"  

  file_cpptraj_inp  = "cpptraj_LIE_ener"
  file_cpptraj_out  = "cpptraj_LIE_ener.out"
  file_cpptraj      = open (file_cpptraj_inp,"w")

  file_cpptraj.write( "  {}  >  {}  << EOF ".format(cpptraj,file_cpptraj_out) + cc )
  file_cpptraj.write( " parm  {}  {} ".format(top_NoWat,parm_NoWat) + cc )
  file_cpptraj.write( cc )

  print ( "    ... Procesing Trajectory : ", traj_SupFirst_withoutWat )
  file_cpptraj.write(" trajin  {}  1 last 1  parm {}".format( traj_SupFirst_withoutWat,parm_NoWat ) + cc )

  file_cpptraj.write(' createcrd CRD1 ' + cc )
  file_cpptraj.write(' run ' + cc )            
  file_cpptraj.write( cc )            

  count_ligand  = num_res_prot
  num_res_tot   = num_res_prot + 1
  for ligand in range ( num_lig ) :
    count_ligand += 1
    str_count_ligand = str(count_ligand)
    name_nc_lig  = "lig_" + str_count_ligand  + ".nc"
    name_CAL_out = "lig_" + str_count_ligand  + "_LIE.dat"
    file_cpptraj.write( " crdaction CRD1 lie  :{}  :1-{}  out  {}".format(str_count_ligand,num_res_prot,name_CAL_out) + cc )

  file_cpptraj.write(' run ' + cc )            
  file_cpptraj.write( cc )            
  file_cpptraj.write( 'EOF' + cc )
  file_cpptraj.close()  

  run_cpptraj = "./" + file_cpptraj_inp
  os.system  ('chmod u+x {} '.format(run_cpptraj) )
  os.system (run_cpptraj) 

  if clean :
    os.system  ('rm {} '.format(run_cpptraj) )

  return 

def GET_dis_LIG_DCM (top_name,ini_read,end_read) :
  print ( " 8. Distance LIG - PROTEIN Center of Masses  " )
  cc = '\n'
  parm_OneLigNW     = "[OneLigNoWat]"
  file_cpptraj_inp  = "cpptraj_Dis_DCM"
  file_cpptraj_out  = "cpptraj_Dis_DCM.out"
  file_cpptraj      = open (file_cpptraj_inp,"w")

  file_cpptraj.write( "  {}  >  {}  << EOF ".format(cpptraj,file_cpptraj_out) + cc )
  file_cpptraj.write( " parm  OneLig.NoWat.{}  {} ".format(top_name,parm_OneLigNW) + cc )
  file_cpptraj.write( cc )

  count_ligand  = num_res_prot
  num_res_tot   = num_res_prot + 1
  for ligand in range ( num_lig ) :
    count_ligand += 1
    name_nc_lig  = "lig_" + str(count_ligand) + ".nc"
    name_CAL_out = "lig_" + str(count_ligand) + "_DCM.dat"
    file_cpptraj.write( " trajin    {}  {}  {}  1  parm  {} ".format(name_nc_lig,ini_read,end_read,parm_OneLigNW) + cc )
    file_cpptraj.write( " distance  :{}  :1-{}  out  {} ".format(num_res_tot,num_res_prot,name_CAL_out) + cc )
    file_cpptraj.write( ' run ' + cc )
    file_cpptraj.write( ' clear trajin ' + cc )
    file_cpptraj.write( cc )
  file_cpptraj.write( 'EOF' + cc )
  file_cpptraj.close()  

  run_cpptraj = "./" + file_cpptraj_inp
  os.system  ('chmod u+x {} '.format(run_cpptraj) )
  os.system (run_cpptraj) 

  return 

def GET_RES_from_top (top_name) :

  protein = ['SER','THR','GLN','ASN','TYR','CYS','CYX','CYM','GLY',  \
             'ALA','VAL','LEU','ILE','MET','PRO','PHE','TRP',        \
             'GLU','GLH','ASP','ASH','LYS','ARG','HIE','HID','HIP',  \
             'PHE','TYR','TRP' ]
  solvent = ['WAT' ]
  ions    = ['Cl-','Na+' ]

  space = ' '
  pos_LIG_line = ''

  sys_top       = open (top_name,"r")

  for line in sys_top :
    if "%FLAG POINTERS" in line :
      line_read     = sys_top.readline ().replace("\n","")
      line_read1    = sys_top.readline ().replace("\n","")
      line_read2    = sys_top.readline ().replace("\n","")
      line_read_sp1 = line_read1.split()
      line_read_sp2 = line_read2.split()
      num_atoms_SYS = int (line_read_sp1 [0])
      num_res_SYS   = int (line_read_sp2 [1])
      num_lines     = num_res_SYS  / 20
      num_lines_int = num_res_SYS // 20
      if (num_lines - num_lines_int) != 0 :
        num_lines_int += 1
      break    

  print ( ' ... Number of TOTAL atoms      : {:6d} '.format(num_atoms_SYS ) )
  print ( ' ... Number of TOTAL residues   : {:6d} '.format(num_res_SYS   ) )

  sys_top       = open (top_name,"r")

  residues = []
  pos_LIG  = []
  num_res_PROT = 0
  num_res_LIG  = 0
  num_res_SOL  = 0
  num_res_ION  = 0
  for line in sys_top :
    if "FLAG RESIDUE_LABEL" in line :
      line_read = sys_top.readline ().replace("\n","")
      for values in range (num_lines_int) :
        line_read    = sys_top.readline ().replace("\n","")
        line_read_sp = line_read.split()
        num_val      = len (line_read_sp)
        for num in range (num_val) :
          residues.append(line_read_sp[num])
      break

  indx = 0
  for res in residues :
    indx += 1
    if res in protein :
      num_res_PROT += 1
    elif res in solvent :
      num_res_SOL += 1
    elif res in ions :
      num_res_ION += 1
    else:
      num_res_LIG += 1
      name_res_LIG = res
      pos_LIG.append(indx)
      pos_LIG_char = str(indx)
      pos_LIG_line += pos_LIG_char + space

  pos_lig_NoWat  = []
  pos_LIG_LNoWat = ''
  for lig in range (num_res_LIG) :
    indx            = num_res_PROT + lig + 1
    pos_LIG_char    = str(indx)
    pos_LIG_LNoWat += pos_LIG_char + space
    pos_lig_NoWat.append(indx)
  pos_LIG_LNoWat +=  '\n'

  print ( ' ... Number of PROTEIN residues : {:6d} '.format(num_res_PROT ) )
  print ( ' ... Number of LIGAND  residues : {:6d} '.format(num_res_LIG  ) )
  print ( ' ... Number of IONS    residues : {:6d} '.format(num_res_ION  ) )
  print ( ' ... Number of SOLVENT residues : {:6d} '.format(num_res_SOL  ) )
  print ( ' ... LIGAND residue name        : {:>6} '.format(name_res_LIG ) )
  print ( ' ... LIGAND Positions           : ' )
  print ( ' ...        {:>6} '.format(pos_LIG_line   ) )
  print ( ' ...        {:>6} '.format(pos_LIG_LNoWat ) )

  return num_res_PROT,num_res_LIG,name_res_LIG,pos_LIG,pos_lig_NoWat 

def cmdlineparse():

    parser = ArgumentParser(description="command line arguments")

    parser.add_argument("-top_name"     , dest="top_name"     , required=True, help=" Name of the Topological file ")
    parser.add_argument("-info_from_top", dest="info_from_top", required=True, help=" TRUE : data from TOP    file ")
    parser.add_argument("-num_res_prot" , dest="num_res_prot" , required=False, help=" Number of Protein Residues ")
    parser.add_argument("-num_lig"      , dest="num_lig"      , required=False, help=" Number of Ligands ")
    parser.add_argument("-name_res_lig" , dest="name_res_lig" , required=False, help=" Residue Name  ")
    parser.add_argument("-Atom_Name"    , dest="Atom_Name"    , required=True,  help="  Atom to analyse (C99) ")
    parser.add_argument("-SuperposeProt", dest="SuperposeProt", required=True,  help="first or reference or xray  ")
    parser.add_argument("-ref_pdb" ,      dest="ref_pdb"      , required=True,  help=" file with the reference = Protein+Ligand ")
    parser.add_argument("-ini_res_sup" , dest="ini_res_sup"   , required=False, help=" INI resudue to superpose ")
    parser.add_argument("-ifi_res_sup" , dest="ifi_res_sup"   , required=False, help=" IFI resudue to superpose ") 
    parser.add_argument("-dir_traj" ,    dest="dir_traj"      , required=True, help=" DIR witn the trajectories  ")
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
    parser.add_argument("-en_R2"    ,    dest="en_R2"          , required=False, help=" Name of the next file to execute en_R2 ")
    parser.add_argument("-queue_system"    ,    dest="queue_system"          , required=False, help=" Command to run the script in bash ")
    parser.add_argument("-cpu"    ,    dest="cpu"          , required=False, help = " CPu we want to use to execute the files ")


    args=parser.parse_args()
    return args

def prog_info () :
  t = datetime.now()
  format_date = "%d-%m-%Y %H:%M:%S"
  ti= t.strftime(format_date)
  print ( "... ... ... ... ... ... " )
  print ( "...   fdMD_OneTraj  ... " )
  print ( "...     ( 2021 )    ... " )
  print ( "...     ........    ... " )
  print ( "...     . REM  .    ... " )
  print ( "... ... ... ... ... ... " )
  print ( " ",ti )
  print ( "... ... ... ... ... ... " )
  print ( " " )

def prog_END  () :
  t = datetime.now()
  format_date = "%d-%m-%Y %H:%M:%S"
  tf= t.strftime(format_date)
  print ( "... ... ... ... ... ... " )
  print ( "...   fdMD_OneTraj  ... " )
  print ( "...     ( 2021 )    ... " )
  print ( "...     ........    ... " )
  print ( "...     . REM  .    ... " )
  print ( "... ... ... ... ... ... " )
  print ( " ",tf )
  print ( "... ... ... ... ... ... " )
  print ( " " )

if __name__ == '__main__':

  x_prot   = []
  y_prot   = []
  z_prot   = []
  at_name  = []
  res_name = []
  res_num  = []
  x_lig    = []
  y_lig    = []
  z_lig    = []
  atom_name_lig = []
  resi_name_lig = ' '

  amber_version = 20
  amber_dir     =  "/programas/amber20/amber20"
  cpptraj       =  amber_dir + "/bin/cpptraj"

  amber_mask    = ":WAT,Na+,Cl-"

  max_snaps     = 30000

  prog_info ()
  main()
  prog_END  ()
