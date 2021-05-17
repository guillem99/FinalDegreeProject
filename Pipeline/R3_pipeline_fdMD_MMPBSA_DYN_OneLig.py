#!/usr/bin/python3

import os, sys
from datetime import datetime, date, time

#########################################
def main() :

  global AMBERHOME,AMBER_ENV,UTILS,TLEAP
  global ter,clean

  prog_info ()

  UTILS = '/home/prog/utils'
  TLEAP = 'tleap18'

  clean = False

  job_name  = 'GC_EBS'
  prefix    = 'lig'    

  do_PB          = False
  do_GB          = True
  force_do_top   = False
  rm_leap        = False 

  info_from_top  = True
  if not info_from_top :
    res_prot     = 000 # 521 Residuos de la proteina

  info_from_nc   = True

  dir_param      = '/home/guillemc/frag_params'
  lig_parm       = 'pd1.parm'
  lig_prep       = 'pd1.prep'

  system_name    = 'proteaseD_pd1'

  startf         = 1   
  nfreq          = 1
  endf           = 50000   

  iter_mmpbsa    = 0

  num_proc       = 1    # used only for slurm ( Almeria )

  max_calc_batch = 1   

  DeltaG_file = 'DeltaG_BySnaps.csv'

  indi  = 1.0
  exdi  = 80.0
  dielc = exdi

  mmpbsa_DYR  = 'mmpbsa'
  fftop       = 'fftop'

  if iter_mmpbsa == 0 :
    ffdir   = 'ffdir'
  else :
    ffdir   = 'ffdir' + '_' + str(iter_mmpbsa)

  AMBERHOME = '/home/prog/amber18'
  AMBER_ENV = '/home/prog/amber18/amber.sh'

########################################
  # Execute a bash command to delete the non-reactive trajectories (file created in R2).
  os.system("bash delete_*")


  ter = '\n'

  if ( not do_PB ) and ( not do_GB ) :
    print ( ' ... Select PB or GB : stop ')

  MAIN_work = os.getcwd()
  MAIN_work_list = os.listdir(MAIN_work)
  print (' ... ... Running ... ... ' )
  print (' *** Working in DYR : {} '.format(MAIN_work) )

  top_exist = False    
  for name_files in MAIN_work_list :
    if 'OneLig.NoWat' in name_files :
      top_file = name_files
      name_files_sp  = name_files.split('OneLig.NoWat.')
      top_name_Wat   = name_files_sp [1]
      top_exist = True
      break 
  if top_exist :
    print (' ... Topological file       = {} '.format(top_file) )
  else :
    print (' ... Topological file does not exist ')
    exit ()

  name_lig_files = []
  pose_name_list = []
  lig_file_num   = []
  for lig_files in MAIN_work_list :
    if 'lig_' in lig_files and '.nc' in lig_files :
      name_lig_files.append(lig_files)
  name_lig_files.sort()
  num_lig_files  = len(name_lig_files)

  for lig_files in name_lig_files :
    lig_files_sp      = lig_files.split('.nc')
    pose_neme_list_sp = lig_files_sp [0]
    lig_file_num_sp   = pose_neme_list_sp.split('_')  
    lig_file_num.append   (lig_file_num_sp [1]) 
    pose_name_list.append (pose_neme_list_sp  ) 
  print (' ... Number of LIG files = {} '.format(num_lig_files) )
  
  com_top = 'complex_com.top'
  rec_top = 'complex_rec.top'
  lig_top = 'complex_lig.top'

  dir_work = os.getcwd()
  dir_work_list = os.listdir(dir_work)

  name_queue_file_iqtc04 = 'en_mmpbsaDYN_iqtc04_all.sh'
  name_queue_file_iqtc06 = 'en_mmpbsaDYN_iqtc06_all.sh'
  name_queue_file_iqtc03 = 'en_mmpbsaDYN_iqtc03_all.sh'
  name_queue_file_iqtc09 = 'en_mmpbsaDYN_iqtc09_all.sh'
  name_queue_file_slurm  = 'en_mmpbsaDYN_slurm_all.sh'

  queue_file_iqtc04 = open(name_queue_file_iqtc04,'w')
  queue_file_iqtc06 = open(name_queue_file_iqtc06,'w')
  queue_file_iqtc03 = open(name_queue_file_iqtc03,'w')
  queue_file_iqtc09 = open(name_queue_file_iqtc09,'w')
  queue_file_slurm  = open(name_queue_file_slurm, 'w')

  tot_poses       = 0
  num_poses       = 0
  num_MMPBSA_Done = 0
  num_DYN_NotDone = 0

  do_leap  = True

  for ind_pose  in range ( num_lig_files ) :
    pose_name = pose_name_list [ind_pose]
    name_file = name_lig_files [ind_pose] 
    lig_num   = lig_file_num   [ind_pose]

    print (' ... Analysing  = {} -- {} '.format(pose_name,lig_num) )
    is_pose_dir = os.path.isdir(pose_name)
    if is_pose_dir == True :
      continue

    tot_poses += 1

    if tot_poses  == 1 :

      name_rst = system_name+'.rst'
      name_top = system_name+'.top'

      pdb_LAST = pose_name + '_LAST.pdb'
      pdb_LAST_exist =  pdb_LAST in MAIN_work_list
      if pdb_LAST_exist :
        print (' ... Using {} to Generate top without vdW_rep '.format(pdb_LAST) )
      else :
        print (' ... NO pdb to Generate top without vdW_rep ' )
        exit () 

      parm_file   = dir_param + '/' + lig_parm
      prep_file   = dir_param + '/' + lig_prep
      parm_NOempty = parm_file != [] 
      prep_NOempty = prep_file != [] 
      if ( not parm_NOempty ) or ( not  prep_NOempty ) : 
        print (' ... NO Parm or prep files ')
        print ( parm_file )
        print ( prep_file )
        exit ()

      leap_dir_exist = ('leap' in MAIN_work_list) and os.path.isdir('leap')
      if leap_dir_exist == False :
        os.mkdir('leap')
        os.chdir('leap')
        print (' ... leap DIR does not exist: Running LeAP ' )         
      else :
        if rm_leap :
          os.system ( 'rm -f -r leap ' )
          os.mkdir('leap')
          os.chdir('leap')
          print (' ... leap DIR removed : Running LeAP ' )         
        else :
          os.chdir('leap')
          LEAP_dir = os.getcwd()
          LEAP_dir_list = os.listdir(LEAP_dir)
          is_rst_file = name_rst in LEAP_dir_list
          is_top_file = name_top in LEAP_dir_list
          rst_NOempty   = is_rst_file and ( name_rst != [] )
          top_NOempty   = is_top_file and ( name_top != [] )
          if rst_NOempty and top_NOempty :
            do_leap = False
            print(' ... LeAP ... Done ')
          else: 
            print(' ... LeAP ... NOT Done ')

      if do_leap :
        os.system ( 'cp {} ./ '.format(parm_file) )
        os.system ( 'cp {} ./ '.format(prep_file) )
        os.system ( 'cp ../{} ./ '.format(pdb_LAST ) )
        Run_LeAP_NoWat (system_name,pdb_LAST,lig_parm,lig_prep)

      if info_from_top :
        res_prot,num_lig,name_res_lig = GET_RES_from_top (name_top)

      if info_from_nc  :
        endf = GET_snaps_from_nc (name_top,name_file)

      os.chdir(MAIN_work)

    os.mkdir(pose_name)
    os.chdir(pose_name)
    os.system ( 'cp ../{}      . '.format(name_file) )
    os.system ( 'cp ../leap/{} . '.format(name_top ) )

    MMPB_work = os.getcwd()
    MMPB_work_list = os.listdir(MMPB_work)
    is_dir_sel = os.path.isdir (mmpbsa_DYR)

    print (' ... DYR ... ',MMPB_work)

    do_fftop  = True
    do_mmpbsa = True

    if ( mmpbsa_DYR in MMPB_work_list ) and is_dir_sel :

      do_mmpbsa = False
      os.chdir(mmpbsa_DYR)
      mmpbsa_dir      = os.getcwd()
      mmpbsa_dir_list = os.listdir(mmpbsa_dir)
      is_fftop_dir = os.path.isdir(fftop)
      is_ffdir_dir = os.path.isdir(ffdir)

      if is_fftop_dir :
        os.chdir(fftop)
        fftop_dir      = os.getcwd()
        fftop_dir_list = os.listdir(fftop_dir)
        if com_top in fftop_dir_list and rec_top in fftop_dir_list and lig_top in fftop_dir_list :
          do_fftop = False
          os.chdir ('..')
          if force_do_top :
            do_fftop = True
            os.system ( 'rm -f -r {}'.format(fftop) )
        else :
          os.chdir ('..')
          os.system ( 'rm -f -r {}'.format(fftop) )

      if is_ffdir_dir  :
        os.chdir(ffdir)
        ffdir_dir      = os.getcwd()
        ffdir_dir_list = os.listdir(ffdir_dir)
        is_DeltaG_file = os.path.isfile(DeltaG_file)
        if is_DeltaG_file :
          print(' ... ... MMPBSA  of  {} done. OK.'.format(pose_name))
          num_MMPBSA_Done += 1
          os.chdir(dir_work)
          continue
        else :
          os.chdir ('..')
          os.system ( 'rm -f -r {}'.format(ffdir) )

      os.chdir ('..') # come back to pose

    num_poses += 1
    print(' ... ... {} : Running the MMPBSA for : {} '.format(num_poses,pose_name))

    if do_mmpbsa :
      os.mkdir(mmpbsa_DYR)
    os.chdir (mmpbsa_DYR)

    if do_fftop  :
      os.mkdir(fftop)
      os.chdir(fftop)
    
      RUN_anteMMPBSA ( name_top,res_prot )

    os.mkdir(ffdir)
    os.chdir(ffdir)
    ffdir_dir = os.getcwd()

    GEN_MMPBSA_inputs        ( do_PB,do_GB,startf,endf,nfreq,indi,exdi   )
    GEN_MMPBSA_inputs_iqtc04 ( system_name,job_name,lig_num,name_file )
    GEN_MMPBSA_inputs_iqtc06 ( system_name,job_name,lig_num,name_file )
    GEN_MMPBSA_inputs_iqtc09 ( system_name,job_name,lig_num,name_file )
    GEN_MMPBSA_inputs_iqtc03 ( system_name,job_name,lig_num,name_file )
    GEN_MMPBSA_inputs_slurm  ( system_name,job_name,lig_num,name_file,num_proc )

    os.chdir(dir_work)

    queue_file_iqtc04.write('cd ./{}/mmpbsa/{}'.format(pose_name,ffdir)+ter)
    queue_file_iqtc04.write('qsub -q iqtc04.q en_{}_mmpbsa_iqtc04'.format(lig_num)+ter)
    queue_file_iqtc04.write('cd ../../../'+ter)

    queue_file_iqtc06.write('cd ./{}/mmpbsa/{}'.format(pose_name,ffdir)+ter)
    queue_file_iqtc06.write('qsub -q iqtc06.q en_{}_mmpbsa_iqtc06'.format(lig_num)+ter)
    queue_file_iqtc06.write('cd ../../../'+ter)

    queue_file_iqtc09.write('cd ./{}/mmpbsa/{}'.format(pose_name,ffdir)+ter)
    queue_file_iqtc09.write('qsub -q iqtc09_g1112.q en_{}_mmpbsa_iqtc09'.format(lig_num)+ter)
    queue_file_iqtc09.write('cd ../../../'+ter)

    queue_file_iqtc03.write('cd ./{}/mmpbsa/{}'.format(pose_name,ffdir)+ter)
    queue_file_iqtc03.write('qsub -q g6.q en_{}_mmpbsa_iqtc03'.format(lig_num)+ter)
    queue_file_iqtc03.write('cd ../../../'+ter)

    queue_file_slurm.write ('cd ./{}/mmpbsa/{}'.format(pose_name,ffdir)+ter)
    queue_file_slurm.write ('sbatch en_{}_mmpbsa_slurm'.format(lig_num)+ter)
    queue_file_slurm.write ('cd ../../../'+ter)

  queue_file_iqtc04.close()
  queue_file_iqtc06.close()
  queue_file_iqtc09.close()
  queue_file_iqtc03.close()
  queue_file_slurm.close()

  os.system('chmod u+x en_mmpbsaDYN_iqtc04_all.sh')
  os.system('chmod u+x en_mmpbsaDYN_iqtc06_all.sh')
  os.system('chmod u+x en_mmpbsaDYN_iqtc09_all.sh')
  os.system('chmod u+x en_mmpbsaDYN_iqtc03_all.sh')
  os.system('chmod u+x en_mmpbsaDYN_slurm_all.sh' )

  print ( ' ')
  print ( ' Total number of Poses = {} '.format (tot_poses       ) )
  print ( ' Num DYN Not Done      = {} '.format (num_DYN_NotDone   ) )
  print ( ' Num MMPBSA not Done   = {} '.format (num_poses       ) )
  print ( ' Num MMPBSA     Done   = {} '.format (num_MMPBSA_Done ) )
  print ( ' ')
  print ( ' ! Remember to RUN de MMPBSA shell ! ' )
  print ( ' ')


  if max_calc_batch == 0 :
    exit ()

  queue_file_iqtc04 = open(name_queue_file_iqtc04,'r')
  queue_file_iqtc06 = open(name_queue_file_iqtc06,'r')
  queue_file_iqtc09 = open(name_queue_file_iqtc09,'r')
  queue_file_iqtc03 = open(name_queue_file_iqtc03,'r')
  queue_file_slurm  = open(name_queue_file_slurm, 'r')

  num_poses_batch = []
  num_batch_send = num_poses // max_calc_batch

  for batch in range ( num_batch_send ) :
    num_poses_batch.append ( max_calc_batch )

  res_poses_send = num_poses %  max_calc_batch
  if res_poses_send != 0 :
    num_batch_send += 1
    num_poses_batch.append ( res_poses_send )

  for batch in range ( num_batch_send ) :
    poses_to_send = num_poses_batch [batch]

# iqtc04
    name_send = 'en_mmmpbsa_iqtc04_all_' + str(batch) + '.sh'
    ANTE_send = open(name_send,'w')

    for pose in range (poses_to_send) :
      line = queue_file_iqtc04.readline ().replace("\n","")
      ANTE_send.write(line+ter)
      line = queue_file_iqtc04.readline ().replace("\n","")
      ANTE_send.write(line+ter)
      line = queue_file_iqtc04.readline ().replace("\n","")
      ANTE_send.write(line+ter)

    ANTE_send.close()
    os.system('chmod u+x {}'.format(name_send))

# iqtc06
    name_send = 'en_mmpbsa_iqtc06_all_' + str(batch) + '.sh'
    ANTE_send = open(name_send,'w')

    for pose in range (poses_to_send) :
      line = queue_file_iqtc06.readline ().replace("\n","")
      ANTE_send.write(line+ter)
      line = queue_file_iqtc06.readline ().replace("\n","")
      ANTE_send.write(line+ter)
      line = queue_file_iqtc06.readline ().replace("\n","")
      ANTE_send.write(line+ter)

    ANTE_send.close()
    os.system('chmod u+x {}'.format(name_send))

# iqtc09
    name_send = 'en_mmpbsa_iqtc09_all_' + str(batch) + '.sh'
    ANTE_send = open(name_send,'w')

    for pose in range (poses_to_send) :
      line = queue_file_iqtc09.readline ().replace("\n","")
      ANTE_send.write(line+ter)
      line = queue_file_iqtc09.readline ().replace("\n","")
      ANTE_send.write(line+ter)
      line = queue_file_iqtc09.readline ().replace("\n","")
      ANTE_send.write(line+ter)

    ANTE_send.close()
    os.system('chmod u+x {}'.format(name_send))

# iqtc03  g6.q 
    name_send = 'en_mmpbsa_iqtc03_all_' + str(batch) + '.sh'
    ANTE_send = open(name_send,'w')

    for pose in range (poses_to_send) :
      line = queue_file_iqtc03.readline ().replace("\n","")
      ANTE_send.write(line+ter)
      line = queue_file_iqtc03.readline ().replace("\n","")
      ANTE_send.write(line+ter)
      line = queue_file_iqtc03.readline ().replace("\n","")
      ANTE_send.write(line+ter)

    ANTE_send.close()
    os.system('chmod u+x {}'.format(name_send))

# slurm 
    name_send = 'en_mmpbsa_slurm_all_' + str(batch) + '.sh'
    ANTE_send = open(name_send,'w')

    for pose in range (poses_to_send) :
      line = queue_file_slurm.readline ().replace("\n","")
      ANTE_send.write(line+ter)
      line = queue_file_slurm.readline ().replace("\n","")
      ANTE_send.write(line+ter)
      line = queue_file_slurm.readline ().replace("\n","")
      ANTE_send.write(line+ter)

    ANTE_send.close()
    os.system('chmod u+x {}'.format(name_send))

  queue_file_iqtc04.close()
  queue_file_iqtc06.close()
  queue_file_iqtc09.close()
  queue_file_slurm.close()

  print("R3 executed correctly")


  return

def GET_snaps_from_nc (com_top,name_file) :

  parm_OneLigNW     = '[OneLigNoWat]'
  file_cpptraj_inp  = "cpptraj_snaps.inp"
  file_cpptraj_out  = "cpptraj_snaps.out"
  run_cpptraj       = open (file_cpptraj_inp,"w")

  run_cpptraj.write(' {} > {} << EOF'.format(cpptraj,file_cpptraj_out)+ter)
  run_cpptraj.write(' parm {} {} '.format(com_top,parm_OneLigNW)+ter)
  run_cpptraj.write(' trajin ../{} 1 last 1 parm {}'.format(name_file,parm_OneLigNW)+ter)
  run_cpptraj.write(' run'+ter)
  run_cpptraj.write('EOF' +ter)

  run_cpptraj.close()
  os.system('chmod u+x {}'.format(file_cpptraj_inp) )
  RUN_cpptraj = "./" + file_cpptraj_inp
  os.system(RUN_cpptraj)

  num_snaps = -1 
  if os.path.isfile(file_cpptraj_out) :
    file_cpptraj = open (file_cpptraj_out,"r")
    for line in file_cpptraj :
      if "frames and processed" in line :
        line_split      = line.split(" ")
        num_snaps = int ( line_split[5] )
        print ( " ... There are {:6d} Total Snapshots ".format(num_snaps) )
        break
  else :
    print ( " --- cpptraj fails to obtain num_snaps : STOP --- ) " )
    exit ()

  if clean :
    os.system  ('rm {} '.format(file_cpptraj_inp) )
    os.system  ('rm {} '.format(file_cpptraj_out) )

  return num_snaps

def Run_LeAP_NoWat (system_name,pdb,parm,prep) :

  leap_script = open('tleap_input.script','w')

  leap_script.write('source leaprc.protein.ff14SB       '+ter)
  leap_script.write('source leaprc.gaff2'+ter)
  leap_script.write(ter)
  leap_script.write('set default pbradii mbondi2'+ter)
  leap_script.write(ter)
  leap_script.write('com_par = loadamberparams {}'.format(parm)+ter)
  leap_script.write('loadamberprep {}'.format(prep)+ter)
  leap_script.write(ter)
  leap_script.write('complex = loadpdb {}'.format(pdb)+ter)
  leap_script.write(ter)
  leap_script.write('saveamberparm complex {}.top {}.rst'.format(system_name,system_name)+ter)
  leap_script.write(ter)
  leap_script.write('quit'+ter)
  leap_script.write(ter)
  leap_script.close()

  leap_run = open('run_tleap18','w')
  leap_run.write('#! /bin/csh'+ter)
  leap_run.write('#'+ter)
  leap_run.write('##########################################'+ter)
  leap_run.write(ter)
  leap_run.write(' {}/{} -f tleap_input.script'.format(UTILS,TLEAP)+ter)
  leap_run.write(ter)
  leap_run.close()

  os.system('chmod u+x run_tleap18')
  os.system('./run_tleap18')

  return

#================================================================
# _GET_RES_from_top
#================================================================
def GET_RES_from_top (top_name) :

  protein = ['SER','THR','GLN','ASN','TYR','CYS','CYX','CYM','GLY',  \
             'ALA','VAL','LEU','ILE','MET','PRO','PHE','TRP',        \
             'GLU','GLH','ASP','ASH','LYS','ARG','HIE','HID','HIP',  \
             'PHE','TYR','TRP' ]
  solvent = ['WAT' ]
  ions    = ['Cl-','Na+' ]

  sys_top      = open (top_name,"r")

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
  indx     = 0
  num_res_PROT = 0
  num_res_LIG  = 0
  num_res_SOL  = 0
  num_res_ION  = 0
  for line in sys_top :
    indx += 1
    if "FLAG RESIDUE_LABEL" in line :
      line_read = sys_top.readline ().replace("\n","")
      for values in range (num_lines_int) :
        line_read    = sys_top.readline ().replace("\n","")
        line_read_sp = line_read.split()
        num_val      = len (line_read_sp)
        for num in range (num_val) :
          residues.append(line_read_sp[num])
      break
  for res in residues :
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

  print ( ' ... Number of PROTEIN residues : {:6d} '.format(num_res_PROT ) )
  print ( ' ... Number of LIGAND  residues : {:6d} '.format(num_res_LIG  ) )
  print ( ' ... Number of IONS    residues : {:6d} '.format(num_res_ION  ) )
  print ( ' ... Number of SOLVENT residues : {:6d} '.format(num_res_SOL  ) )
  print ( ' ... LIGAND residue name        : {:>6} '.format(name_res_LIG ) )

  return num_res_PROT,num_res_LIG,name_res_LIG

#================================================================
# _RUN_anteMMPBSA_
#================================================================
def RUN_anteMMPBSA ( top_file,res_prot ) :
  antePB_file = open('run_anteMMPBSA','w')
  antePB_file.write('#!/bin/bash'+ter)
  antePB_file.write('AMBERHOME={}'.format(AMBERHOME)+ter)
  antePB_file.write('source {}'.format(AMBER_ENV)+ter)
  antePB_file.write(ter)
  antePB_file.write('$AMBERHOME/bin/ante-MMPBSA.py       \\'+ter)
  antePB_file.write('  -p ../../{}                       \\'.format(top_file)+ter)
  antePB_file.write('  -c ./complex_com.top              \\'+ter)
  antePB_file.write('  -r ./complex_rec.top              \\'+ter)
  antePB_file.write('  -l ./complex_lig.top              \\'+ter)
  #antePB_file.write("  -s ':WAT,Na+,Cl-'                 \\"+ter)
  antePB_file.write("  -m ':1-{}'                        \\".format(res_prot)+ter)
  antePB_file.write('  --radii=mbondi2 >> ante-MMPBSA.log'+ter)
  antePB_file.close()

  os.system ( 'chmod u+x run_anteMMPBSA' )
  os.system ( './run_anteMMPBSA' )
  os.system ( 'cp ../../{}  ./complex_com.top '.format(top_file ) )
  os.chdir('..')

  return

#================================================================
# _GEN_MMPBSA_inputs_
#================================================================
def GEN_MMPBSA_inputs ( do_PB,do_GB,startf,endf,nfreq,indi,exdi ) :
  MMPBSA_inp = open('mmpbsa_py.inp','w')
  MMPBSA_inp.write('&general'+ter)
  MMPBSA_inp.write(' startframe={}, endframe={}, interval={},'.format(startf,endf,nfreq) +ter)
  MMPBSA_inp.write(' keep_files=0,'+ter)
  MMPBSA_inp.write('/'+ter)
  if do_PB :
    MMPBSA_inp.write('&pb'+ter)
    MMPBSA_inp.write(' cavity_surften=0.00542, cavity_offset=0.9200,'+ter)
    MMPBSA_inp.write(' indi={}, exdi={},'.format(indi,exdi)+ter)
    MMPBSA_inp.write(' istrng=0.0, radiopt=0, inp=1,'+ter)
    MMPBSA_inp.write(' fillratio=4.0, scale=2.0,'+ter)
    MMPBSA_inp.write('/'+ter)
  if do_GB :
    MMPBSA_inp.write('&gb'+ter)
    MMPBSA_inp.write(' igb=2,'+ter)
    MMPBSA_inp.write(' surften=0.0072,'+ter)
    MMPBSA_inp.write(' saltcon=0.0,'+ter)
    MMPBSA_inp.write('/'+ter)
    MMPBSA_inp.close()
  return

def GEN_MMPBSA_inputs_iqtc04 ( system_name,job_name,num_poses,name_file ) :
# 
# iqtc04
#
  MMPBSA_send = open('en_{}_mmpbsa_iqtc04'.format(num_poses),'w')
  MMPBSA_send.write('#!/bin/bash'+ter)
  MMPBSA_send.write('######################################'+ter)
  MMPBSA_send.write('# Opcions i parametres del SGE        '+ter)
  MMPBSA_send.write('######################################'+ter)
  MMPBSA_send.write('# (1) Nom del treball                 '+ter)
  MMPBSA_send.write('#$ -N {}_{}                           '.format(job_name,num_poses)+ter)
  MMPBSA_send.write('# (2) Recursos sol.licitats           '+ter)
  MMPBSA_send.write('##$ -l h_rt                           '+ter)
  MMPBSA_send.write('##$ -l mem_free                       '+ter)
  MMPBSA_send.write('#$ -pe smp 12                         '+ter)
  MMPBSA_send.write('##$ -l exclusive=true                 '+ter)
  MMPBSA_send.write('# (3) Fitxers de sortida              '+ter)
  MMPBSA_send.write('#$ -cwd                               '+ter)
  MMPBSA_send.write('#$ -o {}.out                          '.format(job_name)+ter)
  MMPBSA_send.write('#$ -e {}.err                          '.format(job_name)+ter)
  MMPBSA_send.write('# (4) Envia un mail                   '+ter)
  MMPBSA_send.write('##$ -m e                              '+ter)
  MMPBSA_send.write('##$ -M jaime.rubio@ub.edu             '+ter)
  MMPBSA_send.write('######################################'+ter)
  MMPBSA_send.write('# Entorn de usuari                    '+ter)
  MMPBSA_send.write('######################################'+ter)
  MMPBSA_send.write('# Es carreguen els moduls             '+ter)
  MMPBSA_send.write('source /etc/profile                   '+ter)
  MMPBSA_send.write('env                                   '+ter)
  MMPBSA_send.write('module load numpy/1.6.1               '+ter)
  MMPBSA_send.write('module load amber/16_ompi             '+ter)
  MMPBSA_send.write('######################################'+ter)
  MMPBSA_send.write('echo " SGE_O_WORKDIR : $SGE_O_WORKDIR"'+ter)
  MMPBSA_send.write('echo " nslots        : $NSLOTS"       '+ter)
  MMPBSA_send.write('echo " TMP DIR       : $TMPDIR"       '+ter)
  MMPBSA_send.write('######################################'+ter)
  MMPBSA_send.write('# Calcul                              '+ter)
  MMPBSA_send.write('######################################'+ter)
  MMPBSA_send.write('export RUN=/home/g6jaime/jaime/run    '+ter)
  MMPBSA_send.write('echo " RUN DIR       : $RUN"          '+ter)
  MMPBSA_send.write('echo "AMBERHOME   : $AMBERHOME "      '+ter)
  MMPBSA_send.write('source  $AMBERHOME/amber.sh           '+ter)
  MMPBSA_send.write('#'+ter)
  MMPBSA_send.write('cd ../fftop                           '+ter)
  MMPBSA_send.write('export TOP_DIR=`pwd`                  '+ter)
  MMPBSA_send.write('cd $SGE_O_WORKDIR                     '+ter)
  MMPBSA_send.write('cd ../../                             '+ter)
  MMPBSA_send.write('export DYN_DIR=`pwd`                  '+ter)
  MMPBSA_send.write('echo "TOP DIR   : $TOP_DIR "          '+ter)
  MMPBSA_send.write('echo "DYN DIR   : $DYN_DIR "          '+ter)
  MMPBSA_send.write('#'+ter)
  MMPBSA_send.write('cd $SGE_O_WORKDIR                     '+ter)
  MMPBSA_send.write('cp $SGE_O_WORKDIR/mmpbsa_py.inp $TMPDIR   '+ter)
  MMPBSA_send.write('cd $TMPDIR                            '+ter)
  MMPBSA_send.write('#'+ter)
  MMPBSA_send.write('mpirun -np 12 $AMBERHOME/bin/MMPBSA.py.MPI -O  \\'+ter)
  MMPBSA_send.write('             -i  ./mmpbsa_py.inp              \\'+ter)
  MMPBSA_send.write('             -o  ./mmpbsa_py.out              \\'+ter)
  MMPBSA_send.write('             -eo ./DeltaG_BySnaps.csv         \\'+ter)
  MMPBSA_send.write('             -sp $DYN_DIR/{}.top      \\'.format(system_name)+ter)
  MMPBSA_send.write('             -cp $TOP_DIR/complex_com.top     \\'+ter)
  MMPBSA_send.write('             -rp $TOP_DIR/complex_rec.top     \\'+ter)
  MMPBSA_send.write('             -lp $TOP_DIR/complex_lig.top     \\'+ter)
  MMPBSA_send.write('             -y  $DYN_DIR/{}          \\'.format(name_file)+ter)
  MMPBSA_send.write('             >>  ./progres.log'+ter)
  MMPBSA_send.write('#'+ter)
  MMPBSA_send.write('rm -f -r _MMPBSA_*'+ter)
  MMPBSA_send.write('rm -f -r reference.frc'+ter)
  MMPBSA_send.write('#'+ter)
  MMPBSA_send.write('cp $TMPDIR/* $SGE_O_WORKDIR/'+ter)
  MMPBSA_send.write('#'+ter)
  MMPBSA_send.close()
  os.system('chmod u+x en_{}_mmpbsa_iqtc04'.format(num_poses))
  return

def GEN_MMPBSA_inputs_iqtc06 ( system_name,job_name,num_poses,name_file) :
# 
# iqtc06
#
  MMPBSA_send = open('en_{}_mmpbsa_iqtc06'.format(num_poses),'w')
  MMPBSA_send.write('#!/bin/bash'+ter)
  MMPBSA_send.write('######################################'+ter)
  MMPBSA_send.write('# Opcions i parametres del SGE        '+ter)
  MMPBSA_send.write('######################################'+ter)
  MMPBSA_send.write('# (1) Nom del treball                 '+ter)
  MMPBSA_send.write('#$ -N {}_{}                           '.format(job_name,num_poses)+ter)
  MMPBSA_send.write('# (2) Recursos sol.licitats           '+ter)
  MMPBSA_send.write('##$ -l h_rt                           '+ter)
  MMPBSA_send.write('##$ -l mem_free                       '+ter)
  MMPBSA_send.write('#$ -pe smp 16                         '+ter)
  MMPBSA_send.write('##$ -l exclusive=true                 '+ter)
  MMPBSA_send.write('# (3) Fitxers de sortida              '+ter)
  MMPBSA_send.write('#$ -cwd                               '+ter)
  MMPBSA_send.write('#$ -o {}.out                          '.format(job_name)+ter)
  MMPBSA_send.write('#$ -e {}.err                          '.format(job_name)+ter)
  MMPBSA_send.write('# (4) Envia un mail                   '+ter)
  MMPBSA_send.write('##$ -m e                              '+ter)
  MMPBSA_send.write('##$ -M jaime.rubio@ub.edu             '+ter)
  MMPBSA_send.write('######################################'+ter)
  MMPBSA_send.write('# Entorn de usuari                    '+ter)
  MMPBSA_send.write('######################################'+ter)
  MMPBSA_send.write('# Es carreguen els moduls             '+ter)
  MMPBSA_send.write('source /etc/profile                   '+ter)
  MMPBSA_send.write('env                                   '+ter)
  MMPBSA_send.write('module load ambertools/16_ompi        '+ter)
  MMPBSA_send.write('######################################'+ter)
  MMPBSA_send.write('echo " SGE_O_WORKDIR : $SGE_O_WORKDIR"'+ter)
  MMPBSA_send.write('echo " nslots        : $NSLOTS"       '+ter)
  MMPBSA_send.write('echo " TMP DIR       : $TMPDIR"       '+ter)
  MMPBSA_send.write('######################################'+ter)
  MMPBSA_send.write('# Calcul                              '+ter)
  MMPBSA_send.write('######################################'+ter)
  MMPBSA_send.write('export RUN=/home/g6jaime/jaime/run    '+ter)
  MMPBSA_send.write('echo " RUN DIR       : $RUN"          '+ter)
  MMPBSA_send.write('echo "AMBERHOME   : $AMBERHOME "      '+ter)
  MMPBSA_send.write('source  $AMBERHOME/amber.sh           '+ter)
  MMPBSA_send.write('#'+ter)
  MMPBSA_send.write('cd ../fftop                           '+ter)
  MMPBSA_send.write('export TOP_DIR=`pwd`                  '+ter)
  MMPBSA_send.write('cd $SGE_O_WORKDIR                     '+ter)
  MMPBSA_send.write('cd ../../                             '+ter)
  MMPBSA_send.write('export DYN_DIR=`pwd`                  '+ter)
  MMPBSA_send.write('echo "TOP DIR   : $TOP_DIR "          '+ter)
  MMPBSA_send.write('echo "DYN DIR   : $DYN_DIR "          '+ter)
  MMPBSA_send.write('#'+ter)
  MMPBSA_send.write('cd $SGE_O_WORKDIR                     '+ter)
  MMPBSA_send.write('cp $SGE_O_WORKDIR/mmpbsa_py.inp $TMPDIR   '+ter)
  MMPBSA_send.write('cd $TMPDIR                            '+ter)
  MMPBSA_send.write('#'+ter)
  MMPBSA_send.write('mpirun -np 16 $AMBERHOME/bin/MMPBSA.py.MPI -O  \\'+ter)
  MMPBSA_send.write('             -i  ./mmpbsa_py.inp              \\'+ter)
  MMPBSA_send.write('             -o  ./mmpbsa_py.out              \\'+ter)
  MMPBSA_send.write('             -eo ./DeltaG_BySnaps.csv         \\'+ter)
  MMPBSA_send.write('             -sp $DYN_DIR/{}.top      \\'.format(system_name)+ter)
  MMPBSA_send.write('             -cp $TOP_DIR/complex_com.top     \\'+ter)
  MMPBSA_send.write('             -rp $TOP_DIR/complex_rec.top     \\'+ter)
  MMPBSA_send.write('             -lp $TOP_DIR/complex_lig.top     \\'+ter)
  MMPBSA_send.write('             -y  $DYN_DIR/{}          \\'.format(name_file)+ter)
  MMPBSA_send.write('             >>  ./progres.log'+ter)
  MMPBSA_send.write('#'+ter)
  MMPBSA_send.write('rm -f -r _MMPBSA_*'+ter)
  MMPBSA_send.write('rm -f -r reference.frc'+ter)
  MMPBSA_send.write('#'+ter)
  MMPBSA_send.write('cp $TMPDIR/* $SGE_O_WORKDIR/'+ter)
  MMPBSA_send.write('#'+ter)
  MMPBSA_send.close()
  os.system('chmod u+x en_{}_mmpbsa_iqtc06'.format(num_poses))
  return

def GEN_MMPBSA_inputs_iqtc09 ( system_name,job_name,num_poses,name_file) :
# 
# iqtc09
#
  MMPBSA_send = open('en_{}_mmpbsa_iqtc09'.format(num_poses),'w')
  MMPBSA_send.write('#!/bin/bash'+ter)
  MMPBSA_send.write('######################################'+ter)
  MMPBSA_send.write('# Opcions i parametres del SGE        '+ter)
  MMPBSA_send.write('######################################'+ter)
  MMPBSA_send.write('# (1) Nom del treball                 '+ter)
  MMPBSA_send.write('#$ -N {}_{}                           '.format(job_name,num_poses)+ter)
  MMPBSA_send.write('# (2) Recursos sol.licitats           '+ter)
  MMPBSA_send.write('##$ -l h_rt                           '+ter)
  MMPBSA_send.write('##$ -l mem_free                       '+ter)
  MMPBSA_send.write('#$ -pe smp 14                         '+ter)
  MMPBSA_send.write('##$ -l exclusive=true                 '+ter)
  MMPBSA_send.write('# (3) Fitxers de sortida              '+ter)
  MMPBSA_send.write('#$ -cwd                               '+ter)
  MMPBSA_send.write('#$ -o {}.out                          '.format(job_name)+ter)
  MMPBSA_send.write('#$ -e {}.err                          '.format(job_name)+ter)
  MMPBSA_send.write('# (4) Envia un mail                   '+ter)
  MMPBSA_send.write('##$ -m e                              '+ter)
  MMPBSA_send.write('##$ -M jaime.rubio@ub.edu             '+ter)
  MMPBSA_send.write('######################################'+ter)
  MMPBSA_send.write('# Entorn de usuari                    '+ter)
  MMPBSA_send.write('######################################'+ter)
  MMPBSA_send.write('# Es carreguen els moduls             '+ter)
  MMPBSA_send.write('source /etc/profile                   '+ter)
  MMPBSA_send.write('env                                   '+ter)
  MMPBSA_send.write('module load amber/18_mod              '+ter)
  MMPBSA_send.write('######################################'+ter)
  MMPBSA_send.write('echo " SGE_O_WORKDIR : $SGE_O_WORKDIR"'+ter)
  MMPBSA_send.write('echo " nslots        : $NSLOTS"       '+ter)
  MMPBSA_send.write('echo " TMP DIR       : $TMPDIR"       '+ter)
  MMPBSA_send.write('######################################'+ter)
  MMPBSA_send.write('# Calcul                              '+ter)
  MMPBSA_send.write('######################################'+ter)
  MMPBSA_send.write('export RUN=/home/g6jaime/jaime/run    '+ter)
  MMPBSA_send.write('echo " RUN DIR       : $RUN"          '+ter)
  MMPBSA_send.write('echo "AMBERHOME   : $AMBERHOME "      '+ter)
  MMPBSA_send.write('source  $AMBERHOME/amber.sh           '+ter)
  MMPBSA_send.write('#'+ter)
  MMPBSA_send.write('cd ../fftop                           '+ter)
  MMPBSA_send.write('export TOP_DIR=`pwd`                  '+ter)
  MMPBSA_send.write('cd $SGE_O_WORKDIR                     '+ter)
  MMPBSA_send.write('cd ../../                             '+ter)
  MMPBSA_send.write('export DYN_DIR=`pwd`                  '+ter)
  MMPBSA_send.write('echo "TOP DIR   : $TOP_DIR "          '+ter)
  MMPBSA_send.write('echo "DYN DIR   : $DYN_DIR "          '+ter)
  MMPBSA_send.write('#'+ter)
  MMPBSA_send.write('cd $SGE_O_WORKDIR                     '+ter)
  MMPBSA_send.write('cp $SGE_O_WORKDIR/mmpbsa_py.inp $TMPDIR   '+ter)
  MMPBSA_send.write('cd $TMPDIR                            '+ter)
  MMPBSA_send.write('#'+ter)
  MMPBSA_send.write('mpirun -np $NSLOTS $AMBERHOME/bin/MMPBSA.py.MPI -O  \\'+ter)
  MMPBSA_send.write('             -i  ./mmpbsa_py.inp              \\'+ter)
  MMPBSA_send.write('             -o  ./mmpbsa_py.out              \\'+ter)
  MMPBSA_send.write('             -eo ./DeltaG_BySnaps.csv         \\'+ter)
  MMPBSA_send.write('             -sp $DYN_DIR/{}.top      \\'.format(system_name)+ter)
  MMPBSA_send.write('             -cp $TOP_DIR/complex_com.top     \\'+ter)
  MMPBSA_send.write('             -rp $TOP_DIR/complex_rec.top     \\'+ter)
  MMPBSA_send.write('             -lp $TOP_DIR/complex_lig.top     \\'+ter)
  MMPBSA_send.write('             -y  $DYN_DIR/{}          \\'.format(name_file)+ter)
  MMPBSA_send.write('             >>  ./progres.log'+ter)
  MMPBSA_send.write('#'+ter)
  MMPBSA_send.write('rm -f -r _MMPBSA_*'+ter)
  MMPBSA_send.write('rm -f -r reference.frc'+ter)
  MMPBSA_send.write('#'+ter)
  MMPBSA_send.write('cp $TMPDIR/* $SGE_O_WORKDIR/'+ter)
  MMPBSA_send.write('#'+ter)
  MMPBSA_send.close()
  os.system('chmod u+x en_{}_mmpbsa_iqtc09'.format(num_poses))
  return


def GEN_MMPBSA_inputs_iqtc03 ( system_name,job_name,num_poses,name_file ) :
# 
# iqtc03  g6.q 
#
  MMPBSA_send = open('en_{}_mmpbsa_iqtc03'.format(num_poses),'w')
  MMPBSA_send.write('#!/bin/bash'+ter)
  MMPBSA_send.write('######################################'+ter)
  MMPBSA_send.write('# Opcions i parametres del SGE        '+ter)
  MMPBSA_send.write('######################################'+ter)
  MMPBSA_send.write('# (1) Nom del treball                 '+ter)
  MMPBSA_send.write('#$ -N {}_{}                           '.format(job_name,num_poses)+ter)
  MMPBSA_send.write('# (2) Recursos sol.licitats           '+ter)
  MMPBSA_send.write('##$ -l h_rt                           '+ter)
  MMPBSA_send.write('##$ -l mem_free                       '+ter)
  MMPBSA_send.write('#$ -pe smp 27                         '+ter)
  MMPBSA_send.write('##$ -l exclusive=true                 '+ter)
  MMPBSA_send.write('# (3) Fitxers de sortida              '+ter)
  MMPBSA_send.write('#$ -cwd                               '+ter)
  MMPBSA_send.write('#$ -o {}.out                          '.format(job_name)+ter)
  MMPBSA_send.write('#$ -e {}.err                          '.format(job_name)+ter)
  MMPBSA_send.write('# (4) Envia un mail                   '+ter)
  MMPBSA_send.write('##$ -m e                              '+ter)
  MMPBSA_send.write('##$ -M jaime.rubio@ub.edu             '+ter)
  MMPBSA_send.write('######################################'+ter)
  MMPBSA_send.write('# Entorn de usuari                    '+ter)
  MMPBSA_send.write('######################################'+ter)
  MMPBSA_send.write('# Es carreguen els moduls             '+ter)
  MMPBSA_send.write('source /etc/profile                   '+ter)
  MMPBSA_send.write('env                                   '+ter)
  MMPBSA_send.write('module load amber/20_cuda9.0_ompi_gcc-5.5.0  '+ter)
  MMPBSA_send.write('######################################'+ter)
  MMPBSA_send.write('echo " SGE_O_WORKDIR : $SGE_O_WORKDIR"'+ter)
  MMPBSA_send.write('echo " nslots        : $NSLOTS"       '+ter)
  MMPBSA_send.write('echo " TMP DIR       : $TMPDIR"       '+ter)
  MMPBSA_send.write('######################################'+ter)
  MMPBSA_send.write('# Calcul                              '+ter)
  MMPBSA_send.write('######################################'+ter)
  MMPBSA_send.write('export RUN=/home/g6jaime/jaime/run    '+ter)
  MMPBSA_send.write('echo " RUN DIR       : $RUN"          '+ter)
  MMPBSA_send.write('echo "AMBERHOME   : $AMBERHOME "      '+ter)
  MMPBSA_send.write('source  $AMBERHOME/amber.sh           '+ter)
  MMPBSA_send.write('#'+ter)
  MMPBSA_send.write('cd ../fftop                           '+ter)
  MMPBSA_send.write('export TOP_DIR=`pwd`                  '+ter)
  MMPBSA_send.write('cd $SGE_O_WORKDIR                     '+ter)
  MMPBSA_send.write('cd ../../                             '+ter)
  MMPBSA_send.write('export DYN_DIR=`pwd`                  '+ter)
  MMPBSA_send.write('echo "TOP DIR   : $TOP_DIR "          '+ter)
  MMPBSA_send.write('echo "DYN DIR   : $DYN_DIR "          '+ter)
  MMPBSA_send.write('#'+ter)
  MMPBSA_send.write('cd $SGE_O_WORKDIR                     '+ter)
  MMPBSA_send.write('cp $SGE_O_WORKDIR/mmpbsa_py.inp $TMPDIR   '+ter)
  MMPBSA_send.write('cd $TMPDIR                            '+ter)
  MMPBSA_send.write('#'+ter)
  MMPBSA_send.write('mpirun -np $NSLOTS $AMBERHOME/bin/MMPBSA.py.MPI -O  \\'+ter)
  MMPBSA_send.write('             -i  ./mmpbsa_py.inp              \\'+ter)
  MMPBSA_send.write('             -o  ./mmpbsa_py.out              \\'+ter)
  MMPBSA_send.write('             -eo ./DeltaG_BySnaps.csv         \\'+ter)
  MMPBSA_send.write('             -sp $DYN_DIR/{}.top      \\'.format(system_name)+ter)
  MMPBSA_send.write('             -cp $TOP_DIR/complex_com.top     \\'+ter)
  MMPBSA_send.write('             -rp $TOP_DIR/complex_rec.top     \\'+ter)
  MMPBSA_send.write('             -lp $TOP_DIR/complex_lig.top     \\'+ter)
  MMPBSA_send.write('             -y  $DYN_DIR/{}          \\'.format(name_file)+ter)
  MMPBSA_send.write('             >>  ./progres.log'+ter)
  MMPBSA_send.write('#'+ter)
  MMPBSA_send.write('rm -f -r _MMPBSA_*'+ter)
  MMPBSA_send.write('rm -f -r reference.frc'+ter)
  MMPBSA_send.write('#'+ter)
  MMPBSA_send.write('cp $TMPDIR/* $SGE_O_WORKDIR/'+ter)
  MMPBSA_send.write('#'+ter)
  MMPBSA_send.close()
  os.system('chmod u+x en_{}_mmpbsa_iqtc03'.format(num_poses))
  return

def GEN_MMPBSA_inputs_slurm  ( system_name,job_name,num_poses,name_file,num_proc ) :
##
## Almeria : Slurm
#
  MMPBSA_send = open('en_{}_mmpbsa_slurm'.format(num_poses),'w')
  MMPBSA_send.write('#!/bin/bash '+ter)
  MMPBSA_send.write('#SBATCH -J {}_{}   '.format(job_name,num_poses)+ter)
  MMPBSA_send.write('#SBATCH slurm.err  '+ter)
  MMPBSA_send.write('#SBATCH slurm.out  '+ter)
  MMPBSA_send.write('#SBATCH -N 1       '+ter)
  MMPBSA_send.write('#SBATCH -n {}      '.format(num_proc)+ter)
  MMPBSA_send.write('#SBATCH --ntasks-per-node=1  '+ter)
  MMPBSA_send.write('#SBATCH --time=6-0           '+ter)
  MMPBSA_send.write('#'+ter)
  MMPBSA_send.write('export AMBERHOME={}'.format(AMBERHOME)+ter)
  MMPBSA_send.write('source {}'.format(AMBER_ENV)+ter)
  MMPBSA_send.write('#'+ter)
  if num_proc > 1 :
    MMPBSA_send.write('mpirun -np {} $AMBERHOME/bin/MMPBSA.py.MPI  -O  \\'.format(num_proc)+ter)
  else :
    MMPBSA_send.write('$AMBERHOME/bin/MMPBSA.py  -O     \\'+ter)
  MMPBSA_send.write('    -i  ./mmpbsa_py.inp          \\'+ter)
  MMPBSA_send.write('    -o  ./mmpbsa_py.out          \\'+ter)
  MMPBSA_send.write('    -eo ./DeltaG_BySnaps.csv     \\'+ter)
  MMPBSA_send.write('    -sp ../../{}.top     \\'.format(system_name)+ter)
  MMPBSA_send.write('    -cp ../fftop/complex_com.top \\'+ter)
  MMPBSA_send.write('    -rp ../fftop/complex_rec.top \\'+ter)
  MMPBSA_send.write('    -lp ../fftop/complex_lig.top \\'+ter)
  MMPBSA_send.write('    -y  ../../{}          \\'.format(name_file)+ter)
  MMPBSA_send.write('     >>  ./progres.log'+ter)
  MMPBSA_send.write('#'+ter)
  MMPBSA_send.write('rm -f -r _MMPBSA_*'+ter)
  MMPBSA_send.write('rm -f -r reference.frc'+ter)
  MMPBSA_send.write('#'+ter)
  MMPBSA_send.close()
  os.system('chmod u+x en_{}_mmpbsa_slurm'.format(num_poses))
  return

#================================================================
# _MAIN_work_
#================================================================
def prog_info () :
  t = datetime.now()
  format_date = "%d-%m-%Y %H:%M:%S"
  ti= t.strftime(format_date)
  print ( "... ... ... ... ... ... ... " )
  print ( "...     Step R3         ... " )
  print ( "... MMPB/GBSA OneLig    ... " )
  print ( "...   ( 23-01-2021 )    ... " )
  print ( "...       ........      ... " )
  print ( "... ... ... ... ... ... ... " )
  print ( " ",ti )
  print ( "... ... ... ... ... ... ... " )
  print ( " " )

if __name__ == '__main__':

# 12-08-2020
# - Possibility to use .top file to determine number of protein residues
# - iqtc09 files generated 
#232-01-2021
# - Reading number os snaps from  a cpptraj calculation

  amber_version = 20
  amber_dir     =  "/programas/amber20/amber" + str(amber_version)
  cpptraj       =  amber_dir + "/bin/cpptraj"



  main()


