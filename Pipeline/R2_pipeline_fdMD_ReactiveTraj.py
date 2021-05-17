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
1. use_pocket: True/False
        : True : The reference structure is defined by residues in a pocket 
2. file_pocket : One line with the residues of the protein defining the pocket
                 ( in AMBER Mask format ) 
3. XRay : False/ True 
        : False :  The reference structure is not the X-Ray Structure
        : True  :  The reference structure is     the X-Ray Structure

# If use_pocket = False AND XRay= False, the LAST structure of each trajectory is
  used as reference for the plots time/distance. Where distance is defined between 
  the position of C99/X99 in the last point and the position of C99/X99 along the trajectory. 

4. protein_pdb : Must be defined and ALL the trajectories must be superposed to this structure.
                 The protein structure to be used to analyse reactivity.
                 Must be the same used in the R1_fdMD_OneTrajREM_????.py calculations
                 Can include the ligand or not. If XRay=True The Ligand must be present. 
                 Take care. A TER line indicate end of protein. For dimers,trimer...this is
                 a problem. Remove the TER between monomers.

5. dis_lig_prot_min : A trajectory in which all the ligand atoms are farther than 
                      this value with respect to any protein atoms is a NON REACTIVE trajectory.

6. num_plots    : Number of different plots to be done.

7. dis_plots    : Limits for the axis Y .

8. snaps_by_one : Number of snapshots per ns .

9. num_ns_anal  : Number of ns to be analysed in depth.  

10. percent_anal : % of points in the analysed ns ( num_ns_anal ) that are allowed to be
                  greater than " dis_lig_prot_min " to stay the trajectory as REACTIVE.

11. prefix_filpdb : Prefix for the files containing the positions of the atoms with wdW repulsion.

12. name_atomrep  : Name of the atom with wdW repulsion.
        
"""
# global variables
lig_BS = dict()

###########MAIN
def main() :
  
  global  x_prot, y_prot, z_prot, at_name, res_name, res_num
  global  x_lig, y_lig, z_lig, atom_name_lig, resi_name_lig
  
  

  prog_info ()
  
  args = cmdlineparse()

  cha_pocket =      args.use_pocket
  cha_XRay   =      args.XRay
  file_pocket=      args.file_pocket 
  protein_pdb=      args.protein_pdb 

  XRay = False
  if cha_XRay   == 'True' :
    XRay = True

  use_pocket = False
  if cha_pocket == 'True' :
    use_pocket = True

  use_last = False 
  if ( not XRay )  and ( not use_pocket)  :
    use_last = True

  dis_lig_prot_min = float ( args.dis_lig_prot_min )
  num_plots    = int ( args.num_plots )
  dis_plots    =     ( args.dis_plots )
  snaps_by_one = int ( args.snaps_by_one )
  num_ns_anal  = int ( args.num_ns_anal  )
  percent_anal = int ( args.percent_anal )

  prefix_filrep=     ( args.prefix_filrep )
  name_atomrep =     ( args.name_atomrep  )
  
  top_name = args.top_name
  num_res_prot = int(args.num_res_prot)
  fragment = (args.fragment)
  dir_traj = (args.dir_traj)
  end_traj = args.end_traj
  
  if use_pocket :
    print ( " Use pocket to define distances  : ",use_pocket )
    print ( " Pocket Definition               : ",file_pocket )
  if XRay :
    print ( " Use XRay structure as reference : ",XRay )
    print ( " Using X-Ray structure           : ",protein_pdb )
  if use_last :
    print ( " Use Last Protein structure      : ",use_last )
    print ( " Using a Referece structure      : ",protein_pdb )

  print ( " Atom name to analyse            : ",name_atomrep )

  print ( " Minimum distance Ligand-Prot    : ",dis_lig_prot_min )
  print ( " Num Plots to be done            : ",num_plots )       
  print ( " Dis Plots  [x,y,..]             : ",dis_plots )       

  print ( " Snaps by nanosecond             : ",snaps_by_one )    
  print ( " Num ns to analyse               : ",num_ns_anal  )    
  print ( " PerCent out of the limits       : ",percent_anal )    

  print ( " Prefix added to pdb files       : ",prefix_filrep )

  dis_plots = dis_plots.split(",")
  for ndim in range (len(dis_plots)) :
    dis_plots[ndim] = float (dis_plots[ndim] )

  initial_dir = os.getcwd()
  print ( "  " )
  print ( " Initial Dir = ",initial_dir )
  print ( "  " )
#
# Directory for pictures
#
  dir_fig_dist = 'figures_distances'
  exist_dir = os.path.isdir(dir_fig_dist)
  if exist_dir :
    if os.listdir(dir_fig_dist) != [] :
      os.chdir(dir_fig_dist)
      os.system ( "rm * " )
      os.chdir('..')
    os.rmdir(dir_fig_dist)
  os.mkdir(dir_fig_dist)
#
# Directory for DCM pictures
#
  dir_fig_DCM  = 'figures_DCM'
  exist_dir = os.path.isdir(dir_fig_DCM)
  if exist_dir :
    if os.listdir(dir_fig_DCM) != [] :
      os.chdir(dir_fig_DCM)
      os.system ( "rm * " )
      os.chdir('..')
    os.rmdir(dir_fig_DCM)
  os.mkdir(dir_fig_DCM)
#
# Directory for DCM pictures
#
  dir_fig_LIE  = 'figures_LIE'
  exist_dir = os.path.isdir(dir_fig_LIE)
  if exist_dir :
    if os.listdir(dir_fig_LIE) != [] :
      os.chdir(dir_fig_LIE)
      os.system ( "rm * " )
      os.chdir('..')
    os.rmdir(dir_fig_LIE)
  os.mkdir(dir_fig_LIE)
#
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

  num_at_prot, num_at_lig, pos_name_atomrep, resi_name_lig, at_CA = read_prot (protein_pdb,name_atomrep)
  print ( "... There are = ",num_at_prot," PROTEIN Atoms " )
  if num_at_lig != 0 :
    print ( "... There are = ",num_at_lig, " LIGAND  Atoms " )
    print ( "... With Residue name = ",resi_name_lig )
    print ( "... And ",name_atomrep, " at position ",pos_name_atomrep+1 )
  else :
    if XRay :
      print ( "... Using XRAY = True ... but LIGAND is NOT present " )
      exit ()
#
# Directory for pictures  X-Ray
#
  if use_pocket or ( num_at_lig != 0 and XRay ) :
    dir_fig_xray = 'figures_distances_XRay'
    exist_dir = os.path.isdir(dir_fig_xray)
    if exist_dir :
      if os.listdir(dir_fig_xray) != [] :
        os.chdir(dir_fig_xray)
        os.system ( "rm * " )
        os.chdir('..')
      os.rmdir(dir_fig_xray)
    os.mkdir(dir_fig_xray)

#file_delete_nc = ("Pre_Grafic_{}.out".format(nombres[KUAL]))

  file_delete_nom = "delete_Non_Reactive"
  file_delete = open (file_delete_nom,"w")

  if use_pocket :
    print ( " ... Pocket ACTIVE for reference definition " ) 
    if XRay :
      print ( " ... Pocket ACTIVE thus XRay turned OFF " )
      XRay = False
    res_pocket  = read_pocket ( file_pocket ) 
    x_ref,y_ref,z_ref = calc_cent_geom (res_pocket,num_at_prot,x_prot,y_prot,z_prot,at_name,res_num)
    print ( "    ... Pocket REF  {:7.3f} {:7.3f} {:7.3f} ".format(x_ref,y_ref,z_ref)) 

  if num_at_lig != 0 and XRay : 
    x_ref = x_lig[pos_name_atomrep]
    y_ref = y_lig[pos_name_atomrep]
    z_ref = z_lig[pos_name_atomrep]
    print ( "... XRay REF  {:7.3f} {:7.3f} {:7.3f} ".format(x_ref,y_ref,z_ref)) 

  y_ave_min      = 99999.9
  y_ave_min_xray = 99999.9
  lig_ave_min      = 0

  lig_ave_min_xray = 0
  num_react_xray   = 0
  num_reac         = 0
  num_LIG          = 0
  react_xray       = []
  traj_react_xray  = []
  reac_LIG         = []
  list_LIG         = []
  reactive_LIGNUM = []

  print_ref = True
  lista = os.listdir( os.getcwd() )
  for pdb in lista:
    if prefix_filrep in pdb:
      print("                     ")
      print("... Pocessing File  :",pdb)

      pdb_split = pdb.split("_")
      lig_num = pdb_split [1]
      list_LIG.append(lig_num)
      num_LIG += 1

      x  = []
      y  = []
      z  = []
      num_snaps = read_pdb_c99(pdb,name_atomrep,x,y,z)

      num_ns = num_snaps / snaps_by_one
      num_snaps_NOanal = num_snaps - snaps_by_one * num_ns_anal

# always use last to annalise 

      x_last = x[-1]
      y_last = y[-1]
      z_last = z[-1]
      if  print_ref  :
        print ( "    ... Last selected ALWAYS  {:7.3f} {:7.3f} {:7.3f} ".format(x_last,y_last,z_last)) 
        print_ref = False

      min_lig_prot, num_near_ref = dis_min_prot (num_at_prot,dis_lig_prot_min,x_last,y_last,z_last)

      dis_at = []
      d_max = 0.0
      for i in range( num_snaps ) :
        x_tm = x [i]
        y_tm = y [i]
        z_tm = z [i]
        x2 = ( x_tm - x_last ) * ( x_tm - x_last )
        y2 = ( y_tm - y_last ) * ( y_tm - y_last )
        z2 = ( z_tm - z_last ) * ( z_tm - z_last )
        d = sqrt ( x2 + y2 + z2 )
        if d > d_max :
          d_max = d
        dis_at.append(d)

      print ( "    ... D_max ( {}_last --> {}_dyn  )    = {:10.4f} ".format(name_atomrep,name_atomrep,d_max) )
      print ( "    ... D_min ( {}_last --> Protein )     = {:10.4f} ".format(name_atomrep,min_lig_prot) )
      print ( "    ... Number of PROTEIN atoms near LAST  = {:6d}  ".format(num_near_ref) )

      if use_pocket or ( num_at_lig != 0 and XRay ) :
        dis_at_xray = []
        d_min_xray = 99999.0
        for i in range( num_snaps ) :
          x_tm = x [i]
          y_tm = y [i]
          z_tm = z [i]
          x2 = ( x_tm - x_ref ) * ( x_tm - x_ref )
          y2 = ( y_tm - y_ref ) * ( y_tm - y_ref )
          z2 = ( z_tm - z_ref ) * ( z_tm - z_ref )
          d = sqrt ( x2 + y2 + z2 )
          if d < d_min_xray :
            d_min_xray = d
            i_min_xray = i
          dis_at_xray.append(d)
          dis_xray2last = d
        print ( "    ... D_min ( X-Ray --> {}_dyn )        = {:10.4f} in snapshot {} ".format(name_atomrep,d_min_xray,i_min_xray) )
        print ( "    ... D_min ( X-Ray --> LAST_dyn )       = {:10.4f}  ".format(dis_xray2last) )
        if dis_xray2last  <= dis_lig_prot_min :
          num_react_xray += 1
          react_xray.append(dis_xray2last)
          traj_react_xray.append(lig_num)

      if XRay :
        reftraj_to_delet =  'first_lig_' + lig_num + '.nc'
        reftraj_to_delet =  ' rm -f ' + reftraj_to_delet + '\n'
        file_delete.write( reftraj_to_delet )

# FIRST : a minimum distance to Protein

      if min_lig_prot >= dis_lig_prot_min :

        print("    --> File :",pdb, ' Is a NON reactive trajectory withing ',dis_lig_prot_min, ' A ')

        name_to_delete     = 'lig_' + lig_num + '*'
        name_nc_to_delete  = 'lig_' + lig_num + '.nc' 
        name_pdb_to_delete = 'lig_' + lig_num + prefix_filrep + '.pdb' 
        name_end_to_delete = 'lig_' + lig_num + '_LAST.pdb'
        name_csv_to_delete = 'lig_' + lig_num + '_disat.csv'
        name_sa_to_delete  = 'lig_' + lig_num + '_SASA.out'

        print("       --> Files To Delete :",name_to_delete)

        name_nc_to_delete  = ' rm -f ' + name_nc_to_delete  + '\n'
        name_pdb_to_delete = ' rm -f ' + name_pdb_to_delete + '\n'
        name_end_to_delete = ' rm -f ' + name_end_to_delete + '\n'
        name_csv_to_delete = ' rm -f ' + name_csv_to_delete + '\n'
        name_sa_to_delete  = ' rm -f ' + name_sa_to_delete  + '\n'

        file_delete.write( name_nc_to_delete  )
        file_delete.write( name_pdb_to_delete )
        file_delete.write( name_end_to_delete )
        file_delete.write( name_csv_to_delete )
        file_delete.write( name_sa_to_delete  )

      else :
        num_no_react_pro   = 0
        num_snaps_work =  snaps_by_one * num_ns_anal
        for i in range ( num_snaps_NOanal,num_snaps ) :
          x_tm = x [i]
          y_tm = y [i]
          z_tm = z [i]
          min_snap_prot,num_near_atom = dis_min_prot (num_at_prot,dis_lig_prot_min,x_tm,y_tm,z_tm)
          if min_snap_prot >= dis_lig_prot_min :
            num_no_react_pro  += 1
        print ( "    ... {}-PRO : There are {} snaps Farther than {} in a total of {} ".format(name_atomrep,num_no_react_pro,dis_lig_prot_min,num_snaps_work ) )
        per_cent_pro = 100 * num_no_react_pro / num_snaps_work

        num_no_react_c99   = 0
        for i in range ( num_snaps_NOanal,num_snaps ) :
          if dis_at[i] >= dis_lig_prot_min :
            num_no_react_c99  += 1
        print ( "    ... {}-{} : There are {} snaps Farther than {} in a total of {} ".format(name_atomrep,name_atomrep,num_no_react_c99,dis_lig_prot_min,num_snaps_work ) )
        per_cent_c99 = 100 * num_no_react_c99 / num_snaps_work

        if  ( per_cent_pro  >= percent_anal ) or ( per_cent_c99  >= percent_anal ) :
          if  per_cent_pro  >= percent_anal :
            print('    --> File : {} Is a NON-REACTIVE trajectory with {} % OUT [ Protein ]'.format(pdb,per_cent_pro) )
          if  per_cent_c99  >= percent_anal :
            print('    --> File : {} Is a NON-REACTIVE trajectory with {} % OUT [   {}    ]'.format(pdb,per_cent_c99,name_atomrep) )

          name_to_delete     = 'lig_' + lig_num + '*'
          name_nc_to_delete  = 'lig_' + lig_num + '.nc' 
          name_pdb_to_delete = 'lig_' + lig_num + prefix_filrep + '.pdb' 
          name_end_to_delete = 'lig_' + lig_num + '_LAST.pdb'
          name_csv_to_delete = 'lig_' + lig_num + '_disat.csv'

          print("       --> Files To Delete :",name_to_delete)

          name_nc_to_delete  = ' rm -f ' + name_nc_to_delete  + '\n'
          name_pdb_to_delete = ' rm -f ' + name_pdb_to_delete + '\n'
          name_end_to_delete = ' rm -f ' + name_end_to_delete + '\n'
          name_csv_to_delete = ' rm -f ' + name_csv_to_delete + '\n'

          file_delete.write( name_nc_to_delete  )
          file_delete.write( name_pdb_to_delete )
          file_delete.write( name_end_to_delete )
          file_delete.write( name_csv_to_delete )
        else :
          print('    --> File : {} Is a REACTIVE trajectory withing {} A of Protein'.format(pdb,dis_lig_prot_min) )
          num_reac += 1
          reac_LIG.append(pdb)
          reactive_LIGNUM.append(lig_num) 
         
#
# SASA calculation 
      name_file_SASA  = 'lig_' + lig_num + '_SASA.out'
      sasa,pc_sasa = GET_sasa (name_file_SASA)
      print ( "    ... Ligand {:5} : SASA = {:10.1f}  %SASA = {:10.1f} ".format(lig_num,sasa,pc_sasa) )
# ... END  SASA ...

      x_din = np.linspace( 0, num_ns, num_snaps )
      y_dis = np.array ( dis_at )
    
      print ( "    ... Analysing ",num_ns_anal," ns with ",snaps_by_one * num_ns_anal," snaps in a total of ",num_snaps)
      y_ave, sd_y, y_max   = stat_x_ns (num_snaps,num_snaps_NOanal,y_dis )
      print ( "        ... y_ave ({:10.3f} ),sd( {:10.3f} ) y_max( {:10.3f} ) [ LAST ]" .format(y_ave, sd_y, y_max )) 
      if y_ave <= y_ave_min :
        y_ave_min = y_ave
        lig_ave_min = lig_num

      if use_pocket or ( num_at_lig != 0 and XRay ) :
        y_dis_xray = np.array ( dis_at_xray )
        y_ave_xray, sd_y_xray, y_max_xray   = stat_x_ns (num_snaps,num_snaps_NOanal,y_dis_xray )
        print ( "        ... y_ave ({:10.3f} ),sd( {:10.3f} ) y_max( {:10.3f} ) [ XRay ]" .format(y_ave_xray, sd_y_xray, y_max_xray )) 
        if y_ave_xray <= y_ave_min_xray :
          y_ave_min_xray   = y_ave_xray
          lig_ave_min_xray = lig_num


      file_excel_disat = 'lig_'+lig_num+'_disat.csv'
      file_excel       = open (file_excel_disat,"w")

      for rowx in range ( num_snaps ) :
        strw = ( "{:15.9} {:15.9}".format(x_din[rowx],y_dis[rowx]) )
        strw = str ( strw ) 
        strw  = strw + '\n'
        file_excel.write(strw)

      file_plot        = 'lig_'+lig_num
    
      for num_plt in range ( num_plots ) :

        file_plot_disat  = 'lig_'+lig_num+'_disat'+str(dis_plots[num_plt])+'.png'
        distan_to_plot = dis_plots[num_plt]

        #"""
        plt.xlim(0,num_ns)
        plt.ylim(0,distan_to_plot)

        plt.title(file_plot)
        plt.xlabel("t(ns)")
        plt.ylabel("Distance (Angstrom)")
#       plt.plot (x_din,y_dis,'g', label='_nolegend_')
        plt.plot (x_din,y_dis,'g', label='Dist_Last')
        plt.legend()
        plt.savefig (file_plot_disat, format='png',dpi=600)
        plt.show()
        plt.close()

        move = 'mv ' + file_plot_disat + ' ./' + dir_fig_dist + '/'
        os.system ( move  )
        #""" 
    
      if use_pocket or ( num_at_lig != 0 and XRay ) :

        file_excel_disat_xray = 'lig_'+lig_num+'_disat_xray.csv'
        file_excel_xray       = open (file_excel_disat_xray,"w")

        for rowx in range ( num_snaps ) :
          strw = ( "{:15.9} {:15.9}".format(x_din[rowx],y_dis_xray[rowx]) )
          strw = str ( strw ) 
          strw  = strw + '\n'
          file_excel_xray.write(strw)

        file_plot        = 'lig_'+lig_num
        for num_plt in range ( num_plots ) :

          distan_to_plot = dis_plots[num_plt]
          file_plot_disat  = file_plot + '_disat_xray' + str(distan_to_plot) + '.png'

          #"""
          plt.xlim(0,num_ns)
          plt.ylim(0,distan_to_plot)

          plt.title(file_plot)
          plt.xlabel("t(ns)")
          plt.ylabel("Distance (Angstrom)")
#         plt.plot (x_din,y_dis,'g', label='_nolegend_')
          plt.plot (x_din,y_dis_xray,'g', label='Dist_X-Ray')
          plt.legend()
          plt.savefig (file_plot_disat, format='png',dpi=600)
          plt.show()
          plt.close()

          move = 'mv ' + file_plot_disat + ' ./' + dir_fig_xray + '/'
          os.system ( move  )
        #""" 

  os.system ('chmod u+x {}'.format(file_delete_nom) )
  run_delete_nom = "./" + file_delete_nom
  os.system('chmod u+x {} '.format(run_delete_nom))

  print ( "..." ) 
  print ( "... There are  = ",num_reac," REACTIVE trajectories " ) 
  if use_pocket or ( num_at_lig != 0 and XRay ) :
    if num_react_xray != 0 :
      print ( '... There are  = {:3d} trajectories NEAR Pocket-XRay ( {} A )'.format(num_reac_xray,dis_lig_prot_min) ) 
      for i in range (num_react_xray) :
        print ( '... ...  LIG : {_3d} -->  D(X-Ray-Last) : {} '.format(react_xray[i],traj_react_xray[i]) ) 
    else :
      print ( '... NO trajectories NEAR Pocket-XRay ( {} A )'.format(dis_lig_prot_min)  ) 
  print ( "..." ) 

  if num_reac != 0 :
    print ( ">>> Processing = ",num_reac," REACTIVE trajectories <<<" ) 

    pos_LIG_line = ''
    space        = ' '
    reac_LIG.sort()
    list_LIG.sort()
    for lig in range (num_reac) :
      pdb          = reac_LIG [ lig ]
      pdb_split    = pdb.split("_")
      lig_num      = pdb_split [1]
      pos_LIG_line += lig_num + space
    print ( ' ... Active LIGANDS  : ' )
    print ( ' ...        {:>6} '.format(pos_LIG_line   ) )
    print ("                     ")

    print ( " 1 --- Determine Binding Sites ---")
    for lig in range(num_reac) :

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

      pdb          = reac_LIG [ lig ]
      pdb_split    = pdb.split("_")
      lig_num      = pdb_split [1]
      name_pdb_end = 'lig_' + lig_num + '_LAST.pdb'
      print("... Pocessing File  :",name_pdb_end)
      num_at_prot, num_at_lig, pos_name_atomrep, resi_name_lig, at_CA = read_prot (name_pdb_end,name_atomrep)
# Binding Site
      calc_BindingSite (lig_num,num_at_prot,num_at_lig,dis_lig_prot_min,at_CA)

# Distance to the Protein Geometric Center
    print ( " 2 --- Plot Distances LIG-ProteinCM  ---")
    for lig in range(num_LIG) :
      lig_num             = list_LIG [lig]
      name_file_data_DCM  = 'lig_'+lig_num+'_DCM.dat'
      print("... Pocessing File  :",name_file_data_DCM)
      file_data_DCM       = open (name_file_data_DCM,'r')
      plot_DisLigProtCenter (lig_num,file_data_DCM,num_ns,num_snaps,dir_fig_DCM)

# Evolution of LIE interaction energy            
    print ( " 3 --- Plot LIE Interaction Energy   ---")
    lig_num             = list_LIG [0]
    name_file_data_LIE  = 'lig_'+lig_num+'_LIE.dat'
    file_data_LIE       = open (name_file_data_LIE,'r')
    num_lines = 0
    for lin in file_data_LIE :
      num_lines += 1
    num_lines -= 1
    if num_lines != num_snaps :
      num_snaps_LIE = num_lines
    else :
      num_snaps_LIE = num_snaps
    for lig in range(num_LIG) :
      lig_num             = list_LIG [lig]
      name_file_data_LIE  = 'lig_'+lig_num+'_LIE.dat'
      print("... Pocessing File  :",name_file_data_LIE)
      file_data_LIE       = open (name_file_data_LIE,'r')
      plot_LIE_IntEnergy  (lig_num,file_data_LIE,num_ns,num_snaps_LIE,dir_fig_LIE)

    # Print the next step.
    print ( " 4 --- Create LAST PDBs and run SASA ---")
    
    # Calcualte the number of pictures to analyse.
    num_pictures_analyse = num_snaps - (num_ns_anal*snaps_by_one)
    
    # If there is no reactive ligand.
    if len(reac_LIG) == 0:
      # End the main().
      return
    
    # For each ligand in the reactie trajectories.
    for lig in range(num_reac):
        # get the ligand number of each reactive trajectory.
        lig_num = reactive_LIGNUM[lig]

        # Get the ligand name of the .nc (trajectroy) file.
        name_nc_lig = "lig_" + str(lig_num) + ".nc"
        
        # Print which file we are going to process.
        print("Processing file : lig_" + lig_num + "_last_100_SASA.out")
        #os.system("cd lig_" + lig_num + "_SASA")
        
        # Create the input file to create the LAST pdb.
        file_cpptraj_inp  = "Jinptrj_ligands_SASA"
        
        # Create the output file.
        file_cpptraj_out  = "Jinptrj_ligands_SASA.out"
        
        # Open the input file.
        file_cpptraj      = open (file_cpptraj_inp,"w")
        
        # Define the parm of One ligand and no water.
        parm_OneLigNW = "[OneLigNoWat]"
        
        # Define a variable to add the final line jump.
        cc = "\n"
        
        # Create line to execute cpptraj and save the information in the output file and execute until EOF.
        line = " " + cpptraj + " >  " + file_cpptraj_out + " << EOF " + "\n"      
        # Write line into the input file.
        file_cpptraj.write(line)
        
        # Create parm line to read the topology.
        line = " parm " + "OneLig.NoWat." + top_name + " " + parm_OneLigNW + "\n"
        # Write line into the input file.
        file_cpptraj.write(line)  

        # Add one empty line
        line = " " + "\n"
        # Write line into the input file.
        file_cpptraj.write(line) 
        
        # Get the trajin file and put the number of snaps we want to analyse (last 100ns).
        line_inp = " trajin " + name_nc_lig + " " + str(num_pictures_analyse) + " " + str(num_snaps) + " " + str(snaps_by_one) + " parm " + parm_OneLigNW + " " + cc
        # Write line into the input file.
        file_cpptraj.write(line_inp)
        
        # Get the trajout line and the name where we save the SASA last pdb with the last 100ns snaps.
        line_out = " trajout lig_" + lig_num + "_SASA_last_100" + ".pdb pdb " + cc
        # Write line into the input file.
        file_cpptraj.write(line_out)          
        
        # Line to run the input file.
        line_run = " run " + cc + cc
        # Write line into the input file.
        file_cpptraj.write(line_run)
        
        # Line that tells to the cpptraj that is the end.
        line = "EOF"
        # Write line into the input file.
        file_cpptraj.write(line)
        
        # Close the file and save the written data.
        file_cpptraj.close()  

        # Create the running command to run the file created before.
        run_cpptraj = "./" + file_cpptraj_inp
        
        # Change the mode of the file created before (to be able to execute it)
        os.system('chmod u+x {} '.format(run_cpptraj) )
        
        # Run the file created before.
        os.system(run_cpptraj) 
    
        # Create the input file to Calculate the SASA.
        file_cpptraj_inp  = "Kinptrj_SASA_calculation"
        
        # Create the output file to create the SASA.
        file_cpptraj_out  = "Kinptrj_SASA_calculation.out"
        
        # Open the input file.
        file_cpptraj      = open (file_cpptraj_inp,"w")
        
        # Define the parm of One ligand and no Water.
        parm_OneLigNW = "[OneLigNoWat]"
        
        # Define a jump line (to add to the end of each line).
        cc = "\n"
        
        # Create line to execute cpptraj and save the information in the output file and execute until EOF.
        line = " " + cpptraj + " >  " + file_cpptraj_out + " << EOF " + "\n"
        # Write line into the input file.
        file_cpptraj.write(line)
        
        # Addd the line to read the topology file.
        line = " parm " + "OneLig.NoWat." + top_name + " " + parm_OneLigNW + "\n"
        # Write line into the input file.
        file_cpptraj.write(line)  
        
        # Add an empty line.
        line = " " + "\n"
        # Write line into the input file.
        file_cpptraj.write(line) 

        # Define a variable to know how many residues the complex has.
        num_res_tot = num_res_prot + 1

        # Define the name of the output SASA file to save the data.
        name_out_sasa = "lig_" + str(lig_num) + "_last_100" + "_SASA.out"

        # Craete the input line of trajin with the pdb created before (with all the snaps of last 100ns).
        line_inp = " trajin lig_" + lig_num + "_SASA_last_100" + ".pdb" + " parm " +         parm_OneLigNW + " " + cc
        
        # Write line into the input file.
        file_cpptraj.write(line_inp)
            
        # Add the command to calculate the SASA of the Complex.
        line = " molsurf  out  " + name_out_sasa + ' :1-' + str(num_res_tot) + "\n"
        # Write line into the input file.
        file_cpptraj.write(line)  
        
        # Add the command to calculate the SASA of the Receptor.
        line = " molsurf  out  " + name_out_sasa + ' :1-' + str(num_res_prot) + "\n"
        # Write line into the input file.
        file_cpptraj.write(line)
        
        # Add the command to calculate the SASA of the Ligand.
        line = " molsurf  out  " + name_out_sasa + ' :' + str(num_res_tot) + "\n"
        # Write line into the input file.
        file_cpptraj.write(line)  
        
        # Add the run command to executed all.
        line = " run " + "\n"
        # Write line into the input file.
        file_cpptraj.write(line)  
        
        # Add an empty line.
        line = " " + "\n"
        # Write line into the input file.
        file_cpptraj.write(line)  
        
        # Line that tells to the ccptraj that it is the end.
        line = "EOF" + "\n"
        # Write line into the input file.
        file_cpptraj.write(line)  
 
        # Close the input file.
        file_cpptraj.close() 
        
        # Define the variable that has the execution bash command.
        run_cpptraj = "./" + file_cpptraj_inp
        
        # Change the mode of the input file to be able to execute it.
        os.system  ('chmod u+x {} '.format(run_cpptraj) )
        
        # Execute the input file to calcualte the SASA for the Complex, Receptor and Ligand.
        os.system (run_cpptraj) 
        
        #os.system("mv *_SASA ")
    # Print that the SASA calculation has ended.
    print ( ">>> END Processing SASA calculation <<<" ) 
    
    # Print the next step done.
    print ( " 5 --- Find Hydrogen bonds during last 100ns and all trajectory ---")
    
    # Define the number of pictures we need to analyse.
    num_pictures_analyse = num_snaps - (num_ns_anal*snaps_by_one)        
    
    # For each ligan din the number of reactive trajectories.
    for lig in range(num_reac):
        # Define a variable with the ligand number.
        lig_num = reactive_LIGNUM[lig]
        
        # Define a variable with the ligand file .nc.
        name_nc_lig = "lig_" + str(lig_num) + ".nc"
        
        # Print which ligand we are processing.
        print("Processing file : lig_" + lig_num + "_avg-hb_atoms_last100.dat")
              
        # Create the input file for processing the Hydrogen bonds.
        file_cpptraj_inp  = "Linptrj_ligands_HB"
        
        # Create the output file for processing the hydrogen bonds.
        file_cpptraj_out  = "Linptrj_ligands_HB.out"
        
        # Open the input file for processing the Hydrogen Bonds.
        file_cpptraj      = open (file_cpptraj_inp,"w")
        
        # Define the parm that has One ligand and No water.
        parm_OneLigNW = "[OneLigNoWat]"
        
        # Define the variable that allow us to add a jump line at the end of each line.
        cc = "\n"
        
        # Create line to execute cpptraj and save the information in the output file and execute until EOF.
        line = " " + cpptraj + " >  " + file_cpptraj_out + " << EOF " + "\n"
        # Write line into the input file.
        file_cpptraj.write(line)
        
        # Add the line to read the topology.
        line = " parm " + "OneLig.NoWat." + top_name + " " + parm_OneLigNW + "\n"
        # Write line into the input file.
        file_cpptraj.write(line)  
        
        # Add one empty line.
        line = " " + "\n"
        # Write line into the input file.
        file_cpptraj.write(line) 
        
        # Add the line needed to define the input of the cpptraj (read the input file with only the snaps we want to analyse).
        line_inp = " trajin " + name_nc_lig + " " + str(num_pictures_analyse) + " " + str(num_snaps) + " 1 parm " + parm_OneLigNW + " " + cc
        # Write line into the input file.
        file_cpptraj.write(line_inp)
        
        # Get the hydrogen bonds and save them into the file called lig_"lig_num"_avg-hb_atoms_last100.dat.
        line_out = " hbond :* avgout lig_" + lig_num + "_avg-hb_atoms_last100.dat printatomnum nointramol" + cc
        # Write line into the input file.
        file_cpptraj.write(line_out)          
        
        # Line that tells to the cpptraj qhat has to run.
        line_run = " run " + cc + cc
        # Write line into the input file.
        file_cpptraj.write(line_run)
        
        # Line that tells to the cpptraj that it's the end.
        line = "EOF"
        # Write line into the input file.
        file_cpptraj.write(line)
        
        # Close the input file.
        file_cpptraj.close()  

        # Define the variable that has the execution bash command.
        run_cpptraj = "./" + file_cpptraj_inp
        
        # Change the mode of the input file to be able to execute it.
        os.system('chmod u+x {} '.format(run_cpptraj) )
        
        # Run the input file using bash.
        os.system(run_cpptraj) 
        
        # Run a grep command to only get the Hydrogen bonds that imply the fragment analysed.
        os.system("grep " + fragment + " lig_" + lig_num + "_avg-hb_atoms_last100.dat > lig_" + lig_num + "_avg-hb_fragment_last100.dat")
        
        # Print the next file to be processed.
        print("Processing file : lig_" + lig_num + "_avg-hb_atoms_all.dat")
        
        ## Calculate Hydrogen bonds for all the trajectory
        
        # Create the input file to Calculate the Hydrogen bonds during all the trajectory.
        file_cpptraj_inp  = "Minptrj_ligands_HB_all"
        
        # Create the output file.
        file_cpptraj_out  = "Minptrj_ligands_HB_all.out"
        
        # Open th einput file.
        file_cpptraj      = open (file_cpptraj_inp,"w")
        
        # Define the parm of One ligand and no Water.
        parm_OneLigNW = "[OneLigNoWat]"
        
        # Define the variable that allow us to add a jump line at the end of each line.
        cc = "\n"
        
        # Create line to execute cpptraj and save the information in the output file and execute until EOF.
        line = " " + cpptraj + " >  " + file_cpptraj_out + " << EOF " + "\n"
        # Write line into the input file.
        file_cpptraj.write(line)
        
        # Add the line to read the topology.
        line = " parm " + "OneLig.NoWat." + top_name + " " + parm_OneLigNW + "\n"
        # Write line into the input file.
        file_cpptraj.write(line)  
        
        # Add an empty line.
        line = " " + "\n"
        # Write line into the input file.
        file_cpptraj.write(line) 

        # Define the number of residues we have in the complex.
        num_res_tot = num_res_prot + 1
        
        # Create the input line that allow us to run along all the trajectory.
        line_inp = " trajin lig_" + lig_num + ".nc "  + cc
        # Write line into the input file.
        file_cpptraj.write(line_inp)
            
        # Create the line needed to calculate the hydorgen bonds for each atom number and save them into the file called lig_"lig_num"_avg-hb_atoms_all.dat
        line = " hbond :* avgout lig_" + lig_num + "_avg-hb_atoms_all.dat printatomnum nointramol" + "\n"
        # Write line into the input file.
        file_cpptraj.write(line)  
        
        # Line that tells to cpptraj what has to be executed.
        line = " run " + "\n"
        # Write line into the input file.
        file_cpptraj.write(line)  
        
        # Add an empty line.
        line = " " + "\n"
        # Write line into the input file.
        file_cpptraj.write(line)  
        
        # Line that tells that is the end of the execution.
        line = "EOF" + "\n"
        # Write line into the input file.
        file_cpptraj.write(line)  
 
        # Close input file.
        file_cpptraj.close() 
        
        # Define the variable that has the execution bash command.
        run_cpptraj = "./" + file_cpptraj_inp
        
        # Change the mode of the input file to be able to run with a bash command.
        os.system  ('chmod u+x {} '.format(run_cpptraj) )
        
        # Run the input file using bash.
        os.system (run_cpptraj) 
        
    # Print the line that tells us that we have ended the Hydrogen bonds Calculation.
    print ( ">>> END Processing Hydrogen Bonds Calculation <<<" ) 
         
    # Print the next step.
    print ( " 6 --- Superpose the Binding Sites pockets and calculate the RMSD of each ligand alone ---")

    # Define the number of pictures that we want to analyse.
    num_pictures_analyse = num_snaps - (num_ns_anal*snaps_by_one)

    # Create a folder called pictures_RMSD_ligand.
    os.system("mkdir pictures_RMSD_ligand") 
    
    # For each ligand in the number of reactive trajetcories.
    for lig in range(num_reac):
        # Get the ligand number and define a variable with the ligand number.
        lig_num = reactive_LIGNUM[lig]

        # Print which ligand/file we are processing.
        print("Processing file : lig_" + lig_num + "_Superpose_BS.nc")

        # Create an input file to Superpose the BS.
        file_cpptraj_inp  = "Ninptrj_Superpose_BS"
        
        # Create the output file to Superpose the BS.
        file_cpptraj_out  = "Ninptrj_Superpose_BS.out"
        
        # Open the input file.
        file_cpptraj      = open (file_cpptraj_inp,"w")
        
        # Define the parm that ha sNo water and the BS.
        parm_W = "[NoWat_BS]"
        
        # Define a variable that allow us to add the jump of line at the end of each added line.
        cc = "\n"

        # Create line to execute cpptraj and save the information in the output file and execute until EOF.
        line = " " + cpptraj + " >  " + file_cpptraj_out + " << EOF " + "\n"
        # Write line into the input file.
        file_cpptraj.write(line)
        
        # Add the line to read the topology (by cpptraj).
        line = " parm "+ top_name + " " + parm_W + "\n"
        # Write line into the input file.
        file_cpptraj.write(line)
        
        # Add one empty line.
        line = " " + "\n"
        # Write line into the input file.
        file_cpptraj.write(line)

        # Iterate the number of trajectories we have (1 -> 200ns, 2 -> 400ns, 3 -> 600ns, etc).
        for i in range(1, int(end_traj) + 1):
          # Define the name of the input coordinates file (.nc).
          name_nc_1 = dir_traj + top_name[0:-4] + "_" + str(i) + "_dyn.nc"
          
          # Using the input file define the trajin command to be able to read all the trajectory (depending on the ns we have).
          line_inp = " trajin " + name_nc_1 + " 1 last 1 parm " + parm_W + " " + cc
          # Write line into the input file.
          file_cpptraj.write(line_inp)
        
        # Add the line that tells to cpptraj to do the autoimage.
        line_inp = " autoimage " + cc
        # Write line into the input file.
        file_cpptraj.write(line_inp)
       
        # Add the line to remove Waters, Na+ and Cl- ions of all the trajectory.
        line_inp = " strip :WAT,Na+,Cl-" + cc
        # Write line into the input file.
        file_cpptraj.write(line_inp) 

        # Add the line to strip all the ligands that are not the one we are processing. 
        line_inp = " strip !:1-" + str(num_res_prot) + "," + lig_num + " outprefix NoWat_BS " + cc
        # Write line into the input file.
        file_cpptraj.write(line_inp)
        
        # Define a variable to add the residues involved in the interaction of the BS. The BS of the ligand we are using.
        result = ":"
        
        # Iterate all the ligands
        for l in lig_BS:
          # If the ligand is the one analysed.
          if (l == ("lig_" + lig_num)):
            # Iterate the all the residues the ligand analysed has in his BS.
            for residue in lig_BS[l]:
              # Add the residues in the result variable.
              result += residue
              
              # Add a comma in the result variable.
              result += ","
        
        # Catch only the result till the final element -1 due that it has one last comma (that we don't want to have).
        result = result[:-1]
        
        # Create the input line we want to execute with cpptraj to superpose the BS during all the trajectory.
        line_inp = " rms first out lig_" + lig_num + "_Superpose_BS.dat " + result + "@CA " + cc
        # Write line into the input file.
        file_cpptraj.write(line_inp)

        # Save the trajectory using trajout command of cpptraj.
        line_out = " trajout lig_" + lig_num + "_Superpose_BS.nc" + cc
        # Write line into the input file.
        file_cpptraj.write(line_out)
        
        # Tell to the cpptraj what has to run.
        line_run = " run " + cc + cc
        # Write line into the input file.
        file_cpptraj.write(line_run)
        
        # Line that tells that is the end of the file.
        line = "EOF"
        # Write line into the input file.
        file_cpptraj.write(line)

        # Close the input file.
        file_cpptraj.close()

        # Define the variable that has the execution bash command.
        run_cpptraj = "./" + file_cpptraj_inp
        
        # Change the mode of the input file to be able to execute using bash command.
        os.system('chmod u+x {} '.format(run_cpptraj) )
        
        # Execute the input file using bash.
        os.system(run_cpptraj)

        # Create the input file to get the final pdb.
        file_cpptraj_inp  = "Pinptrj_traj_LAST_RMSD"
        
        # Create the output file to get the final pdb.
        file_cpptraj_out  = "Pinptrj_traj_LAST_RMSD.out"
        
        # Open th einput file.
        file_cpptraj      = open (file_cpptraj_inp,"w")
        
        # Define the parm that has No water and the BS.
        parm_W = "[NoWat_BS]"
        
        # Define the variable that allow us to add the jump line at the end of each line.
        cc = "\n"

        # Create line to execute cpptraj and save the information in the output file and execute until EOF.
        line = " " + cpptraj + " >  " + file_cpptraj_out + " << EOF " + "\n"
        # Write line into the input file.
        file_cpptraj.write(line)
        
        # Add the line that allow cpptraj to reac the topology.
        line = " parm NoWat_BS."+ top_name + " " + parm_W + "\n"
        # Write line into the input file.
        file_cpptraj.write(line)
        
        # Add an empty line.
        line = " " + "\n"
        # Write line into the input file.
        file_cpptraj.write(line)

        # Add the line that allow us to read the input trajectory and get only the last snap of the trajectory (using trajin of cpptraj).
        line_inp = " trajin lig_" + lig_num + "_Superpose_BS.nc " + str(num_snaps) + " " + str(num_snaps) + " 1 parm [NoWat_BS] "  + cc
        # Write line into the input file.
        file_cpptraj.write(line_inp)

        # Save the Last snapshot into a pdb using trajout of cpptraj.
        line_out = " trajout lig_" + lig_num + "_traj_LAST.pdb pdb " + cc
        # Write line into the input file.
        file_cpptraj.write(line_out)

        # Tell to cpptraj what has to run.
        line_run = " run " + cc + cc
        # Write line into the input file.
        file_cpptraj.write(line_run)
        
        # Tell to cpptraj what is the ned of the file.
        line = "EOF"
        # Write line into the input file.
        file_cpptraj.write(line)

        # Close the input file.
        file_cpptraj.close()

        # Define the variable that has the execution bash command.       
        run_cpptraj = "./" + file_cpptraj_inp
        
        # Change the mode of the input file to be able to run the input file using bash.
        os.system('chmod u+x {} '.format(run_cpptraj) )
        
        # Run the input file using bash.
        os.system(run_cpptraj)

        # Print which file we are processing.
        print("Processing file : lig_" + lig_num + "_Superpose_BS_RMSD_ligand_alone.dat")

        # Create the input file to calculate the rmsd of each ligand alone in respect to the his BS.
        file_cpptraj_inp  = "Oinptrj_RMSD_ligand_alone"
        
        # Create the output file.
        file_cpptraj_out  = "Oinptrj_RMSD_ligand_alone.out"
        
        # Open the input file.
        file_cpptraj      = open (file_cpptraj_inp,"w")
        
        # Define the parm that has No water and the BS.
        parm_W = "[NoWat_BS]"
        
        # Define the variable that allow us to add the jump line at the end of each line.
        cc = "\n"
        
        # Create line to execute cpptraj and save the information in the output file and execute until EOF.
        line = " " + cpptraj + " >  " + file_cpptraj_out + " << EOF " + "\n"
        # Write line into the input file.
        file_cpptraj.write(line)
        
        # Add the line that tells to cpptraj which is the topology needed.
        line = " parm "+ "NoWat_BS." + top_name + " " + parm_W + "\n"
        # Write line into the input file.
        file_cpptraj.write(line)
        
        # Add an empty line.
        line = " " + "\n"
        # Write line into the input file.
        file_cpptraj.write(line)
        
        # Define the input name that contains the data of all the trajectory superposed.
        name_inp = "lig_" + lig_num + "_Superpose_BS.nc"

        # Define the input line that allows to go over all the trajectory.
        line_inp = " trajin " + name_inp + " 1 last 1 parm " + parm_W + " " + cc
        # Write line into the input file.
        file_cpptraj.write(line_inp)

        # Define the input line that tells to cpptraj which pdb has to be used as a reference (the last pdb of the trajectory (extracted before)).
        line_ref = " reference lig_" + lig_num + "_traj_LAST.pdb " + cc
        # Write line into the input file.
        file_cpptraj.write(line_ref)
        
        # Define which residue is from the ligand.
        num_res_ligand = num_res_prot + 1
        
        # Define the line that allows to cpptraj calculate the RMSD of the ligand along the trajectory in comparision of his final binding site. 
        line_out = " rms reference out lig_" + lig_num + "_Superpose_BS_RMSD_ligand_alone.dat :" + str(num_res_ligand) + "&!@H= nofit " + cc
        # Write line into the input file.
        file_cpptraj.write(line_out)
 
        # Line that tells to cpptraj what has to run.
        line_run = " run " + cc + cc
        # Write line into the input file.
        file_cpptraj.write(line_run)
        
        # Line that tells to cpptraj that is the end.
        line = "EOF"
        # Write line into the input file.
        file_cpptraj.write(line)

        # Close the input file.
        file_cpptraj.close()
        
        # Define the variable that has the execution bash command.
        run_cpptraj = "./" + file_cpptraj_inp

        # Change the mode of the input file to be able to run the input file using bash.
        os.system('chmod u+x {} '.format(run_cpptraj) )
        
        # Run the input file using bash.
        os.system(run_cpptraj)

        # Define in which file we have saved the data of the RMSD of the ligand against the Binding site.
        file_name = 'lig_' + lig_num + "_Superpose_BS_RMSD_ligand_alone.dat"
        
        # Call the function plot_RMSD_superposed_ligand to Create the plot of the RMSD along the trajectory.
        plot_RMSD_superposed_ligand(lig_num, file_name, snaps_by_one)
        
        # Move the picture saved in the initial directory to the pictures_RMSD_ligand folder.
        os.system("mv lig_" + str(lig_num) + "_RMSD_superposed_ligand.png pictures_RMSD_ligand/")
        
        # Remove the coordinates file created for each ligand.
        os.system("rm lig_" + str(lig_num) + "_Superpose_BS.nc")
    
    # Print that step 6 has ended.
    print ( ">>> END Processing Superposing BS and calculating RMSD of ligands alone <<<" )
  
  # Change the directory to the initial one.
  os.chdir(initial_dir)
  
  # Run the deletion of the non-Reactive trajectories.
  os.system(run_delete_nom)
 
  # Read final arguments
  en_R3 = args.en_R3
  queue_system = args.queue_system
  
  # If queue system is qsub 
  if (queue_system == "qsub*"):
    # Read the final argument
    cpu = args.cpu
    
    # Print the command that is going to be executed.
    print(queue_system + " -q " + cpu + " " + str(en_R3))
    
    # Execute next file R2
    os.system(queue_system + " -q " + cpu + " " + str(en_R3))
  # If queue system is not qsub (is sbatch).
  else:
    # Print the command that is going to be executed.    
    print(queue_system + " " + str(en_R3))
    
    # Execute next file R2
    os.system(queue_system + " " +  str(en_R3))

  # Return nothing of the main function.
  
  return

# From one list of data convert to the average of ns (using snaps_by_one(ns))
def snapshots_to_ns(Total_Energy_GB, snaps_by_one):
  # Define the number of snapshots we have from the total length of the list of snapshots.
  num_snapshots = len(Total_Energy_GB)
  
  # Define the number of ns depending on the number of snapshots.
  num_ns = int(num_snapshots/snaps_by_one)
  
  # Create on list to save the ns data.
  llista_ns = list()
  
  # Iterate the number of ns we need to get.
  for ns in range(0, num_ns):
    # Define a temporary list (to convert from snapshots to ns).
    llista_tmp = list()
    
    # Iterate all the snapshots that produce one ns.
    for snap in range(0, snaps_by_one):
      # Add those values of the list of snapshots into the temporary list.
      llista_tmp.append(Total_Energy_GB[0])
      
      # Remove the first element of the list (the one added to the temporary list).
      Total_Energy_GB.pop(0)

    # Once we have all the snapshots in the temporary list, make the average of that list (sum and divide).
    number = sum(llista_tmp)/len(llista_tmp)
    
    # Add that number into the list of ns.
    llista_ns.append(number)
    
  # Return the ns list.
  return llista_ns


# Function that Create a plot for each ligand.
def plot_RMSD_superposed_ligand(lig_num, file_name, snaps_by_one):
  # Define a boolean True to avoid the header of the file that needs to be readed.
  first = True
  
  # Create a list of snapshots data.
  snapshots_data = list()
  
  # Open the file that has all the values of each snapshots of the RMSD.
  with open(file_name, "r") as f:
    # For each line in the file
    for line in f:
      # If first is true
      if first:
        # first is false and continue with the next line.
        first = False
      else: # If first is false.
        # Split  the line and save into list t. 
        t = line.split()
        
        # Add the value of RMSD in the snpahsots_data list.
        snapshots_data.append(float(t[1]))
 
  # Convert all the snapshots to ns and save the data as Y-axis.
  y_data = snapshots_to_ns(snapshots_data, snaps_by_one)  
  
  # Create an arange with the same length as y_data (in increasing value 0,1,2,3,4...).
  x_data = np.arange(0, len(y_data), 1)
  
  # Transform the arange to a list.
  x_data = list(x_data) 
  
  ## Create a plot of the RMSD of the superposed ligand in respect to the BS.
  
  # Define the X-lim.
  plt.xlim(0, len(y_data))
  
  # Define the Y-lim
  plt.ylim(0, max(y_data))
  
  # Define the title of the plot.
  plt.title(str(lig_num) + " RMSD of superposed ligand")
  
  # Define the xlabel.
  plt.xlabel("t(ns)")
  
  # Define the ylabel.
  plt.ylabel("RMSD of ligand")
  
  # Define the data of the plot.
  plt.plot(x_data, y_data,'g', label='RMSD ligand alone')
  
  # Define the legend.
  plt.legend()
  
  # Save the plot.
  plt.savefig("lig_" + str(lig_num) + "_RMSD_superposed_ligand.png", format='png',dpi=600)
  
  # Close the plot.
  plt.close()
  
  # Return nothing for the funciton plot_RMSD_superposed_ligand.
  
  return

def plot_LIE_IntEnergy  (lig_num,file_data_LIE,num_ns,num_snaps,dir_fig_LIE ):
  file_LIE_energy   = 'lig_'+lig_num+'_LIE.csv'
  file_excel_LIE    = open (file_LIE_energy,"w")
  y_LIE = READ_dat_LIE ( file_data_LIE,num_snaps )
  x_LIE = np.linspace( 0, num_ns, num_snaps )

  for rowx in range ( num_snaps ) :
    strw = ( "{:15.9} {:15.9}".format(x_LIE[rowx],y_LIE[rowx]) )
    strw = str ( strw ) 
    strw  = strw + '\n'
    file_excel_LIE.write(strw)

  file_plot        = 'lig_'+lig_num
  file_plot_disat  = file_plot + '_LIE.png'
#     
  plt.xlim(0,num_ns)
  plt.ylim(0,Max_LIE_energy)

  plt.title(file_plot)
  plt.xlabel("t(ns)")
  plt.ylabel(" Interaction Energy (kcal/mol)")
# plt.plot (x_din,y_dis,'g', label='_nolegend_')
  plt.plot (x_LIE,y_LIE,'g', label='Interaction Energy')
  plt.legend()
  plt.savefig (file_plot_disat, format='png',dpi=600)
  plt.show()
  plt.close()

  move = 'mv ' + file_plot_disat + ' ./' + dir_fig_LIE + '/'
  os.system ( move  )
# 
  return

def plot_DisLigProtCenter (lig_num,file_data_DCM,num_ns,num_snaps,dir_fig_DCM ):
  file_DisLigProtCenter = 'lig_'+lig_num+'_DCM.csv'
  file_excel_DCM        = open (file_DisLigProtCenter,"w")
  y_DCM = READ_dat_DCM ( file_data_DCM,num_snaps )
  x_DCM = np.linspace( 0, num_ns, num_snaps )

  for rowx in range ( num_snaps ) :
    strw = ( "{:15.9} {:15.9}".format(x_DCM[rowx],y_DCM[rowx]) )
    strw = str ( strw ) 
    strw  = strw + '\n'
    file_excel_DCM.write(strw)

  file_plot        = 'lig_'+lig_num
  file_plot_disat  = file_plot + '_DCM' + str(distan_to_plot_DCM) + '.png'
#     
  plt.xlim(0,num_ns)
  plt.ylim(0,distan_to_plot_DCM)

  plt.title(file_plot)
  plt.xlabel("t(ns)")
  plt.ylabel("Distance (Angstrom)")
# plt.plot (x_din,y_dis,'g', label='_nolegend_')
  plt.plot (x_DCM,y_DCM,'g', label='Dist_PCM')
  plt.legend()
  plt.savefig (file_plot_disat, format='png',dpi=600)
  plt.show()
  plt.close()

  move = 'mv ' + file_plot_disat + ' ./' + dir_fig_DCM + '/'
  os.system ( move  )
# 
  return

def READ_dat_DCM (file_data_DCM,num_snaps) :
  y_DCM = []
  num_lines = num_snaps
  line_title = file_data_DCM.readline ().replace("\n","")
  for line in range(num_lines) :
    line_data   = file_data_DCM.readline ().replace("\n","")
    distance_sp = line_data.split()
    y_DCM.append(float(distance_sp[1]))
  y_DCM = np.array ( y_DCM ) 
  return y_DCM
  

def READ_dat_LIE (file_data_LIE,num_snaps) :
  y_LIE = []
  num_lines = num_snaps
  line_title = file_data_LIE.readline ().replace("\n","")
  for line in range(num_lines) :
    line_data   = file_data_LIE.readline ().replace("\n","")
    distance_sp = line_data.split()
    ele_energy  = float(distance_sp[1])
    vdW_energy  = float(distance_sp[2])
    tot_energy  = ele_energy + vdW_energy
    y_LIE.append(tot_energy)
  y_LIE = np.array ( y_LIE ) 
  return y_LIE

def calc_BindingSite (lig_num,num_at_prot,num_at_lig,dis_lig_prot_min,at_CA) :

  xm = ym = zm = 0 
  for at in range(num_at_lig):
    xm += x_lig[at]
    ym += y_lig[at]
    zm += z_lig[at]
  x_gc = xm / num_at_lig
  y_gc = ym / num_at_lig
  z_gc = zm / num_at_lig

  dis_intra_MAX = 0
  for atA in range(num_at_lig-1):
    xA = x_lig[atA]
    yA = y_lig[atA]
    zA = z_lig[atA]
    for atB in range(atA+1,num_at_lig):
      xB = x_lig[atB]
      yB = y_lig[atB]
      zB = z_lig[atB]
      disAB = sqrt( (xA-xB)*(xA-xB) + (yA-yB)*(yA-yB) + (zA-zB)*(zA-zB) )
      if disAB > dis_intra_MAX :
        dis_intra_MAX = disAB
  print ( "    ... Ligand GC   {:7.3f} {:7.3f} {:7.3f} Max Intra-Dist {:7.3f} ".format(x_gc,y_gc,z_gc,dis_intra_MAX)) 

  res_inBS = []
  dis_LIM = dis_intra_MAX + dis_lig_prot_min
  for atP in range(num_at_prot):
    xP = x_prot [atP]
    yP = y_prot [atP]
    zP = z_prot [atP]
    res_BS = res_num [atP]
    if res_BS in res_inBS :
      continue
    disPL = sqrt( (xP-x_gc)*(xP-x_gc) + (yP-y_gc)*(yP-y_gc) + (zP-z_gc)*(zP-z_gc) )
    if disPL < dis_LIM :      
      for atL in range(num_at_lig):
        xL = x_lig[atL]
        yL = y_lig[atL]
        zL = z_lig[atL]
        disPL = sqrt( (xP-xL)*(xP-xL) + (yP-yL)*(yP-yL) + (zP-zL)*(zP-zL) )
        if disPL < dis_lig_prot_min :
          res_inBS.append(res_BS)
          break  

  num_resBS = len(res_inBS)
  print ( "    ... There are  {:4} Protein Residues in this Binding Site ".format(num_resBS)) 

  name_fileBS = 'lig_' + lig_num + '_BS.dat'
  fileBS      = open (name_fileBS,"w")
  fileBS.write( str(num_resBS) + '\n' )

  llista = ""
  line_wrt = ''
  for res in range (num_resBS) :
    line_wrt += str(res_inBS[res]) + ' '
    llista += (str(res_inBS[res])) + ","
    
  print("lig_" + str(lig_num))
  print(llista[:-1])
  lig_BS["lig_" + str(lig_num)] = llista[:-1].split(",")
  print("ARA")
  fileBS.write( line_wrt + '\n' )
  fileBS.close()

  print ( "    ... RESIDUES : " )
  print ( line_wrt ) 
  print ( "... ... " )

  for at in range (num_resBS) :
    res_pos = res_inBS [at]
    CA_pos  = at_CA    [res_pos-1]
    x_CA    = x_prot   [CA_pos-1]
    y_CA    = y_prot   [CA_pos-1]
    z_CA    = z_prot   [CA_pos-1]
    print ( "    ... CA({}) :  {:7.3f} {:7.3f} {:7.3f} ".format(res_pos,x_CA,y_CA,z_CA)) 

  return 

""" GET_sasa """

def GET_sasa (name_file_SASA) :
  file_sasa  = open (name_file_SASA,"r")
  line_title = file_sasa.readline ().replace("\n","")
  line_data  = file_sasa.readline ().replace("\n","")
  surf = line_data.split()
  surf_complex  = float(surf[1])
  surf_receptor = float(surf[2])
  surf_ligand   = float(surf[3])
  sasa = (surf_complex - surf_receptor + surf_ligand ) / 2.0
  pc_sasa = sasa / surf_ligand * 100.0
  return sasa,pc_sasa

""" read_pocket """

def read_pocket (data_inp):
  residues = []
  file_pocket =  open (data_inp)
  num_lines = 0
  for lines in file_pocket :
    num_lines += 1
  if num_lines != 1 :
    print ( " Only ONE line is allowed : STOP ")
    exit ()
  file_pocket =  open (data_inp)
  line = file_pocket.readline ().replace("\n","")
  num_char = len (line)
# print ( " Num_Char = ",num_char )
  double_point = line.find(":")
# print ( "  :       = ",double_point )
  if double_point == -1 :
    print("Please use : ")
    exit ()
  line_nop = line[double_point+1:num_char]
# print( line_nop )
  coma = line.find(",")
  if coma == -1 :
    hairline = line.find("-")
    if hairline == -1 :
      residues.append ( int ( line_nop ) )
    else:
      range_res = line_nop.split("-")
      for i in range(len(range_res)):
        range_res[i] = int ( range_res[i] )
#     print ( range_res )
      for i in range(range_res[0],range_res[1]+1):
        residues.append ( i )
    print ( residues )
  else:
    all_range_res = line_nop.split(",")
#   print ( all_range_res )
    for i in range ( len(all_range_res) ) :
      if len (all_range_res[i]) == 0 :
        print ( " Error in Residue definition ")
        exit ()
    for i in range ( len(all_range_res) ) :
      range_res = all_range_res[i].split("-")
      for i in range(len(range_res)):
        range_res[i] = int ( range_res[i] )
      if len (range_res) == 2 :
        for j in range(range_res[0],range_res[1]+1):
          residues.append ( j )
      else:
          residues.append ( range_res[0] )

# print ( residues )
  return residues

""" Geometry Center for the pocket """

def calc_cent_geom ( res_pocket,num_at_prot,x_prot,y_prot,z_prot,at_name,res_num) :
  num_res_pocket = len(res_pocket)
# print ( "**",num_res_pocket)
  x_ref = 0.0
  y_ref = 0.0
  z_ref = 0.0
  for i in range(num_res_pocket):
    for n_at in range (num_at_prot):
      if  res_num[n_at] == res_pocket[i] and at_name[n_at] == 'CA' :
#       print (  res_num[n_at],res_pocket[i],at_name[n_at],x_prot[n_at] )
        x_ref += x_prot[n_at]
        y_ref += y_prot[n_at]
        z_ref += z_prot[n_at]
  x_ref /= num_res_pocket
  y_ref /= num_res_pocket
  z_ref /= num_res_pocket
  return x_ref,y_ref,z_ref

""" stat_x_ns  """     

def stat_x_ns (num_at,num_snaps_NOanal,y_dis ) :
  sy    =  0.0
  sd_y  =  0.0
  y_max =  0.0
  num_snaps_work =  num_at - num_snaps_NOanal
  for ns2anal in range ( num_snaps_NOanal,num_at ) :
    y   = y_dis [ns2anal] 
    sy  = y + sy
#   print ( ns2anal, y) 
    if y > y_max :
      y_max = y
  y_ave = sy / num_snaps_work
  for ns2anal in range ( num_snaps_NOanal,num_at ) :
    y  = y_dis[ns2anal] 
    sd_y = sd_y + ( y - y_ave ) * ( y - y_ave )
  sd_y = sqrt ( sd_y / num_snaps_work )
# print ( " y_ave ({:.3f} ),sd( {:.3f} ) y_max( {:.3f} )" .format(y_ave, sd_y, y_max )) 
  return y_ave, sd_y, y_max

""" dis_min_prot """

def dis_min_prot (num_at_prot,dis_lig_prot_min,x_cal,y_cal,z_cal):
  d_min = 100 * dis_lig_prot_min  
  num_nearest = 0
  for i in range( num_at_prot ) :
    x_tm = x_prot [i]
    y_tm = y_prot [i]
    z_tm = z_prot [i]
    x2 = ( x_tm - x_cal ) * ( x_tm - x_cal )
    y2 = ( y_tm - y_cal ) * ( y_tm - y_cal )
    z2 = ( z_tm - z_cal ) * ( z_tm - z_cal )
    d = sqrt ( x2 + y2 + z2 )
    if d < d_min :
      d_min = d
    if d <= dis_lig_prot_min :
      num_nearest += 1
# print ( "    ... D_min = {:.4f} ".format(d_min) )
  return d_min,num_nearest

""" read_prot """

def read_prot (pdb,name_atomrep):
  ABC = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']
  abc = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z']
  protein = ['SER','THR','GLN','ASN','TYR','CYS','CYX','CYM','GLY',  \
             'ALA','VAL','LEU','ILE','MET','PRO','PHE','TRP',        \
             'GLU','GLH','ASP','ASH','LYS','ARG','HIE','HID','HIP',  \
             'PHE','TYR','TRP' ]
  solvent = ['WAT' ]
  ions    = ['Cl-','Na+' ]

# print ( '    ... Reading     : {}  '.format(pdb) ) 
  file_pdb  = open (pdb)
  num_lines = 0
  for lines in file_pdb :
    num_lines += 1
# print ( "    ... There are = ",num_lines," Lines " )
  file_pdb  = open (pdb)
  num_atoms     = 0
  num_atoms_lig = 0
  resi_name_lig   = 'UNK'
  at_CA    = []
  for i in range ( num_lines )  :

    line = file_pdb.readline ().replace("\n","") 
    if 'ATOM' in line  :
      values=line.split()
      res_name_dum   = values [3]
      if res_name_dum in protein :
        num_atoms += 1
        is_CA  = values [2]
        if is_CA == 'CA' :
          at_CA.append(num_atoms)
        at_name.append (is_CA)
        res_name.append (values [3])
        chain_prot = values [4]
        if (chain_prot in ABC) or (chain_prot in abc) :
          res_num.append (int (values [5]))
          x_prot.append (float(values [6]))
          y_prot.append (float(values [7]))
          z_prot.append (float(values [8]))
        else:
          res_num.append (int (values [4]))
          x_prot.append (float(values [5]))
          y_prot.append (float(values [6]))
          z_prot.append (float(values [7]))
      else :
        if res_name_dum not in ions :
          values=line.split()
#         print ( values )
          atom_name_lig.append (values [2])
          resi_name_lig   = values [3]
          chain_lig = values [4]
          if (chain_lig in ABC) or (chain_lig in abc) :
            x_lig.append (float(values [6]))
            y_lig.append (float(values [7]))
            z_lig.append (float(values [8]))
          else :
            x_lig.append (float(values [5]))
            y_lig.append (float(values [6]))
            z_lig.append (float(values [7]))
          num_atoms_lig += 1

  posi_name_atomrep = 0
  if  num_atoms_lig != 0 :
    for j in range ( num_atoms_lig )  :
      if atom_name_lig[j] == name_atomrep :      
        posi_name_atomrep = j

  """
  print ( "... ... There are = ",num_atoms," Protein Atoms " )
# for i in range ( num_atoms )  :
#   print ( " ATOM {:.3f} {:.3f} {:.3f} ".format(x_prot[i],y_prot[i],z_prot[i])) 
  if num_atoms_lig != 0 :
    print ( "... ... There are = ",num_atoms_lig," Ligand Atoms " )
    for i in range ( num_atoms_lig )  :
      print ( " ATOM {:.3f} {:.3f} {:.3f} ".format(x_lig[i],y_lig[i],z_lig[i])) 

# exit ()
  """
  return num_atoms,num_atoms_lig,posi_name_atomrep,resi_name_lig,at_CA

""" read_pdb_c99 """

def read_pdb_c99 (pdb,name_atomrep,x,y,z):
  file_pdb  = open (pdb)
  num_lines = 0
  for lines in file_pdb :
    num_lines += 1
# print ( "... ... There are = ",num_lines," Lines " )
  file_pdb  = open (pdb)
  num_atoms = 0
  for i in range ( num_lines )  :
    line = file_pdb.readline ().replace("\n","") 
    if name_atomrep in line :
      x_chg = ''
      y_chg = ''
      z_chg = ''
      for i in range(30,38) :
        x_chg += line[i]
      for i in range(38,46) :
        y_chg += line[i]
      for i in range(46,54) :
        z_chg += line[i]
      x.append (float(x_chg))
      y.append (float(y_chg))
      z.append (float(z_chg))
      num_atoms += 1
  print ( "    ... There are = ",num_atoms," Atoms(Structures) " )
  return num_atoms

""" END read_pdb_c99 """

# ... READ datafiles

def cmdlineparse():

  parser = ArgumentParser(description="command line arguments")

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
  
  parser.add_argument("-top_name"     , dest="top_name"     , required=True, help=" Name of the Topological file ")
  parser.add_argument("-num_res_prot" , dest="num_res_prot" , required=False, help=" Number of Protein Residues ")
  parser.add_argument("-fragment" , dest="fragment" , required=False, help=" Name of the fragment used ")
  parser.add_argument("-dir_traj" , dest="dir_traj" , required=False, help=" Directory where is located the trajectory ")
  parser.add_argument("-end_traj" , dest="end_traj" , required=False, help=" Number of final trajectory ")
  parser.add_argument("-en_R3"    ,    dest="en_R3"          , required=False, help = " Name of the next file to execute en_R3 ")
  parser.add_argument("-queue_system"    ,    dest="queue_system"          , required=False, help = "Command to run the script  with bash")
  parser.add_argument("-cpu"    ,    dest="cpu"          , required=False, help = " CPu we want to use to execute the files ")


  args=parser.parse_args()
  return args

def prog_info () :
  t = datetime.now()
  format_date = "%d-%m-%Y %H:%M:%S"
  ti= t.strftime(format_date)
  print ( "... ... ... ... ... ... ... " )
  print ( "...  fdMD_ReactiveTraj  ... " )
  print ( "...     ( 2020 )        ... " )
  print ( "... ... ... ... ... ... ... " )
  print ( " ",ti )
  print ( "... ... ... ... ... ... " )
  print ( " " )

if __name__ == '__main__':

  distan_to_plot_DCM =   100 
  Max_LIE_energy     =  -100.0
  
  amber_version = 20
  amber_dir     =  "/programas/amber20/amber" + str(amber_version)
  cpptraj       =  amber_dir + "/bin/cpptraj"


  main()

