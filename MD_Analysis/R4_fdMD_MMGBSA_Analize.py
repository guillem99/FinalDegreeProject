#/usr/bin/env python3

"""
 .......................
  fdMD_MMGBSA_Analize.py
 .......................

The structure of directories MUST be :
1) Dir where you run this file
- prefix_Dir_1
- prefix_Dir_2
- .....
Into these directories you MUST have
---- ffdir
---- fftop
"""
#######################################################################
#                    DATA to Modify:
#######################################################################
# Sanpshots used to calculate the DeltaG.
snaps_by_one = 50

# Snapshots used to calculate the RMSD and the LIE Energy.
snaps_by_one_RMSD = 50

# ns we want to analyse and calculate the Average.
ns_ave       = 100

# Prefix of the ligands and files of each ligand.
prefix       ='lig_'

# Time we want to average and sum at each step to calculate the residence time. 
time_residence_avg = 10

# Boolean to decide if we want to analyse.
anal_error   = True

# Maximum distance we want to analyse.
max_posi     = 5.0

#######################################################################
#                     END Data to modify
#######################################################################

# Import libraries needed.
import numpy as np
import os
import time
from math import sqrt
import  matplotlib
matplotlib.use('Agg')
import  matplotlib.pyplot as plt
##from argparse import ArgumentParser

# Define Strings
gb_mmpbsa = "GENERALIZED BORN:"
pb_mmpbsa = "POISSON BOLTZMANN:"
delta     = "DELTA Energy Terms"

# Define booleans to decide if you want to calculate the LIE Energy avg.
LIE_Energy = True

# Define Boolean to decide if you want to calculate the Hydrogen Bonds.
Hydrogen_Bonds = True

# Define Boolean to decide if you want to calculate the SASA Energy.
SASA_avg_calculation = True


def solve_error ( val_ener ) :
  num_val = val_ener.size
  val_ener_ord = val_ener.copy()
  ind_ener_ord = np.argsort(val_ener_ord)[::-1]
  val_ener_ord = np.sort   (val_ener_ord)[::-1]

  zero = 0.0
  ind_count = 0
  value = val_ener_ord [ind_count]
  ind_zero = -1
  while value > zero and ind_count < num_val :
    ind_zero  += 1 
    ind_count += 1
    value = val_ener_ord [ind_count]
  if ind_zero >= 0 :
    for ind in range (ind_zero) :
      if val_ener_ord [ind] >= max_posi :
        val_ener [ind_ener_ord[ind]] = 0.1* np.random.random(1) - 0.1
        print ( ind, val_ener[ind_ener_ord[ind]] )
  return val_ener

def ave_data (energy) :
  is_div = num_ns%ns_ave
  if is_div == 0 :
    num_interval = int ( num_ns // ns_ave )
    print ( "num_interval",num_interval)
  else :
    num_interval = 1

  print ('num_ns,ns_ave,num_interval',num_ns,ns_ave,num_interval)
  ave_interval = []
  num_ns_interval = int ( num_ns / num_interval )
  data_interval    = num_ns_interval * snaps_by_one
  ind = 0
  for interval in range ( num_interval ) :
#   print ( " >>> interval = ",interval )
    ave = 0.0
    for data in range ( data_interval ) :
#     print ( " >>> >>> data = ",data )
      ave = ave + energy [ind]
#     print ( " >>> >>> >>> ind = ",ind )
      ind += 1
    ave /= data_interval
    ave_interval.append ( ave ) 
  print ( ave_interval )

  ave_graph = []
  for interval in range ( num_interval ) :
    for data in range ( data_interval ) :
      ave_graph.append ( ave_interval [interval] )

  num_data_new = len ( ave_graph )
  print ( "num_data_new = ",num_data_new )

  return ave_interval, ave_graph


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

# Function to calculate the residence time given a list of ns.
def tiempo_de_residencia(Total_Energy_GB):
  # Total_Energy_GB = average value of GB for each ns
  # Convert to a list the list of ns. 
  Total_Energy_GB     = list(Total_Energy_GB)

  # Save the original list into another list.
  Total_Energy_GB_ini = Total_Energy_GB[:]

  # Save the total_time_residence to the time we need to average at each step.
  total_time_residence = time_residence_avg
  
  # Create a list to calculate the average of all.
  llista_avg_all = list()

  # As we need the residence time of the last period, we do the reverse of the list (Start from the end).
  Total_Energy_GB.reverse()

  # Define the number of intervals we have during all the ns.
  num_intervals = int(len(Total_Energy_GB)/time_residence_avg)

  # Print the Time avg we are using.
  print ( ' ... Time_avg      = ',time_residence_avg)

  # Print the number of intervals that will be checked.
  print ( ' ... Num Intervals = ',num_intervals)

  # Create a list of the Last interval 
  llista_avg1 = list()

  # Iterate all the values of the residence time we want to check (f.ex. 0 -> 10).
  for i in range(0, time_residence_avg):
    # Get the DeltaG of the ns iterated in the list of the avg.
    llista_avg1.append(Total_Energy_GB[0])

    # Get the DeltaG of the ns iterated in the list of the avg during all the trajectory.
    llista_avg_all.append(Total_Energy_GB[0])
  
    # Remove that DeltaG of the ns iterated.
    Total_Energy_GB.pop(0)

  # Create the avg of the list of the last interval.
  avg0 = sum(llista_avg1)/len(llista_avg1)

  # Save the result avg during the last interval as result_avg_100.
  result_avg_100 = avg0
  
  # Print the result avg of the last interval.
  print ( ' ... ... Average DeltaG last interval = {:6.3f} '.format(avg0) )
  
  # Calculate sigma of the last interval.
  suma = 0

  # Iterate the list of the Last interval.
  for ener in llista_avg1:
    # We sum all the difference values between the value and the avg of the list and to the power of 2.
    suma += (ener-avg0)**2
  
  # We calculate the avg of the sum done before.
  sigma2  = suma/(len(llista_avg1))

  # From the sigma2 value, we apply the square root and obtain the value of sigma.
  sigma   = np.sqrt(sigma2)

  # We define a value of three sigma to test the residence time.
  THREEsigma= 3 * sigma

  # We define the min value found with 3 sigma.
  min_value  = avg0 - THREEsigma

  # We define the max value found with 3 sigma.
  max_value  = avg0 + THREEsigma

  # Print the Sigma value used of the last interval.
  print ( ' ... ... Sigma   value  last interval = {:6.3f}'.format(sigma) )

  # Print the value used to test all inervals (3*sigma).
  print ( ' ... ... THREE  *  Sigma  last interval = {:6.3f}'.format(THREEsigma) )

  # Print the minimum and maximum values we found with 3*sigma and the avg of the last residence time interval.
  print ( ' ... ... [ Min Max ] interval = [{:6.3f},{:6.3f}]'.format(min_value,max_value) )

  # Iterate all the intervals found before (not including the first (last) one).
  for t in range(1, num_intervals ) :
    # Define a list to calculate the avg of each interval.
    llista_avg1 = list()
    
    # Iterate the ns we want to analyse in each step (in each interval).
    for i in range(0, time_residence_avg):
      # Add the value of DeltaG in the avg list of the interval.
      llista_avg1.append(Total_Energy_GB[0])

      # Add the value of DeltaG in the avg list of all the trajectory.
      llista_avg_all.append(Total_Energy_GB[0])

      # Remove the value of Delta G (ns) of the list that includes all the non analysed ns. 
      Total_Energy_GB.pop(0)

    # Calculate the avg of the interval we are iteraing.
    avg1 = sum(llista_avg1)/len(llista_avg1)
    
    # If the avg of the interval is between the minimum and the maximum value.
    if min_value <= avg1 <= max_value:
      # Add the residence time (ns of the interval) to the total residence time.
      total_time_residence += time_residence_avg
    else: # If not
      # Calculate the minimum value of the total residence time and the ns wanted to average.
      lenght_ave = min(total_time_residence, ns_ave)

      # Get the list of values to calculate the avg of the length wanted (length_ave).
      llista_avg_last = Total_Energy_GB_ini[-lenght_ave:]
      
      # Calculate the avg of the length wanted (length_ave).
      result_avg_100 = sum(llista_avg_last)/len(llista_avg_last)

      # Break the loop and return
      break

  # Return the total residence time and the average in that total residence time.
  return total_time_residence, result_avg_100

# Define the function to read the RMSD file.
def read_RMSD_Superposed_ligand(file_name): # The only argument makes reference to the RMSD file.
  # Define a first boolean to True (to avoid the first line of the file (header)).
  first = True

  # Define the list of the Snapshots we are going to read in the RMSD file.
  snapshots_data = list()

  # Open the RMSD file.
  with open(file_name, "r") as f:
    # For each line.
    for line in f:
      # If is the first line.
      if first:
        # Don't do anything (because it is the header) (Put first to False).
        first = False
      # If is not the first line.
      else:
        # Split the data of the line and save to t list.
        t = line.split()
       
        # Save the RMSD value we want to to the list of snapshots defined before.
        snapshots_data.append(float(t[1]))

  # Return the list of snapshots with RMSD values.
  return snapshots_data


initial_dir = os.getcwd()
print ( " Working in the Directory : ",initial_dir )

exist_dir = os.path.isdir('pictures')
if exist_dir : 
  if os.listdir('pictures') != [] :
    os.chdir('pictures')
    os.system ( "rm * " )
    os.chdir('..')       
  os.rmdir('pictures')
os.mkdir('pictures')

dir_all = os.listdir(os.getcwd())
num_dir = 0
for dir in dir_all:
  if os.path.isdir(dir):
    if prefix in dir :
      num_dir += 1
      
      print(" {} {}".format(num_dir,dir))

print ( " There are {} Directories : ".format(num_dir)) 

file_ave  = open ( "ave_gbpb.res","w" )

for dir in dir_all:
  if os.path.isdir(dir):
    if prefix in dir :
      name_dir = dir.split("_")
      molec =  name_dir  [0]
      pose  =  name_dir  [1]
      line_to_write = str (molec) + '_' + str (pose) 
      print ( molec,pose )
      wd = os.getcwd()
      os.chdir(dir)
      print ( " >>> Working in the Directory : ",dir )
      NOTmmgbsa = True  
      dir_mmgbsa_all = os.listdir(os.getcwd())
      for dir_mmgbsa in dir_mmgbsa_all:
        if os.path.isdir(dir_mmgbsa):
          if 'mmpbsa' in dir_mmgbsa :
            NOTmmgbsa = False 
            os.chdir(dir_mmgbsa)
            print ( " >>> >>> YES mmpbsa DIR ")
            break
      if NOTmmgbsa :
        print ( " >>> >>> NO mmpbsa DIR ")
        os.chdir(initial_dir)
        continue
#     dir_mmpbsa = ("mmpbsa_{}".format(dir))
#     print ( " >>> >>> mmpbsa_ Directory : ",dir_mmpbsa )
#     os.chdir(dir_mmpbsa)
      NOffdir = True  
      dir_ffdir_all = os.listdir(os.getcwd())
      for dir_ffdir in dir_ffdir_all:
        if os.path.isdir(dir_ffdir):
          if 'ffdir' in dir_ffdir :
            NOffdir = False
            os.chdir(dir_ffdir)
            print ( " >>> >>> >>> ffdir_ Directory : ",dir_ffdir )
            Files_ffdir = os.listdir(os.getcwd())
            if 'DeltaG_BySnaps.csv' not in Files_ffdir :
              print ( " >>> >>> >>> MMGBSA not done ")
              os.chdir("../..")
              break

            os.system("cp DeltaG_BySnaps.csv ./tempfile.csv")
#           os.system("ls -l")
#           
            
# 1. GENERALIZED BORN
#
            file_mmpbsa = open ( "tempfile.csv","r+" )
            gb_data  = False
            gb_delta = False
            lin_gb = 0
            for line in file_mmpbsa :
              lin_gb += 1
              if gb_mmpbsa in line :
                gb_data = True
                gp_pos  = lin_gb
                print ( " >>> >>> >>> >>> GENERALIZED BORN " )
              if gb_data and not gb_delta:
                if delta in line :
                  gb_delta = True
                  gb_delta_pos  = lin_gb
                  print ( " >>> >>> >>> >>> >>> GB Delta ",gb_delta_pos )
                  break
#
# 1. POISSON BOLTZMANN 
#
            file_mmpbsa = open ( "tempfile.csv","r+" )
            pb_data  = False
            pb_delta = False
            lin_pb = 0
            for line in file_mmpbsa :
              lin_pb += 1
              if pb_mmpbsa in line :
                pb_data = True
                pp_pos  = lin_pb
                print ( " >>> >>> >>> >>> POISSON BOLTZMANN " )
              if pb_data and not pb_delta:
                if delta in line :
                  pb_delta = True
                  pb_delta_pos  = lin_pb
                  print ( " >>> >>> >>> >>> >>> GB Delta ",pb_delta_pos )
                  break
#
# 2. GENERALIZED BORN
#
            file_mmpbsa = open ( "tempfile.csv","r+" )
            if gb_delta :
              Total_Energy_GB = []
              for lin_null in range ( 1,gb_delta_pos+1 ) :
                line = file_mmpbsa.readline ().replace("\n","")
#               print ( " line ",line )
              line_info = file_mmpbsa.readline ().replace("\n","")
              names_ener =  line_info.split(",")
              num_term = len (names_ener)
              for num_t in range(num_term):
#               print ( names_ener[num_t])
                if  names_ener[num_t] == "DELTA TOTAL" :
                  ind_ETot_GB = num_t
#             print ( " Total Energy = ", names_ener[ind_ETot_GB])
              num_ener = 1
              line_ener = file_mmpbsa.readline ().replace("\n","")
              while line_ener != "" :
                energies = line_ener.split(",")
                Total_Energy_GB.append( float (energies [ind_ETot_GB] ))
                line_ener = file_mmpbsa.readline ().replace("\n","")
                num_ener += 1
              num_ener -= 1
              print ( " There are {} GB snapshots  ".format(num_ener) )
              num_ns = num_ener / snaps_by_one
              
              y_gb = np.array ( Total_Energy_GB )
              if anal_error :
                y_gb = solve_error ( y_gb )
              
              # Convert Snapshots to ns of DeltaG (Energy of Binding)
              llista_ns = snapshots_to_ns(y_gb, snaps_by_one)

              # Call the residence time funciton to calculate the residence time with DeltaG
              time, avg_100 = tiempo_de_residencia(llista_ns)

              # Print the residence time.
              print("Tiempo de residencia = ", time)

              # Print the Delta G Average.
              print("Average Delta G of last 100 ns = ", avg_100)

              # Get the file name that has the RMSD of the ligand we are iterating.
              file_name = wd + "/" + line_to_write[0:7] + "_Superpose_BS_RMSD_ligand_alone.dat"

              # Read the RMSD file for the ligand we are iterating.
              llista_snapshots_RMSD = read_RMSD_Superposed_ligand(file_name)             
              
              # Convert from Snapshots to ns the RMSD and calculate the Residence time with RMSD
              llista_ns_RMSD = snapshots_to_ns(llista_snapshots_RMSD, snaps_by_one_RMSD)

              # Calculate the residence time using the RMSD instead of using DeltaG values.
              time_RMSD, avg_100_RMSD = tiempo_de_residencia(llista_ns_RMSD)

              # Print the Residence time using RMSD.
              print("Tiempo de residencia RMSD = ", time_RMSD)

              # Print the RMSD Average.
              print("Average RMSD of last 100 ns = ", avg_100_RMSD)
             
              ave_gb, y_ave_gb = ave_data (y_gb)
              y_ave_gb = np.array ( y_ave_gb )
              line_to_write = line_to_write + ' GB '
              for i in range ( len(ave_gb) ) :
                line_to_write = line_to_write + str(ave_gb[i]) + ' '

#
# LIE Energy
#
            # If LIE_Energy is true.
            if LIE_Energy:
             # Define a list to catch the Evdw values.
             EVDW = list()

             # Define a list to catch the EElec values.
             EELEC = list()

             # Define a list to catch the LIE values.
             LIE = list()

             # Define a boolean to avoid the first line.
             first = True

             # Open the file that has the LIE Energy.
             with open(wd + "/" + line_to_write[0:7] + "_LIE.dat" ,"r") as f:
               # For each line in the LIE file.
               for line in f:
                 # If is the first line.
                 if first:
                   # Set first to False.
                   first = False
                 # If is not the first line.
                 else:
                   # SAve the data of the line in a list called t. 
                   t = line.split()
                   
                   # Add the value of Evdw in his list.
                   EVDW.append(t[2])

                   # Add the value of EElec in his list.
                   EELEC.append(t[1])

                   # Add the value of LIE in his list.
                   LIE.append(float(t[1]) + float(t[2]))

             # Convert the snapshots readed until now to ns (of the LIE Energy).
             LIE_ns = snapshots_to_ns(LIE, snaps_by_one_RMSD)

             # Catch the last ns we want to average.
             llista_LIE = LIE_ns[-ns_ave:]

             # Calculate the residence time using the list caught before.
             result_avg_LIE = sum(llista_LIE)/len(llista_LIE)

             # Print the average of LIE Energy of the last ns to analyse.
             print("Average of LIE ENergy of last 100 ns = ", result_avg_LIE)
             
#
# Hydrogen Bonds calculation 
#
            # If Hydrogen bonds is true.
            if Hydrogen_Bonds:
             # Define that we have 0 Hydrogen bonds.
             number_hyd_bonds = 0

             # Open the file that has all the Hydrogen Bonds.
             with open(wd + "/" + line_to_write[0:7] + "_avg-hb_fragment_last100.dat","r") as f:
               # For each line in the file.
               for line in f:
                 # Split the data of the line and save it into t list.
                 t = line.split()
 
                 # If the Hydrogen bond is found in more than the 50% of the last ns to anlyse and the average distance is less than 3 Amstrongs.
                 if (float(t[4]) > 0.50 and float(t[5]) <= 3.00):
                   # Add one hydrogen bonds to the total.
                   number_hyd_bonds += 1

             # Finally, print how many Hydrogen bonds we have.
             print("Number of hydrogen bonds = ", number_hyd_bonds)                
 
#
# SASA avg last 100 ns calculation
#
            # If SASA Calculation is true.
            if SASA_avg_calculation:
              # Create a list for the values of SASA of the Complex.
              llista_Complex = list()

              # Create a list for the values of SASA of the Receptor.
              llista_Receptor = list()

              # Create a list for the values of SASA of the Ligand.
              llista_Ligand = list()

              # Open the file that includes all the values of SASA.
              with open(wd + "/" + line_to_write[0:7] + "_last_100_SASA.out", "r") as f:
                # Read the line of the title and don't do anything with it.
                line_title = f.readline().replace("\n","")

                # Iterate the SASA from 0 to the total of ns to analyse we want (ns_ave).
                for p in range(0, ns_ave):
                  # Read the line p.
                  line_data = f.readline().replace("\n","")
                 
                  # Split the data of line p and save it into surf list.
                  surf = line_data.split()
 
                  # If the length of the line is greater than 2.
                  if len(surf) > 2:
                    # Get the value of the Complex surface of the ns iterated.
                    surf_complex  = float(surf[1])

                    # Get the value of the Receptor surface of the ns iterated.
                    surf_receptor = float(surf[2])

                    # Get the value of the Ligand surface of the ns iterated.
                    surf_ligand   = float(surf[3])

                    # If some surface is equal to -1, that means that has been wrongly calculated, so we continue.
                    if surf_complex == -1.0000 or surf_receptor == -1.0000 or surf_ligand == -1.0000:
                      # Continue with the next line.
                      continue
                    # If not.
                    else:
                      # Add the Complex SASA value to the list of the Complex.
                      llista_Complex.append(surf_complex)

                      # Add the Receptor SASA value to the list of the Receptor.
                      llista_Receptor.append(surf_receptor)

                      # Add the Ligand SASA value to the list of the Ligand.
                      llista_Ligand.append(surf_ligand)
                  # If the length is not greater than 2.
                  else:
                    # Print that there is an error reading the file.
                    print("Error while reading the " + line_to_write[0:7] + "_last_100_SASA.out file.")
                    
                    # Break the loop.
                    break

              # Calculate the average of the SASA Complex.
              surf_complex = sum(llista_Complex)/len(llista_Complex)

              # Calculate the average of the SASA Receptor.
              surf_receptor = sum(llista_Receptor)/len(llista_Receptor)
              
              # Calculate the average of the SASA Ligand.
              surf_ligand = sum(llista_Ligand)/len(llista_Ligand)

              # Calculate the SASA Energy.
              sasa = (surf_complex - surf_receptor + surf_ligand ) / 2.0

              # Calculate the percentage of SASA of the ligand. 
              pc_sasa = sasa / surf_ligand * 100.0

              # Print the SASA avg value.
              print("SASA avg last 100 ns = ", sasa)

              # Print the %SASA avg value.
              print("% SASA avg last 100 ns = ", pc_sasa)

  
#
# 2. POISSON BOLTZMANN
#
            file_mmpbsa = open ( "tempfile.csv","r+" )
            if pb_delta :
              Total_Energy_PB = []
              for lin_null in range ( 1,pb_delta_pos+1 ) :
                line = file_mmpbsa.readline ().replace("\n","")
#               print ( " line ",line )
              line_info = file_mmpbsa.readline ().replace("\n","")
              names_ener =  line_info.split(",")
              num_term = len (names_ener)
              for num_t in range(num_term):
#               print ( names_ener[num_t])
                if  names_ener[num_t] == "DELTA TOTAL" :
                  ind_ETot_PB = num_t
#             print ( " Total Energy = ", names_ener[ind_ETot_PB])
              num_ener = 1
              line_ener = file_mmpbsa.readline ().replace("\n","")
              while line_ener != "" :
                energies = line_ener.split(",")
                Total_Energy_PB.append( float (energies [ind_ETot_PB] ))
                line_ener = file_mmpbsa.readline ().replace("\n","")
                num_ener += 1
              num_ener -= 1
              print ( " There are {} PB snapshots  ".format(num_ener) )
              num_ns = num_ener / snaps_by_one
              y_pb = np.array ( Total_Energy_PB )
              if anal_error :
                y_pb = solve_error ( y_pb )
              ave_pb, y_ave_pb = ave_data (y_pb)
              y_ave_pb = np.array ( y_ave_pb )
              line_to_write = line_to_write + ' PB '
              for i in range ( len(ave_gb) ) :
                line_to_write = line_to_write + str(ave_pb[i]) + ' '

            line_to_write = line_to_write + '\n'
            file_ave.write(line_to_write)

            title  = str(molec) + '_' + str(pose) 
            file_plot = '../../../pictures/' + dir + '.png'

            label_gb  = 'GB'
            label_pb  = 'PB'

            style_gb = '--'
            style_pb = ':'

            x_din = np.linspace( 0, num_ns, num_ener )

            plt.ylim()
            plt.xlim(0,num_ns )
            plt.xlabel("ns")
            plt.ylabel(r'$\Delta$G (kcal/mol)')

            plt.title(title)

            if  gb_delta and pb_delta :
              plt.plot (x_din,y_gb,color='b', label=label_gb)
              plt.plot (x_din,y_pb,color='r', label=label_pb)
              plt.plot (x_din,y_ave_gb,color='y',linestyle=style_gb,label=label_gb)
              plt.plot (x_din,y_ave_pb,color='y',linestyle=style_pb,label=label_pb)
            elif     gb_delta and not pb_delta :
              plt.plot (x_din,y_gb,'b', label=label_gb)
              plt.plot (x_din,y_ave_gb,color='y',linestyle=style_gb,label=label_gb)
            elif not gb_delta and     pb_delta :
              plt.plot (x_din,y_pb,'r', label=label_pb)
              plt.plot (x_din,y_ave_pb,color='y',linestyle=style_pb,label=label_pb)

            plt.legend()
            plt.savefig (file_plot, format='png',dpi=600)
            plt.show()
            plt.close ()

            os.system("rm ./tempfile.csv")
#           os.system("ls -l")
            os.chdir("../..")

#     os.chdir("..")
      os.chdir("..")
 
exit ()

