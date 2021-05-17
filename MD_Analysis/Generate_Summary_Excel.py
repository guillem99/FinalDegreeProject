#!/usr/bin/env python3
import numpy as np
import os, sys, stat
from datetime import datetime, date, time
from math import sqrt
import  matplotlib
matplotlib.use('Agg')
import  matplotlib.pyplot as plt
from argparse import ArgumentParser
import openpyxl

'''
###########################################
###              USAGE                  ### 
###                                     ###
###########################################

 python3 Generate_Summary_Excel_+SASA.py -sim_time 1000 -molecule_receptor mpro -info_R2 Info_R2_fdMD_ReactiveTraj_26082020_modified  -ligand ebse -num_ns_anal 100 -info_R4 Info_R4_fdMD_MMGBSA_Analize_modified_25012021 -num_dyns 4

## You have to be in the 1st dynamic of the ligand you want to create the _all_dyns_summary
## for instance, /scr/guillemc/mpro_pd1/dyn_1/ 

'''

##################
##     MAIN     ##
##################

def main():

    print(" Using AMBER : ", amber_dir)
    # Get the directory we are working with.
    DIR = os.getcwd()
    
    # List the files we have in that directory.
    list_DIR = os.listdir(DIR)
    
    # Print the initial directory path.
    print(" Initial Dir : ", DIR)
    print("  ")

    
    #######################
    ### Parse arguments ###
    #######################
    # Catch all the arguments we have.
    args = cmdlineparse()
    
    # Catch each argument alone.
    ligand = args.ligand
    molecule_receptor = args.molecule_receptor
    Simulation_time = args.sim_time
    info_R2 = args.info_R2
    info_R4 = args.info_R4
    num_ns_anal = args.num_ns_anal
    num_dyns = int(args.num_dyns)
  
    ############################
    ###     Create Excel     ###
    ############################
    # Print a line of reminder.
    print("!!!  Remember I'm using Info_R2 and Info_R4 files, so, please, don't change the prints of those files  !!!")
    print("")
    
    # Print start line
    print(">>> Starting ")
    print("")
    
    # Print the first step.
    print("1. Create Excel File with headers of the table.")
    
    # Create a woorkbook using Excel library (openpyxl).
    wb = openpyxl.Workbook()
    
    # Create a sheet called Sheet in that Excel workbook.
    sheet = wb.create_sheet("Sheet")
    
    # Activate the workbook to be able to add data in the excel file.
    hoja = wb.active
    
    # Iterate the number of dyns we have.
    for i in range(0, num_dyns):
        # For each Molecular Dynamic crate the header (which include all the descriptors).
        hoja["A" + str(i*15 + 1)] = "dyn_" + str(i+1)
        hoja["A" + str(i*15 + 2)] = "Simulation time"
        hoja["B" + str(i*15 + 2)] = "Ligand"
        hoja["C" + str(i*15 + 2)] = "Molecule"
        hoja["D" + str(i*15 + 2)] = "Ligand number"
        hoja["E" + str(i*15 + 2)] = "avg LIE of last 100ns"
        hoja["F" + str(i*15 + 2)] = "Residence time with RMSD"
        hoja["G" + str(i*15 + 2)] = "Avg Delta G last 100ns"
        hoja["H" + str(i*15 + 2)] = "Residence time"
        hoja["I" + str(i*15 + 2)] = "Binding Site"
        hoja["J" + str(i*15 + 2)] = "Sasa LAST snapshot"
        hoja["K" + str(i*15 + 2)] = "% Sasa LAST Snapshot"
        hoja["L" + str(i*15 + 2)] = "Sasa ligand last 100ns"
        hoja["M" + str(i*15 + 2)] = "% Sasa ligand last 100ns"
        hoja["N" + str(i*15 + 2)] = "Number of hydrogen bonds"
        hoja["O" + str(i*15 + 2)] = "Protein Residues"
        
        # Add the simulation time, ligand and molecule_receptor for each MD.
        hoja["A"+ str(i*15 + 3)] = Simulation_time 
        hoja["B"+ str(i*15 + 3)] = ligand 
        hoja["C"+ str(i*15 + 3)] = molecule_receptor
    
    # Add the formula of the LIE Energy (as a reminder).
    hoja["E1"] = "LIE = EELEC + EVDW"
    
    ##############################################
    ###  Catch REACTIVE ligands and fill Excel ###
    ##############################################
    
    # Print the second step.
    print("2. Catch the REACTIVE ligands and fill the table")
    
    # Create the dictionary in which we are going to save the different binding sites.
    Binding_sites = dict()
    
    # Create a dictionary to save the Binding sites of each dyn.
    Binding_each_dyn = dict()
    
    # Iterate the number of dyns we have.
    for i in range(0, num_dyns):
        # Define the Directory of the different dyns.
        DIR = DIR[:-1] + str(i+1)
        
        # Create a list of the reactive ligands for each dyn.
        Reactive_ligands = list()
        
        # Create a get_in boolean to know where we have our data.
        get_in = False
        
        # Open the info_R2 file.
        with open(DIR + "/" + info_R2,"r") as f:
            # For each line in the file.
            for line in f:
                # Split the data of each line. 
                t = line.split()
               
                # if get_in is true 
                if get_in:
                    # get_in to false.
                    get_in = False
                    
                    # Save the Reactive trajectories in the list of Reactive ligands.
                    Reactive_ligands = t[1:] # Reactive_ligands = ["621", "632", etc]
                    
                # For each element in the line.    
                for ele in t:
                    # If the element is equal to "Active"
                    if ele == "Active":
                        # Get_in is true (to know the Reactive trajectories).
                        get_in = True

        # Save Get_in as False.
        get_in = False
        
        # Create a counter variable to keep track of which Binding site we are.
        count = 0
        
        # Open file Info_R2.
        with open(DIR + "/" + info_R2,"r") as f:
            # For each line in the file.
            for line in f:
                # Split the data and save as a list (t).
                t = line.split()
                # if Get_in is true
                if get_in:
                    #Put Get_in as False.
                    get_in = False
                    
                    # Add the Reactive_ligands numbers and add in which dynamic we are as a key in the dictionary Binding_sites and the list of residues as a value in the Binding_sites dictionary.
                    Binding_sites["lig_" + str(Reactive_ligands[count]) + "_dyn_" + str(i+1)] = t
                    
                    # Save also the Reactive ligands numbers in the Binding_each_dyn dictionary (having the dyn we have as a key f.ex. dyn_1 or dyn_2 ).
                    Binding_each_dyn[DIR[-5:]] = Reactive_ligands
                    
                    # Sum one unit to count.
                    count += 1
                
                # For each element in the line
                for ele in t:
                    # If the element in the line is equal to "RESIDUES"
                    if ele == "RESIDUES":
                        # get_in is True (to be able to know the residues of each Binding site).
                        get_in = True   

        # Define that the count is 3 because we start adding cases in line 3 (due to the first lines are the header).
        count = 3
        
        # Iterate all the binding sites we have in each dynamic.
        for dyn in Binding_each_dyn:
            # If the dynamic is equal to the Directory we are.
            if dyn == DIR[-5:]:
                # Iterate all the binding sites for the Dynamic we are.
                for j in Binding_each_dyn[dyn]:
                    # Add in column D the ligand we have as reactive f.ex. lig_612.
                    hoja["D" + str(i*15 + count)] = "lig_" + j
                    
                    # For each Binding site in the dictionary
                    for t in Binding_sites:
                        # If the BS is equal to the ligand number + the dyn we are.
                        if t == ("lig_" + j + "_" + DIR[-5:]):
                            # Define a res string to save all the Residues as a string.
                            res = ""
                            # Iterate the values of the Binding site dictionary of t key.
                            for y in Binding_sites[t]:
                                # Add the element one by one in the string called res plus one separation space.
                                res += str(y) + " "
                            # Add the Residues of the ligand in column O.
                            hoja["O" + str(i*15 + count)] = res
                            
                            # Break to exit the loop of the Binding sites and follow with the iteration of (for j in Binding_each_dyn[dyn]:)
                            break

                    # Add one value to the count.
                    count += 1
                    
    ################################################
    ### Compare binding sites and fill the Excel ###
    ################################################
    # Print the Comparision of the binding sites.
    print("3. Compare Binding Sites.")    
     
    all_BS = list()
    
    # Iterate the all the BS we have 
    for lig in Binding_sites:
        # Add each BS to the list of all BS.
        all_BS.append(Binding_sites[lig])
        
    # Define that the count is equal to 1
    count = 1

   # Create a dictionary called dict_BS where we will save the BS we have. Because accross different dyns we could have the same BS (that is compared using the interaction RESIDUES of each ligand)
    dict_BS = {}
    
    # Iterate until the length of all_BS is != of 0.
    while len(all_BS) > 0:
        # Define the length of the list all_BS.
        length = len(all_BS)
        
        # Define the included list to know which ligands have the same BS as the reference.
        included = []
        
        # Define a boolean as true called first 
        first = True
        
        # Iterate all the BS with an index called i.
        for i in range(0, length):
            # If is the first case.
            if first:
                # Save First as False.
                first = False
                
                # include the first case (first BS) in the included list.
                included.append(all_BS[i])
                
            # Iterate all the BS with an index called j.    
            for j in range(0, length):
                # If the indices are different
                if i != j:
                    # Compare the length od the two BS and catch the minimum length.
                    minimo = min(len(all_BS[i]), len(all_BS[j]))
                    
                    # Compare the two Binding sites and save the common residues in compare.
                    compare = set(all_BS[i]) & set(all_BS[j])
                    
                    # If the length of compare is equal to half of the minimum length of the two BS.
                    if len(compare) > (minimo/2):
                        # Include that BS as included. 
                        included.append(all_BS[j])
            # Break the i loop.
            break

        # Save the BS in the dict_BS dictionary as "BS1", "BS2", etc as a key and the list of ligands involved in the included list as a value.
        dict_BS["BS" + str(count)] = included
        
        # Iterate the included BS.
        for i in included:        
            # Remove each included BS of the all_BS list.
            all_BS.remove(i)
            
        # Add one value to the count.
        count += 1
        
    # Define i = -1.
    i = -1
    
    # Define a test variable as 0
    test = 0
    
    # Iterate each dynamic in the dictionary of Binding sites of each dynamic.
    for dyn in Binding_each_dyn:
        # Add one value to i (to know in which dynamic we are)
        i += 1
        
        # Define count (to count the lines we are adding at each time).
        count = 2
        
        # Iterate the ligands we have for each dynamic.
        for lignd in Binding_each_dyn[dyn]:
            # Add one value to count.
            count += 1

            # Iterate the BS in the dict_BS dictionary.
            for BS in dict_BS:
                # For each ligands involved in the BS analysed.
                for equal_BS in dict_BS[BS]:
                    # If the ligand BS is equal to the Binding sites we have in the line.
                    if (equal_BS == Binding_sites["lig_" + lignd + "_" + dyn]):
                        # Add one value to test variable.
                        test += 1
                        
                        # Save in column I the BS we have for the ligand in the line we are.
                        hoja["I" + str(i*15 + count)] = BS
                        
                        # Break the foor loop.
                        break
                    # If not
                    else:
                        # Contine with the next BS.
                        continue
                    # If it is equal to the BS, break the for BS in dict_BS loop.
                    break
    
    
    ###################################################
    ### Catch LIE Energy, residence time and DeltaG ###
    ###################################################
    # Print the next step.
    print("4. Catch LIE Energy, residence time and Delta G and fill the table")


    ##### CATCH LIE Energy, Delta G and Residence time of Info_R4_..._modified
    
    # For each dynamic going through an index "i".
    for i in range(0, num_dyns):
        # Get the directory of the dyn we are.
        DIR = DIR[:-1] + str(i+1)
        
        # Define a dictionary of the ligands.
        ligands = dict()
        
        # Open th efile Info_R4.
        with open(DIR + "/" + info_R4,"r") as f:            
            # For each line in the file.
            for line in f:
                # Split the line and save the line in list t.
                t = line.split()
                # If t has Working and lig_
                if t[1] == "Working" and t[6][0:4] == "lig_":
                    # Save into ligands directory the ligand we have found (as a key) and an empty list as a value.
                    ligands[t[6]] = []
                    
                    # Save the lig we have found into variable lig.
                    lig = t[6]
                    
                # If t has Average
                if (t[0] == "Average"):
                    # If t has LIE
                    if (t[2] == "LIE"):
                        # Save the value of LIE Average Energy in one variable called LIE_Energy_100ns 
                        LIE_Energy_100ns = str(t[9])
                        
                        # Put that variable in the list of ligands defined before.
                        ligands[lig].append(LIE_Energy_100ns)
                    
                    # If t has Delta and != 100
                    if (t[1] == "Delta" and  t[6] != "100"):
                        # Save the value of Delta_G of the last Xns in Delta_G_last variable.
                        Delta_G_last = str(t[8])
                        
                        # Put the varibale in the list of ligands defined before.
                        ligands[lig].append(Delta_G_last)
                    
                    # If t has Delta and 100
                    if (t[1] == "Delta" and t[6] == "100"):
                        # Save the value of Delta_G of the last 100 ns in Delta_G_last variable.
                        Delta_G_100ns = str(t[8])
                        
                        # Put the varibale in the list of ligands defined before.
                        ligands[lig].append(Delta_G_100ns)
                
                # If t has Tiempo and residencia.
                if (t[0] == "Tiempo" and t[2] == "residencia"):
                    # If t has RMSD
                    if (t[3] == "RMSD"):
                        # Save the residence time calculated using RMSD in time_residence variable.
                        time_residence = str(t[5])
                        
                        # Put the variable in the list of ligands defined before.
                        ligands[lig].append(time_residence)
                    else: # If not
                        # Save the residence time calculated using DeltaG in time_residence variable.
                        time_residence = str(t[4])
                        
                        # Put the variable in the list of ligands defined before.
                        ligands[lig].append(time_residence)

        # Define count = 3 (because the 3rd firsts columns are the header).
        count = 3
        
        # Print the directory we are.
        print(DIR)
        
        # For each lig in the ligands key.
        for lig in sorted(ligands):
            # Save the data calculated (Delta G, residence time of Delta G, residence time of RMSD and avg LIE Energy).
            # column E = avg LIE Energy.
            hoja["E" + str(i*15 + count)] = round(float(ligands[lig][3]),2)
            
            # column F = Residence time using RMSD.
            hoja["F" + str(i*15 + count)] = round(float(ligands[lig][2]),2)
           
            # column G = average Delta G.
            hoja["G" + str(i*15 + count)] = round(float(ligands[lig][1]),2)
            
            # column H = Residence time using DeltaG.
            hoja["H" + str(i*15 + count)] = round(float(ligands[lig][0]),2)
            
            # Add one value to count (go to next line).
            count += 1


    ###################################################
    ###      Calculate SASA and fill the table      ###
    ###################################################
    # Print the next step.
    print("5. Calculate SASA of LAST Snapshot (Solvent-accessible surface area) and fill the table")
  
    # Iterate for each molecular dynamic using index i.
    for i in range(0, num_dyns):
        # Get the directory of each dyn.
        DIR = DIR[:-1] + str(i+1)
        
        # Create a ligands dictionary.
        ligands = dict()
        
        # Create an enter boolean as False.
        enter = False
        
        # Open th efile Info_R2.
        with open(DIR + "/" + info_R2,"r") as f:            
            # For each line in the file.
            for line in f:
                # Split the line and save the line as a list (t).
                t = line.split()
                
                # If the length of the line is greater than 6, follow.
                if len(t) > 6:
                    # If t has REACTICE and trajectory.
                    if t[6] == "REACTIVE" and t[7] == "trajectory":
                        # Put enter boolean as True
                        enter = True
                    
                    # If t has  Ligand and SASA and enter is True.
                    if enter and t[1] == "Ligand" and t[4] == "SASA":
                        # Save into ligands the ligand number as a key and a list of SASA and %SASA as a list.
                        ligands["lig_" + t[2]] = [t[6],t[9]]
                        
                        # Put enter boolean as false.
                        enter = False

        # Define count as 3 (first 3 lines are the header).
        count = 3
        
        # Iterate the ligands in the dictionary called ligands.
        for lig in sorted(ligands):
            # Get SASA
            num1 = (round(float(ligands[lig][0]),2))
            
            # Get % SASA of the ligand.
            num2 = (round(float(ligands[lig][1]),2))
            
            # If one those values of SASA are less than -20.
            if num1 < -20 or num2 < -20:
                # Put into the excel file that has occurred an error when reading the SASA.
                hoja["J" + str(i*15 + count)] = "Error... Go to file " + lig + "_SASA.out and change the Surface that is -1.0000"
                
                # Put that has occurred an error when reading the SASA.
                hoja["K" + str(i*15 + count)] = "Error"
            else: # If not
                # Define a string to include the SASA.
                str_num1 = str(num1)
                
                # Deinfe a string to include the % SASA.
                str_num2 = str(num2)
                
                # Save into the Excel the SASA in column J.
                hoja["J" + str(i*15 + count)] = str_num1
                
                # SAve into the Excel the %SASA in column K.
                hoja["K" + str(i*15 + count)] = str_num2
            
            # Add one value to count (go to next line). 
            count += 1


    ###################################################
    ###      Calculate SASA and fill the table      ###
    ###################################################
    
    # Print the next step.  
    print("6. Calculate SASA last 100ns (Solvent-accessible surface area) and fill the table")
  
    # Iterate the dyns we have using index i.
    for i in range(0, num_dyns):
        # Get the directory of the dyn.
        DIR = DIR[:-1] + str(i+1)
        
        # Create a a ligands dictionary.
        ligands = dict()
        
        # Open info_R4 file. 
        with open(DIR + "/" + info_R4,"r") as f:            
            # For each line in the file.
            for line in f:
                # Split each line and save as a list in t.
                t = line.split()
                
                # If the length of the line is greater than 0.
                if len(t) > 0:
                    # If t has Working and lig.
                    if t[1] == "Working" and t[6][0:4] == "lig_":
                        # Deifne a list for each ligand.
                        ligands[t[6]] = []
                        
                        # Save the ligand number as lig variable.
                        lig = t[6]
                    
                    # If t has SASA, avg and last.
                    if t[0] == "SASA" and t[1] == "avg" and t[2] == "last":                        
                        # Get the SASA of the last 100ns.
                        SASA = str(t[6])
                        
                        # Save the SASA in the list created before for each ligand.
                        ligands[lig].append(SASA)
                    
                    # If t has %, SASA, avg and last.
                    if t[0] == "%" and t[1] == "SASA" and t[2] == "avg" and t[3] == "last":
                        # Get the %SASA of the last 100ns.
                        PC_SASA = str(t[7])
                        
                        # Save the %SASA in the list created before for each ligand.
                        ligands[lig].append(PC_SASA)
                             
        # Define count = 3 (due that the firsts lines are the header).        
        count = 3
        
        # Iterate each ligand in the list of ligands.
        for lig in sorted(ligands):
            # Save into the Excel file the SASA in column L.
            hoja["L" + str(i*15 + count)] = round(float(ligands[lig][0]),2)
            
            # Save into the Excel file the %SASA in column M.
            hoja["M" + str(i*15 + count)] = round(float(ligands[lig][1]),2)
            
            # Add one value to count (go to the next line).
            count += 1

   
    ###################################################
    ###      Calculate SASA and fill the table      ###
    ###################################################
    
    # Print the next step.
    print("7. Calculate Number of hydrogen bonds and fill the table")
    
    # Iterate the dyns we have using index i.
    for i in range(0, num_dyns):
        # Get the directory of the dyn we are iterating.
        DIR = DIR[:-1] + str(i+1)
        
        # Create a dicitonary called ligands.
        ligands = dict()
        
        # Open the file Info_R4.
        with open(DIR + "/" + info_R4,"r") as f:            
            # Iterate each line of the file.
            for line in f:
                # Split the data of each line and save it as a lits called t.
                t = line.split()
                
                # If the length of t is greater than 0.
                if len(t) > 0:
                    # If t has Working and lig.
                    if t[1] == "Working" and t[6][0:4] == "lig_":
                        # Create a list of the ligand on the ligands dictionary.
                        ligands[t[6]] = []
                        
                        # Save the ligand number into lig variable.
                        lig = t[6]
                        
                    # If t has Number, hydrogen and bonds.
                    if t[0] == "Number" and t[2] == "hydrogen" and t[3] == "bonds":
                        # Get the number of hydrogen bonds we have for that ligand.
                        hydrogen_bonds = str(t[5])
                        
                        # Put the hydrogen bonds into the ligands dictionary.
                        ligands[lig].append(hydrogen_bonds)

        # Define the count to 3 (because the firsts lines are the header).
        count = 3
        
        # Iterate all the ligands.
        for lig in sorted(ligands):
            # Save into the Excel file the number of hydrogen_bonds obtained.
            hoja["N" + str(i*15 + count)] = int(ligands[lig][0])
            
            # Add one value to count (go to the next line):
            count += 1
    
    ############################
    ###      Save Excel      ###
    ############################
    # Print the next step.
    print("7. Save the final Excel with all the fields.")

    # Save the workbook as an Excel file.
    wb.save(str(ligand) + "_all_dyns_summary" + ".xlsx")
    
    # Print the End of the Program.
    print("File created")
    print("")
    print(">>> End ")





def cmdlineparse():
    # Parse the argument needed to read the information.
    
    parser = ArgumentParser(description="command line arguments")
    parser.add_argument("-ligand" , dest="ligand"   , required=True, help=" To know which ligand we are using")
    parser.add_argument("-sim_time" , dest="sim_time"   , required=True, help=" Simulation time of that dynamic")
    parser.add_argument("-num_dyns" , dest="num_dyns"   , required=True, help="Number of dynamics we have in the trajectory")
    parser.add_argument("-molecule_receptor" , dest="molecule_receptor"   , required=True, help=" Name of the receptor molecule")
    parser.add_argument("-info_R2" , dest="info_R2"   , required=True, help="Information file Info_R2_fdMD_ReactiveTraj_26082020")
    parser.add_argument("-info_R4" , dest="info_R4"   , required=True, help="Information file Info_R4_fdMD_MMGBSA_Analize_25012021")
    parser.add_argument("-num_ns_anal" , dest="num_ns_anal"   , required=True, help="Number of ns to analize, to calculate the avg of the LIE energy")
    
    args=parser.parse_args()
    return args

# Initial information.
def prog_info() :
  t = datetime.now()
  format_date = "%d-%m-%Y %H:%M:%S"
  ti= t.strftime(format_date)
  print ( "... ... ... ... ... ... " )
  print ( "...   fdMD_Summary  ... " )
  print ( "...     ( 2021 )    ... " )
  print ( "...     ........    ... " )
  print ( "... ... ... ... ... ... " )
  print ( " ",ti )
  print ( "... ... ... ... ... ... " )
  print ( " " )

# Ending Information.
def prog_END() :
  t = datetime.now()
  format_date = "%d-%m-%Y %H:%M:%S"
  tf= t.strftime(format_date)
  print ( "... ... ... ... ... ... " )
  print ( "...   fdMD_Summary  ... " )
  print ( "...     ( 2021 )    ... " )
  print ( "...     ........    ... " )
  print ( "... ... ... ... ... ... " )
  print ( " ",tf )
  print ( "... ... ... ... ... ... " )
  print ( " " )

if __name__ == '__main__':

  # Call the Initial information.
  prog_info()
  
  # Call the main information.
  main()
  
  # Call the final information.
  prog_END()
