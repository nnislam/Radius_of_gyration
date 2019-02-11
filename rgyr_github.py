#!/anaconda3/bin/pythonw

######################################
## Radius of Gyration ##
## Naeyma Islam ##
## UNO, 6th Feb 2018 ##
######################################

#----- Importing the libraries -----#

import numpy as np
import pandas as pd
import math as mth
import matplotlib.pyplot as plt
from scipy.interpolate import spline

data_file_name = "out.data"       # .data file name along with location path
trj_file_name = "pmaa.lammpstrj"  # .lammpstrj file name along with location path
rgyr_out_file_name = "rgyr.out"   # output file name along with path to store rgyr values

#----- Reading .data file (from lammps output) -----#
cols = ['TYPE', 'MASS']
mass_of_atoms = pd.DataFrame(columns=cols)
mass_one_row = []

print ("Start reading Data file....\n")

try:
    data_file = open(data_file_name, "r")  #opens the .data file in the working directory in 'reading mode'
    search_for_text = 'Masses'      # Text to search for
    found_text = 0                  # Search text handler
    for line in data_file:          # runs the loop through each line of the file
        #print (line)
        if found_text == 1:         # when search text is found, starting from next iteration of the loop
            list = line.split()     # splits the current line and save values to list[] string array
            if list and list[0].isdigit():   # check whether the line is not empty and the first word is a number
                #mass_of_atoms.append([int(list[0]), float(list[1])]) # storing the values into the array
                mass_one_row = [int(list[0]), float(list[1])]
                new_row_df = pd.DataFrame ([mass_one_row], columns=cols)
                mass_of_atoms = mass_of_atoms.append(new_row_df)
            else:
                found_text = 0      # resetting the search text handle
                break               # stopping the loop iteration since no more processing is needed

        if search_for_text in line: # check for the search text
            found_text = 1          # Text found in the line
            data_file.readline()    # skipping the empty lines
    
    data_file.close()                  # closing the opened file

except ValueError:
    print ("Cannot find the file or wrong file format. Please check and try again.")
    system.exit(os.EX_CONFIG)

#print (mass_of_atoms)     # checking the loaded values

print ("Reading of Data file completed.\n")

print ("Start processing Trajectory file....\n")

#----- Reading .lammpstrj file (from lammps output) -----#

cols = ['TIMESTEP', 'ATOM_ID', 'TYPE', 'x', 'y', 'z', 'xbox', 'ybox', 'zbox', 'MASS']
trj_data = pd.DataFrame(columns=cols)    # Dataframe for the Timestep data
trj_data_new = pd.DataFrame(columns=cols)
trj_one_row = []
xbox = 0
ybox = 0
zbox = 0
atom_type_to_excl = 1  # Water

cols_rgyr = ['TIME', 'RGYR']
rgyr_data = pd.DataFrame(columns=cols_rgyr)  # Dataframe to store the calculated Rgyr values along with Time
new_rgyr = []

try:
    trj_file = open(trj_file_name, "r")  #opens the .data file in the working directory in 'reading mode'
    
    out_file = open (rgyr_out_file_name, "w")   # opening the out file to write rgyr values
    
    search_for_timestep = 'ITEM: TIMESTEP'      # Text to search for
    search_for_number = 'ITEM: NUMBER'
    search_for_box = 'ITEM: BOX'
    search_for_atom = 'ITEM: ATOMS'
    found_timestep = 0                  # Search text handler
    found_number = 0
    found_box = 0
    found_atom = 0
    total_atoms = -1
    atom_count = 0
    time_step = -1

    for line in trj_file:          # runs the loop through each line of the file
        #print (line)
                
        if found_timestep == 1 and time_step == -1:         # when search text is found, starting from next iteration of the loop
            list = line.split()     # splits the current line and save values to list[] string array
            if list and list[0].isdigit():   # check whether the line is not empty and the first word is a number
                time_step = int(list[0])     # storing the timstep value into the variable
                print ("->Processing Frame: %d" % time_step)
        
        if found_number == 1 and total_atoms == -1:         # when search text is found, starting from next iteration of the loop
            list = line.split()     # splits the current line and save values to list[] string array
            if list and list[0].isdigit():   # check whether the line is not empty and the first word is a number
                total_atoms = int(list[0])     # storing the timstep value into the variable
                
        if found_timestep == 1 and found_box >= 1:         # when search text is found, starting from next iteration of the loop
            list = line.split()     # splits the current line and save values to list[] string array
            #print (list)
            if found_box == 1:
                xbox = float(list[1]) - float(list[0])  # X distance for Box
                found_box = found_box + 1
            elif found_box == 2:
                ybox = float(list[1]) - float(list[0])  # Y distance for Box
                found_box = found_box + 1
            elif found_box == 3:
                zbox = float(list[1]) - float(list[0])  # Z distance for Box
                found_box = 0
            #print (xbox, ybox, zbox)

        if found_timestep == 1 and found_atom == 1:         # when search text is found, starting from next iteration of the loop
            atom_count = atom_count + 1
            list = line.split()     # splits the current line and save values to list[] string array
            #print (list)
            
            if list and list[0].isdigit():   # check whether the line is not empty and the first word is a number
                atom_type = int(list[1])
                if atom_type != atom_type_to_excl:      # Process all atoms except Atom_Type_to_exclude (water)

                    x = float(list[2])   
                    y = float(list[3])  
                    z = float(list[4])  
                    atom_mass = float(mass_of_atoms.loc[mass_of_atoms['TYPE'] == atom_type].MASS)
                    
                    trj_one_row = [time_step, int(list[0]), int(list[1]), x, y, z, xbox, ybox, zbox, atom_mass]   # loading one row data
                    new_row_df = pd.DataFrame ([trj_one_row], columns=cols, index=[int(list[0])])          # making dataframe for loaded row
                    trj_data = trj_data.append(new_row_df)                           # add new row into the Dataframe
                    trj_one_row = []                                                 # resetting row placeholder

        if total_atoms == atom_count:
            # When reading of one Timestep/Frame is done, calculate the Rgyr for that Frame
            # print ("Calculating Radius of Gyration for each Frame....\n")

            #trj_data = trj_data.merge(mass_of_atoms, on=['TYPE'])      # Adding MASS value for each ATOM
            #trj_data = trj_data.sort_values(by=['TIMESTEP','ATOM_ID']) # sorting the dataset

            #print (trj_data)

            #i = 0
            j = 1
            x_new = []
            y_new = []
            z_new = []
            x_new.insert(0, trj_data.iloc[0,3])
            y_new.insert(0, trj_data.iloc[0,4])
            z_new.insert(0, trj_data.iloc[0,5])
            
            new_row = [trj_data.iloc[0,0], trj_data.iloc[0,1], trj_data.iloc[0,2], x_new[0], y_new[0], z_new[0], trj_data.iloc[0,6], trj_data.iloc[0,7], trj_data.iloc[0,8], trj_data.iloc[0,9]]   # loading one row data
            new_df = pd.DataFrame ([new_row], columns=cols)          # making dataframe for loaded row
            trj_data_new = trj_data_new.append(new_df)               # add new row into the Dataframe
            
            for i in range(int(trj_data.shape[0] - 1)):
                #xi = trj_data.iloc[i,3]
                xj = trj_data.iloc[j,3]
                new_xj = xj - xbox * round ((xj - x_new[j-1])/xbox)
                x_new.insert(j, new_xj)

                #if ((xj - xi) < -xbox/2) or ((xj - xi) > xbox/2):        # boundary checking for x
                #trj_data.at[j, 'x'] = xj - xbox * round((xj - xi)/xbox)

                #yi = trj_data.iloc[i,4]
                yj = trj_data.iloc[j,4]
                new_yj = yj - ybox * round ((yj - y_new[j-1])/ybox)
                y_new.insert(j, new_yj)

                #if ((yj - yi) < -ybox/2) or ((yj - yi) > ybox/2):        # boundary checking for y
                #trj_data.at[j, 'y'] = yj - ybox * round((yj - yi)/ybox)

                #zi = trj_data.iloc[i,5]
                zj = trj_data.iloc[j,5]
                new_zj = zj - zbox * round ((zj - z_new[j-1])/zbox)
                z_new.insert(j, new_zj)

                #if ((zj - zi) < -zbox/2) or ((zj - zi) > zbox/2):        # boundary checking for z
                #trj_data.at[j, 'z'] = zj - zbox * round((zj - zi)/zbox)

                #print (xi, xj)
                #i = i + 1
                
                new_row = [trj_data.iloc[j,0], trj_data.iloc[j,1], trj_data.iloc[j,2], x_new[j], y_new[j], z_new[j], trj_data.iloc[j,6], trj_data.iloc[j,7], trj_data.iloc[j,8], trj_data.iloc[j,9]]   # loading one row data
                new_df = pd.DataFrame ([new_row], columns=cols)          # making dataframe for loaded row
                trj_data_new = trj_data_new.append(new_df)               # add new row into the Dataframe
                #trj_one_row = []
                
                #print (trj_data)
                #print (trj_data_new)
                
                j = j + 1
            
            #print (trj_data)
            #print ("after boundary checking")
            #print (trj_data_new)
            total_mass = trj_data_new['MASS'].agg(np.sum)  # calculating mass
            #print (total_mass)

                # finding rcm (x, y, z) for each timstep
            rcmx = ((trj_data_new.x * trj_data_new.MASS).agg(np.sum))/total_mass
            rcmy = ((trj_data_new.y * trj_data_new.MASS).agg(np.sum))/total_mass
            rcmz = ((trj_data_new.z * trj_data_new.MASS).agg(np.sum))/total_mass
            #print (rcmx, rcmy, rcmz)

                # Radius of Gyration for the frame
            frame_sum = (((trj_data_new.x - rcmx).agg(np.square) + (trj_data_new.y - rcmy).agg(np.square) +
          (trj_data_new.z - rcmz).agg(np.square)) * trj_data_new.MASS).agg(np.sum)

            rgyr = (mth.sqrt (frame_sum/total_mass))*0.1
            #print (rgyr)
            
            time_data = (time_step * 5)/1000000
            #print (time_data)

            out_file.write(str(time_data) + "\t" +  str(rgyr) + "\n")      # saving rgyr value into the file

            new_rgyr = pd.DataFrame ([[time_data, rgyr]], columns=cols_rgyr)
            rgyr_data = rgyr_data.append(new_rgyr)
            
            print ("---Time: %f, Rgyr: %f" %(time_data, rgyr))

              # reset everything and get ready to process next Frame/Timestep
            found_timestep = 0      # resetting the search text handles and time_step value, box values
            found_number = 0
            found_box == 0
            found_atom = 0
            time_step = -1
            total_atoms = -1
            atom_count = 0
            xbox = 0
            ybox = 0
            zbox = 0
            trj_data = trj_data.iloc[0:0]
            trj_data_new = trj_data_new.iloc[0:0]
            #break                  # stopping the loop iteration since no more processing is needed

        if search_for_timestep in line: # check for the search text
            found_timestep = 1          # Text found in the line
            #trj_file.readline()        # skipping the empty lines
            
        if search_for_number in line:      # check for the search text
            found_number = 1               # Text found in the line
            
        if search_for_box in line:      # check for the search text
            found_box = 1               # Text found in the line

        if search_for_atom in line:     # check for the search text
            found_atom = 1              # Text found in the line
    
    out_file.close()                    # closing the out file
    trj_file.close()                    # closing the opend trj file
    
except ValueError:
    print ("Cannot find the file or wrong file format. Please check and try again.")
    system.exit(os.EX_CONFIG)

#print (rgyr_data)

#print (trj_data.head(5))          # checking the loaded values (samples)
print ("\nProcessing of Trajectory file completed.\n")

print ("Ploting the Radius of Gyration values....\n")

plt.plot(rgyr_data['TIME'], rgyr_data['RGYR'])

plt.ylabel('Rg(nm)')
plt.xlabel('Time(ns)')
plt.plot
plt.show()

# smooth graph
#x = np.array (rgyr_data['TIME'])
#y = np.array (rgyr_data['RGYR'])

#x_smooth = np.linspace (x.min(), x.max(), 100)
#y_smooth = spline(x, y, x_smooth)

#plt.ylabel('Rg(nm)')
#plt.xlabel('Time(ns)')
#plt.plot (x_smooth, y_smooth)
#plt.show()
