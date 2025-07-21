#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
Created on 2024/03/28 10:27:22

filename :   reformat_rsqsim_catalog.py
@author  :   msc 
@contact :   msolares@uoregon.edu
'''

#%%# reformat_rsqsim_catalog.py -------------------------------------------------------------#

###############################################################################
# Reformats the RSQSIM catalog for use with FakeQuakes.
###############################################################################


### IMPORT LIBRARIES ----------- #

### DEFINE VARIABLES ----------- #

### DEFINE FUNCTIONS ----------- #


#%%############################################################################


def extract_event_ids(workin_dir, output_dir, output_filename='rsqsim_evid_list.txt', gns_eventname="event", svo=0):

    """
    Extracts event IDs from files in a specific directory and saves them to a text file.

    Parameters:
        workin_dir (str): Directory containing input files (with .txt extension).
        output_dir (str): Directory to save the output file, if enabled.
        output_filename (str): Name of the output text file for saving event IDs.
        svo (int): Save output option. 1 to save the output, 0 to not save (default).
        gns_eventname (str): Prefix used in filenames to identify event IDs (default is "event").

    Returns:
        list: Sorted list of extracted event IDs.
    """

    import numpy as np
    import glob
    import os
    

    
    # Initialize an empty list to store event IDs
    evidlst = []

    # Find all .txt files in the working directory
    txt_files = glob.glob(os.path.join(workin_dir, "*.txt"))

    print(f"\nGetting event IDs from files ... found [{len(txt_files)}] events\n")
 
    for filepath in txt_files:
        filename = os.path.basename(filepath)
        print 
        try:
            # Extract the part that starts with gns_eventname (e.g., "event1234_xyz.txt")
            if gns_eventname in filename:
                event_part = filename.split('_')[0]
                evidlst.append(event_part)
        except Exception as e:
            print(f"(!) Failed to process file '{filename}': {e}") 

    # Sort IDs numerically (extract number after gns_eventname)
    try:
        sorted_event_ids = sorted(evidlst, key=lambda x: int(x.split(gns_eventname)[1]))
    except Exception as e:
        print(f"(!) Sorting error: {e}")
        sorted_event_ids = sorted(evidlst)

    if svo == 1:
        output_path = os.path.join(output_dir, output_filename)
        try:
            np.savetxt(output_path, sorted_event_ids, fmt='%s')
            print(f" Saved file to: {output_path}")
        except Exception as e:
            print(f"(!) Error saving file: {e}")
   
    elif svo == 0:
        # print(" Output text file not saved (svo=0).")
        pass
   
    else:
        print(" Invalid save option. Output text file not saved.")

    print("done.")
    return sorted_event_ids



def create_log_files(log_file, txt_file, output_dir):
    """
    Create log files for each event ID listed in a text file or list.

    Parameters:
        log_file (str): Path to the source log file.
        txt_file (str or list): Path to the text file containing event IDs,
                                 or a list containing event IDs.
        output_dir (str): Directory to save the new log files.

    Returns:
        None
    """

    import os 

    print("\nCreating .LOG files per rupture using new event id ...\n")

    if isinstance(txt_file, str):
        # If txt_file is a string (path to a file), read event IDs from the file
        with open(txt_file, 'r') as file:
            event_ids = [line.strip() for line in file]
    elif isinstance(txt_file, list):
        # If txt_file is a list, use it directly as event IDs
        event_ids = txt_file
    else:
        raise ValueError(" txt_file must be a path to a text file or a list.")

    # Open the source log file and create log files for each event ID
    with open(log_file, 'r') as original_file:
        for event_id in event_ids:
            # Create a new file with the filename
            with open(os.path.join(output_dir, f"{event_id}.log"), 'w') as new_file:
                # Copy content line by line
                original_file.seek(0)  # Reset file pointer to the beginning
                for line in original_file:
                    new_file.write(line)
    
    print("done.")



def create_rupt_list(event_list, output_dir, output_filename):
    """
    Create a list of rupture files with '.rupt' appended to each event ID and save it as a text file.

    Parameters:
        event_list (list): List of event IDs.
        output_dir (str): Directory where the output file will be saved.
        output_filename (str): Name of the output file. Default is 'rupture.list'.

    Returns:
        None
    """

    import os

    print("\nCreating rupture list file ...\n")

    # Append '.rupt' to each item in the list
    rupt_files = [event_id + '.rupt' for event_id in event_list]

    # Specify the filename for the output file
    output_file = output_filename

    # Combine the output directory and filename
    output_path = os.path.join(output_dir, output_file)

    # Save the modified list as a text file
    with open(output_path, 'w') as file:
        for item in rupt_files:
            file.write("%s\n" % item)

    print("done.")



def convert_from_crs(x, y, z, source_epsg, target_epsg):
    """
    Convert coordinates from a source CRS to a target CRS.

    Parameters:
        x (float): X coordinate.
        y (float): Y coordinate.
        z (float): Z coordinate.
        source_epsg (str): EPSG code of the source CRS.
        target_epsg (str): EPSG code of the target CRS.

    Returns:
        tuple: A tuple containing the converted longitude, latitude, and depth.
    """

    from pyproj import Transformer

    # Define transformer from source CRS to destination CRS
    transformer = Transformer.from_crs(source_epsg, target_epsg)

    # Convert coordinates
    lon, lat, z_dest = transformer.transform(x, y, z)

    # Calculate depth
    depth = z_dest * -0.001  # Converting depth to kilometers

    # Return lon, lat, and depth as a tuple
    return lon, lat, depth



def process_catalog(catalog_file, evid_list_sorted, output_dir, output_filename, svo=0):
    """
    Process a catalog file to extract magnitudes (Mw) and locations (latitude, longitude, depth).

    Parameters:
        catalog_file (str): Path to the catalog file in CSV format.
        evid_list_sorted (list): Sorted list of event IDs corresponding to the catalog entries.
        output_dir (str): Directory where the output file will be saved.
        output_filename (str): Name of the output text file.
        svo (int, optional): Save output option. Set to 1 to save the output file, 0 to not save. Defaults to 0.

    Returns:
        numpy.ndarray or None: Processed data array if `svo` is 1, otherwise None.

    Notes:
        - The catalog file is expected to be in CSV format with columns representing:
            - Magnitude (Mw)
            - X coordinate (longitude)
            - Y coordinate (latitude)
            - Z coordinate (depth)
        - The coordinates are assumed to be in the NZTM (New Zealand Transverse Mercator) projection.
        - The function applies coordinate conversion to convert coordinates from NZTM to WGS 84 (latitude and longitude).
        - The processed data is saved to a text file with columns for event ID, magnitude, latitude, longitude, and depth.
        - If `svo` is set to 1, the output file is saved; otherwise, it's skipped.

    Examples:
        # Process catalog file and save output
        process_catalog("catalog.csv", ['event1', 'event2'], "output_dir", "output.txt", svo=1)

        # Process catalog file without saving output
        process_catalog("catalog.csv", ['event1', 'event2'], "output_dir", "output.txt", svo=0)
    """
    
    import numpy as np
    import os
    
    print(f"\nProcessing catalog in: {catalog_file.split('/')[-2]}/ ...\n")
    
    # Extract mags (Mw) + loc (x,y,z) from gns file (.cvs) for each event using evid from previous routine
    catalog = np.genfromtxt(catalog_file, 
                            skip_header=1, 
                            delimiter=",", 
                            usecols=(3,4,5,6),
                            dtype=[('col1', 'f2'), ('col2', 'f'), ('col3', 'f'), ('col4', 'f')],
    ) 

    # empty list to append elements
    mwlst = []
    lalst = []
    lolst = []
    zzlst = []

    # Define the EPSG codes for the source and target CRS
    source_epsg = "EPSG:2193"  # Source CRS (NZTM)
    target_epsg = "EPSG:4326"  # Destination CRS (WGS 84)

    # Iterate over each row in the array
    for row in catalog:
        m, x, y, z = row
        # Apply the conversion function to the x, y, z coordinates
        lat, lon, depth = convert_from_crs(y, x, z, source_epsg, target_epsg)
        mwlst.append("%.2f"%m)
        lalst.append("%.4f"%lat)
        lolst.append("%.4f"%lon)
        zzlst.append("%.4f"%depth)
        # Print or store the converted coordinates and depth
        print(f" Magnitude: {m:.2f}, Latitude: {lat:.2f}, Longitude: {lon:.2f}, Depth: {depth:.2f} kilometers")

    # create output file by concatinating the two list (now columns) into one file
    ofile_magloc = np.c_[evid_list_sorted, mwlst, lolst, lalst, zzlst]

    if svo == 1:
        # Save list to a text file
        output_path = os.path.join(output_dir, output_filename)
        # Add a header
        header = "EventID, Magnitude, Latitude, Longitude, Depth"
        np.savetxt(output_path, ofile_magloc, delimiter=',', header=header, comments='#', fmt='%s')
        print(f"\n Saved file to: {output_path}")  
    elif svo == 0:
        pass
    else:
        print(" Invalid save option. List will not be saved.")
        
    return ofile_magloc



def process_ruptures(dfile, rsqsim_rupts_dir, fault, mesh, rake_avg=90, dx=0.01, fqs_dir='', output_filename='', output_plots_dir='', svo=0):
    """
    modified from: DMelgar @UO

    Process rupture and create .rupt files in fakequakes format.

    Parameters:
        dfile (str or numpy.ndarray): Data file containing a list of event IDs and magnitudes.
                                       If it's a string, it should be a path to a TXT or CSV file (sep=',').
                                       If it's a numpy array, it should contain event IDs and magnitudes.
        rsqsim_rupt (str): Path to the directory where RSQSIM rupture files are located. 
                       Default is an empty string.
        fault (str): Path to the directory where to find fault model file.
                       Default is an empty string.
        mesh (str): Path to the directory where to find mesh file.
                       Default is an empty string. 
        dx (float): Spacing value to create grid - smaller number, finer grid and more time it takes to complete. 
                       Default is 0.01.                                              
        rake_avg (int): Average rake value in degrees.  
                       Default is 90 for a thrust fault.        
        fqs_dir (str): Path to the directory where rupture output files will be saved.
                       Default is an empty string.
        output_filename (str): Name of the output plot files.
                               Default is an empty string.
        output_plots_dir (str): Directory to save output plots.
                                 Default is an empty string.
        svo (int): Save output option. 1 to save the output, 0 to not save.
                   Default is 0.

    Returns:
        None
    """

    import numpy as np
    import glob
    import matplotlib.pyplot as plt
    from matplotlib.path import Path
    from mudpy.fakequakes import get_rise_times
    from scipy.interpolate import griddata
    

    print("\nCreating .RUPT files ...\n")

    # data file with list [event id, mags] to data 
    if isinstance(dfile, str):
        # If txt_file is a string (path to a file), read event IDs from the file
        array = np.genfromtxt(dfile, delimiter=',', usecols=(0,1), dtype=str)
    elif isinstance(dfile, np.ndarray):
        # If txt_file is a list, use it directly as event IDs
        array = dfile
    else:
        raise ValueError("txt_file must be a path to a text file or a list.")
    
    # create a loop for event id (ev_id) & mags (Mw) 
    for (event, mag) in zip(array[:,0],array[:,1]):
        ev_id = event
        Mw = float(mag)
        print(f"\nEvent id: {ev_id} \tMagnitude: {Mw}\n")

        ## path to rupt file 
        path = glob.glob(f"{rsqsim_rupts_dir}/{ev_id}_with_zeros.txt")
        # print(path) # uncomment for debugging

        ## load data from a text file, with missing values handled as specified
        r=np.genfromtxt(path[0]) 

        ## change outfile extension 
        out_rupt=f"{fqs_dir}{ev_id}.rupt"        
        # print(out_rupt) # uncomment for debugging

        ## reformating rupture file lon,z   
        i=np.where(r[:,0]<0)[0]
        r[i,0]+=360
        r[:,2]*=-0.001

        #### Re sampled model  ####

        ## fault geometry
        f=np.genfromtxt(f"{fault}")        
        ## moment centroid 
        mshout=np.genfromtxt(f"{mesh}")
        #Things I need
        slip=r[:,9]
        onset=r[:,11]

        #Now I need to average into the bigger fault

        #Make regular grid of slip, dx=0.01 by default
        grid_x=np.arange(r[:,0].min(),r[:,0].max(),dx)
        grid_y=np.arange(r[:,1].min(),r[:,1].max(),dx)
        X,Y=np.meshgrid(grid_x,grid_y)
        grid_slip = griddata(r[:,0:2], slip, (X,Y), method='linear')
        grid_onset = griddata(r[:,0:2], onset, (X,Y), method='nearest')

        #Unravel this grid
        grid_x=X.ravel()
        grid_y=Y.ravel()
        grid_slip=grid_slip.ravel()
        grid_onset=grid_onset.ravel()

        #Look at mshout, make a path, get point inside and average
        average_slip=np.zeros(len(f))
        average_onset=np.zeros(len(f))
        for k in range(len(mshout)):
            
            if k % 100 == 0:
                print('Working on point %d of %d' % (k,len(mshout)))
            points=np.r_[np.expand_dims(mshout[k,1:3],0),np.expand_dims(mshout[k,4:6],0),np.expand_dims(mshout[k,7:9],0),np.expand_dims(mshout[k,1:3],0)]
        
            #Check for misbehaving points
            i=np.where(points[:,0] > 360)[0]
            points[i,0] -= 360
            i=np.where(points[:,0] < 0)[0]
            points[i,0] += 360
            
            #Make path
            p=Path(points)
            
            #WHich are in
            i=np.where(p.contains_points(np.c_[grid_x,grid_y])==True)[0]
            if len(i)>0:
                if grid_slip[i].mean() == np.nan:
                    print(k)
                average_slip[k]=grid_slip[i].mean()
                average_onset[k]=grid_onset[i].max()
            
        #Clean up nans
        i=np.where(np.isnan(average_slip) == True)[0]
        average_slip[i]=0

        #Rise time stuff only use subfaults that are non zero
        i_non_zero=np.where(average_slip>0)[0]

        M0 = 10**(1.5*Mw+9.1)
        average_rake=np.deg2rad(rake_avg*np.ones(len(f)))
        rise_times=np.zeros(len(f))
        rise_times[i_non_zero] = get_rise_times(M0,average_slip[i_non_zero],f[i_non_zero,:],[8,14],average_rake[i_non_zero],rise_time_std=0.1)
        rise_times=np.expand_dims(rise_times,1)
        average_onset=np.expand_dims(average_onset,1)


        #ss and ds slip
        ss = average_slip*np.cos(average_rake)
        ds = average_slip*np.sin(average_rake)
        ss=np.expand_dims(ss,1)
        ds=np.expand_dims(ds,1)
        mu=30e9*np.ones((len(f),1))


        #write to file
        fmt='%d\t%12.5f\t%12.5f\t%9.4f\t%9.2f\t%9.2f\t%.1f\t%8.2f\t%8.2f\t%8.2f\t%10.2f\t%10.2f\t%8.2f\t%.2e'
        out=np.c_[f[:,0:7],rise_times,ss,ds,f[:,8:10],average_onset,mu]
        np.savetxt(out_rupt,out,fmt=fmt,header='# No,lon,lat,z(km),strike,dip,rise,slip_duration(s),ss-slip(m),ds-slip(m),ss_len(m),ds_len(m),rupt_time(s),rigidity(Pa)')

        #plot results
        plt.figure()
        plt.subplot(131)
        plt.scatter(out[:,1],out[:,2],c=out[:,9],s=10,cmap="jet")
        plt.colorbar(label='ds slip (m) ')
        # plt.ylim([-43,-40])


        plt.subplot(132)
        plt.scatter(out[:,1],out[:,2],c=out[:,7],s=10,cmap="jet")
        plt.colorbar(label='rise time (s) ')
        ax = plt.gca()
        ax.get_yaxis().set_visible(False)
        # plt.ylim([-43,-40])

                    
        plt.subplot(133)
        plt.scatter(out[:,1],out[:,2],c=out[:,12],s=10,cmap="jet")
        plt.colorbar(label='onset time (s) ')
        ax = plt.gca()
        ax.get_yaxis().set_visible(False)
        # plt.ylim([-43,-40])

        plt.subplots_adjust(wspace= 0.5)

        if svo == 1: 
            Mw_str = str(Mw).replace('.', '_')
            plt.savefig(f"{output_plots_dir}plot_rupt_event{ev_id.split(output_filename)[-1]}_Mw{Mw_str}.png")
        
        # Close the figure
        plt.close()
    
    print("\ndone.")




