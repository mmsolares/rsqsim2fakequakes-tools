#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
Created on 2024/03/28 10:55:11

filename :   main_script.py
@author  :   msc
@contact :   msolares@uoregon.edu
'''

#%%# main_script.py -------------------------------------------------------------#

###############################################################################
# Main script to run the rsqsim2fakequakes.
###############################################################################


### IMPORT LIBRARIES ----------- #

import os
import sys

# This assumes we are in the rsqsim2fakequakes directory 
base_dir = os.path.dirname(os.path.abspath(__file__))

# Build the relative path
find_r2fqs = os.path.join(base_dir, 'src', 'python')
sys.path.append(find_r2fqs)

import reformat_rsqsim_catalog as reformat # type: ignore


### DEFINE VARIABLES ----------- #

###### MODIFY THIS BLOCK OF VARIABLES AS SEE FIT (!)


DAT_user = 'data_example'                                        #! data and auxiliary files directory 
CAT_user = 'rsqsim_catalog_example'                              #! which earthquake catalog directory to use?
RUP_user = 'rsqsim_outputs_example'                              #! where to find rsqsim rupture models?
AUX_user = 'fqs_input_files'                                     #! auxiliary files related to fakequakes
PNM_user = 'fakequakes_project_example'                          #! fakequakes project name; directory where to save fakequakes outputs

########################## GLOBALS & PATHS ####################################
                                     
HOME = base_dir                                                  #: where is home?

# path to data and auxiliary files
DATA_DIR_USER = f"{HOME}/{DAT_user}"                             #: path to where rsqsim outputs/catalog are located
AUXF_DIR_USER = f"{DATA_DIR_USER}/{AUX_user}"                    #: path to where auxiliary files for fakequakes are located
RSQSIM_DIR_RUPT = f"{DATA_DIR_USER}/{RUP_user}"                  #: path to where rsqsim rupture models are located

# path to fakequakes (FQs) directory
FQs_PROJ_USER = f"{HOME}/{PNM_user}"                             #: path to fakequakes project
FQs_PDIR_RUPT = f"{FQs_PROJ_USER}/output/ruptures/"              #: path to FQS datafiles in ruptures dir
FQs_PDIR_SUMM = f"{FQs_PROJ_USER}/output/waveforms/"             #: path to FQS datafiles in summary waveform results

########################## OTHER VARIABLES ####################################

catalog_filename = 'catalog_selected_events.csv'                 #! get rsqsim earthquake catalog
catalog_file = f"{DATA_DIR_USER}/{CAT_user}/{catalog_filename}"  
log_filename = 'wellington8_3.log'                               #! dummy log file to use for RSQSIM purpose
log_file = f"{AUXF_DIR_USER}/{log_filename}"
fault_filename = 'gns.fault'                                     #! fault geometry comaptible with fakequakes (would be in fakequake_project/data/model_info/ & needed to run mudpy/fakequakes code)                        
fault_user = f"{AUXF_DIR_USER}/{fault_filename}"
mesho_filename = 'gns.mshout'                                    #! fault mesh geometry comaptible with fakequakes
mesho_user = f"{AUXF_DIR_USER}/{mesho_filename}"

########################### SWITCH on/off #####################################

### WHICH tasks to run? --------------------------------------------

get_fqs_project_structure = True                                 #! True = create fakequakes project structure for testing
get_rsqsim_catalog_info_for_fqs = True                           #! True = extract RSQSim catalog info to FakeQuakes format

### SAVE outputs?                                               
svo_user = 1                                                     #! 0 = no save, 1 = save outputs 


### DEFINE FUNCTIONS ----------- #

# ... see below for functions



#%%############################################################################


###### (0) INITIALIZE FAKEQUAKES PROJECT DIRECTORY STRUCTURE FOR TESTING (OPTIONAL)

def init_fakequakes_project_for_testing(home, project_name):

    """
    Initializes the directory structure for a Fakequakes-compatible project.

    Args:
        home (str): Base path where the project should be created (e.g., "./" or "/path/to/projects/")
        project_name (str): Name of the project (will create a folder under `home`)

    Returns:
        None
    """
    from os import makedirs
    from os.path import join, exists

    proj_dir = os.path.abspath(join(home, project_name))

    # List of subdirectories to create
    subdirs = [
        "GFs",
        "GFs/static",
        "GFs/dynamic",
        "GFs/matrices",
        "GFs/STFs",
        "data/station_info",
        "data/model_info",
        "data/distances",
        "structure",
        "plots",
        "scripts",
        "forward_models",
        "output/ruptures",
        "output/statics",
        "output/waveforms",
        "logs",
        "analysis",
        "analysis/frequency",
    ]

    print(f"\nCreating project directories under: {proj_dir}\n")

    for subdir in subdirs:
        dir_path = join(proj_dir, subdir)
        if not exists(dir_path):
            makedirs(dir_path)
            print(f"  Created: {dir_path}")
        else:
            print(f"  Already exists: {dir_path}")


###### (1) AFTER running init fuction in here to test or with mudpy/fakequakes BEFORE running GFs

def run_get_catalog_info_for_fqs():

    # Extracts event IDs and can cat to textfile for later use or save as an active variable
    evid_sorted = reformat.extract_event_ids(workin_dir=RSQSIM_DIR_RUPT,
                            output_dir=HOME,
                            gns_eventname='event',  #  Prefix used in filenames for RSQSim rupts
    )


    # Create log files for each event ID and drop in fakequakes tree
    reformat.create_log_files(log_file=log_file,
                            txt_file=evid_sorted, 
                            output_dir=FQs_PDIR_RUPT,
    )


    # Create rupture list file for each event ID and drop in fakequakes tree
    reformat.create_rupt_list(event_list=evid_sorted,
                            output_dir=FQs_PROJ_USER + '/data/',
                            output_filename='ruptures.list',)


    # Extract magnitudes (Mw) and locations (latitude, longitude, depth) from catalog csv file for ref
    evid_plus = reformat.process_catalog(catalog_file=catalog_file, 
                            evid_list_sorted=evid_sorted, 
                            output_dir=FQs_PROJ_USER,
                            output_filename='rsqsim_evid_mag_loc.csv',
                            svo=svo_user,
    )


    # Create .rupt files for each event and drop in fakequakes tree
    reformat.process_ruptures(dfile=evid_plus,
                            rsqsim_rupts_dir=RSQSIM_DIR_RUPT,
                            fault=fault_user,
                            mesh=mesho_user, 
                            dx=0.01,
                            rake_avg=90,
                            fqs_dir=FQs_PDIR_RUPT,
                            output_filename='event', 
                            output_plots_dir=FQs_PROJ_USER + '/plots/',
                            svo=svo_user, 
    )



#%%#

# main execution

if __name__ == '__main__':
    
    # Initialize FakeQuakes project structure for testing
    if get_fqs_project_structure:
        init_fakequakes_project_for_testing(home=HOME, project_name=PNM_user)
        print("\nFakeQuakes project structure has been initialized for testing.\n")

    # Extract RSQSim catalog info and convert to FakeQuakes format
    if get_rsqsim_catalog_info_for_fqs:
        run_get_catalog_info_for_fqs()
        print("\nRSQSim catalog information has been processed and saved in the FakeQuakes project structure.\n")




