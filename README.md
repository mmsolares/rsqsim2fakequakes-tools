# rsqsim2fakequakes-tools

Scripts to convert RSQSim rupture outputs into FakeQuakes inputs for synthetic GNSS waveform generation.

## Overview

[FakeQuakes](https://github.com/UO-Geophysics/MudPy/blob/master/src/python/mudpy/fakequakes.py), part of the [MudPy](https://github.com/UO-Geophysics/MudPy) package, simulates both earthquake ruptures and waveforms. However, to use RSQSim-generated ruptures with FakeQuakes, they must first be reformatted into the appropriate input structure.

This repository provides utilities to convert RSQSim rupture outputs into FakeQuakes-compatible formats, enabling synthetic waveform generation.

It is part of a larger framework used to test the G-FAST early earthquake warning algorithm, using RSQSim-generated ruptures for earthquake scenarios in the Hikurangi subduction zone (New Zealand).

More details will be added soon, including references to related article.

## Project Structure
```
rsqsim2fakequakes-tools/
├── data_example/                         # Example datasets for testing and demonstration
│   ├── fqs_input_files/                  # Required FakeQuakes input files
│   │   ├── gns.fault                     # Example fault model file
│   │   ├── gns.mshout                    # Example fault mesh file
│   │   └── wellington8_3.log             # Example FakeQuakes log file
│   ├── rsqsim_catalog_example/           # Sample RSQSim event catalog (e.g., event times, IDs)
│   │   └── ...
│   └── rsqsim_outputs_example/           # Sample RSQSim rupture outputs (e.g., slip distributions)
│       └── ...
│
├── fakequakes_project_example/           # Dummy FakeQuakes project directory (template)
│   └── ...                               # See full structure below
│
├── src/                                  # Source code
│   └── python/
│       └── reforma_rsqsim_catalog.py     # Core tool for converting RSQSim data to FakeQuakes format
│
├── main.py                               # Main example script for running the conversion
├── LICENSE                               # Repository license
└── README.md                             # Project documentation (this file)
```

A typical FakeQuakes project follows a standardized directory structure, such as:

```
fakequakes_project_example/
├── analysis/
│   └── ...
├── data/
│   ├── distances/
│   ├── model_info/
│   ├── station_info/
│   └── ...
├── forward_models/
├── GFs/
│   └── ...
├── logs/
├── output/
│   ├── ruptures/
│   └── ...
├── plots/
│   └── ...
├── scripts/
├── structure/
├── rsqsim_evid_mag_loc.csv   # output from rsqsim2fakequakes-tools, not FakeQuakes module
```


In this repository, you will find a **dummy directory** that mimics a typical FakeQuakes project. This serves as a template for organizing your own project.

## Converting RSQSim Ruptures

RSQSim-generated rupture files must be reformatted to match the MudPy/FakeQuakes rupture format before waveform simulation can proceed. This repository includes scripts to perform that conversion.

After running the `main.py` script, the following files and directories will be generated:

- Reformatted rupture files (`.rupt`) are saved in:  
  `fakequakes_project_example/outputs/ruptures/`

- For each rupture, a dummy log file (`.log`) with the same content but a different header is created in the same directory.

- A rupture list file (`ruptures.list`) is generated in:  
  `fakequakes_project_example/DATA/`

- Rupture summary plots are saved in:  
  `fakequakes_project_example/DATA/`

These outputs follow the expected FakeQuakes directory structure, enabling seamless integration for synthetic waveform generation.

> ⚠️ Note: Additional steps may be necessary depending on how your RSQSim ruptures are structured. Please refer to the provided examples and script documentation.

## How to Run

To convert RSQSim rupture files into the FakeQuakes format, follow these steps:

1. **Clone this repository** and navigate into it:

   ```
   git clone https://github.com/yourusername/rsqsim2fakequakes-tools.git
   cd rsqsim2fakequakes-tools  
   ```
2. Install dependencies:
    
    Make sure you have the following dependencies installed in your Python environment:

    - numpy, matplotlib, scipy, and pyproj
    
    Built-in modules like os and glob are part of Python and require no installation.

3. MudPy Requirement: 

    This script uses get_rise_times from the mudpy.fakequakes module. You must have MudPy installed and configured.
    
    Refer to the [MudPy installation guide](https://github.com/UO-Geophysics/MudPy/wiki) for instructions.

4. Prepare your input:

    - Create a directory to place your RSQSim rupture files
    - You will also need supporting input files from FakeQuakes, including:
        - `gns.fault`
        - `gns.mshout`
        - `wellington8_3.log`
        
    > These files define the fault geometry, subfault mesh, and a sample log file required for formatting the output correctly.

    - Verify and specify the appropriate input paths and variable directly in the `main.py` script.

5. Run the main script:

    You can run the script from the command line in your Python environment (`python main.py`) or within your preferred editor (e.g., VSCode)
    
    This will:

    - Generate csv file (`rsqsim_evid_mag_loc.csv`) with event information (event ID, magnitude, latitude, longitude and depth).
    - Convert RSQSim ruptures to FakeQuakes format
    - Generate .rupt and .log files
    - Write a ruptures.list
    - Create summary plots for each rupture

6. Check outputs in the fakequakes_project_example/ directory: 
    
    - Reformatted rupture files: outputs/ruptures/
    - Dummy log files: outputs/ruptures/
    - Rupture list: data/ruptures.list
    - Plot summaries: plots/

> ✅ Your converted rupture data is now ready for use with FakeQuakes waveform simulations.

### BONUS: Requirements for Waveform Generation

To generate synthetic GNSS waveforms using FakeQuakes, the following components are required:

- **Fault model** – defines the geometry and segmentation of the fault (e.g., my_fakequake_project/data/model_info/XXXX.fault)
- **Station list** – defines the GNSS station locations (e.g., my_fakequake_project/data/station_info/XXXX.gflist)
- **Velocity model** – needed for Green's function computations (e.g., my_fakequake_project/structure/XXXX.mod)
- **Rupture files** – in FakeQuakes format, converted from RSQSim outputs. (e.g., my_fakequake_project/output/ruptures/subd_event_022.rupt)
- **Log files** - dummy log files in FakeQuakes format (e.g., my_fakequake_project/output/ruptures/subd_event_022.log)
