This package contains the MATLAB implementation and representative NG_PoorE Set1 EIS data for the DT-P2D MSCEF example.

1. DATA FORMAT

Each CSV file must contain three numeric columns without headers:

Column 1: Frequency [Hz]
Column 2: Real(Z) [Ohm]
Column 3: Imag(Z) [Ohm]

The third column must follow the Im(Z) convention used by the MATLAB code. Do not convert it to -Im(Z) before loading.

For another dataset, use the following naming rules:

NG_OptE_Anode_SOC20.csv
NG_OptE_Cathode_SOC20.csv
BG_OptE_Anode_SOC20.csv
BG_OptE_Cathode_SOC20.csv
NG_PoorE_Set1_Anode_SOC20.csv
NG_PoorE_Set1_Cathode_SOC20.csv


2. CODE USAGE AND RESULT CHECK

Run:

run_DT_P2D_MSCEF_example

The default settings are:

root_folder = "AUTO";
save_folder = "AUTO";

DOE = 'NG_PoorE';
set_name = 'Set1';
multi_soc_anode = [20 60];
multi_soc_cathode = [20 60];

num_iter = 100;
num_iter_dist = 30;
start_points = 5;
use_parallel = 1;
use_relative_weight = true;
save_results = true;

When root_folder is set to "AUTO", the code reads the CSV files from the folder containing the main script.

When save_folder is set to "AUTO", the code automatically creates:

Save_example

in the same folder and saves all results there.

To use another data or output folder, replace "AUTO" with the desired folder path.

Example:

root_folder = "D:\DT_P2D_Data";
save_folder = "D:\DT_P2D_Results";

The workflow includes:

1. Loading the selected anode and cathode EIS datasets
2. Preliminary P2D fitting with multiple starting points
3. DT-P2D model calculation using DRT + DDT
4. Electrode-specific impedance element separation
5. Saving the summary table, MATLAB workspace, and SOC-specific figures

The default use_parallel value is 1. Parallel Computing Toolbox is required for parallel multi-start execution. Set use_parallel = 0 when parallel execution is unavailable.

The summary table reports the fitted DT-P2D parameters, separated resistance contributions, and RMSE values for all selected anode and cathode SOCs.

For a quick result check, first open:

1. NG_PoorE_Set1_summary_table.csv
2. The Anode and Cathode SOC-specific .fig files
3. NG_PoorE_Set1_DT_P2D_workspace.mat when detailed MATLAB variables are needed
