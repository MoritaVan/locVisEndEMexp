
"""
-----------------------------------------------------------------------------------------
plot_eyetraces_pur_seq.py
-----------------------------------------------------------------------------------------
Goal of the script:
Plot horizontal and vertical eye trace of each run and the fit
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: subject number (sub-01)
-----------------------------------------------------------------------------------------
Output(s):
h5 files with loads of data on eye traces across seqs
-----------------------------------------------------------------------------------------
To run:
cd /Users/martin/Dropbox/Experiments/pMFexp/stats/
python behav_analysis/plot_eye_traces_pur_seq.py sub-01
-----------------------------------------------------------------------------------------
"""

# Stop warnings
# -------------
import warnings
warnings.filterwarnings("ignore")

# General imports
# ---------------
import os
import sys
import platform
import re
import numpy as np
import ipdb
import json
import h5py
import scipy.io
import matplotlib.pyplot as plt
from seaborn import color_palette 
deb = ipdb.set_trace

# Get inputs
# ----------
subject = sys.argv[1]

# Define analysis parameters
# --------------------------
with open('behavior_settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)

# General settings 
# ----------------
if platform.system() == 'Darwin':
    main_dir = analysis_info['main_dir_mac']
    
elif platform.system() == 'Windows':
    main_dir = analysis_info['main_dir_pc']

elif platform.system() == 'Linux':
    main_dir = analysis_info['main_dir_unix']

runs = np.arange(0,analysis_info['num_run'],1)
eye_mov_seq = analysis_info['eye_mov_seq']
seq_trs = analysis_info['seq_trs']
pursuits_tr = np.arange(0,seq_trs,2)
saccades_tr = np.arange(1,seq_trs,2)

# Load data
# ---------
file_dir = '{exp_dir}/data/{sub}'.format(exp_dir = main_dir, sub = subject)
h5_filename = "{file_dir}/add/{sub}_task-PurVELoc_eyedata.h5".format(file_dir = file_dir, sub = subject)
h5_file = h5py.File(h5_filename,'r')
folder_alias = 'eye_traces'
time_start_seq = np.array(h5_file['{folder_alias}/time_start_seq'.format(folder_alias = folder_alias)])
time_end_seq = np.array(h5_file['{folder_alias}/time_end_seq'.format(folder_alias = folder_alias)])
time_start_trial = np.array(h5_file['{folder_alias}/time_start_trial'.format(folder_alias = folder_alias)])
time_end_trial = np.array(h5_file['{folder_alias}/time_end_trial'.format(folder_alias = folder_alias)])
dir_sequence = np.array(h5_file['{folder_alias}/dir_sequence'.format(folder_alias = folder_alias)])
dir_sequence = np.tile(dir_sequence,2)
eye_data = np.array(h5_file['{folder_alias}/eye_data_seqs_nan_blink'.format(folder_alias = folder_alias)])

# Define colors
# -------------
cmap_steps = 16
col_offset = 0#1/14.0
colmap = np.flipud(np.array(color_palette(palette='Reds',n_colors=cmap_steps)))

# Draw figure
# -----------
axis_width = 0.75
col_idx    = 0
# eye trace analysis per run
for run in runs:
    run_eye_data_logic = eye_data[:,3] == run

    if run > 10:run_txt = '{}'.format(run+1)
    else:run_txt = '0{}'.format(run+1)

    # Define figure folder
    try: os.makedirs('{file_dir}/add/figures/PurVELoc/run-{run_txt}'.format(file_dir = file_dir,run_txt = run_txt))
    except: pass

    
    for sequence in np.unique(eye_data[:,4]).astype(int):
        if sequence > 10:seq_txt = '{}'.format(sequence)
        else:seq_txt = '0{}'.format(sequence)

        seq_eye_data_logic = eye_data[:,4] == sequence

        if sequence > 10: sequence_txt = '{}'.format(sequence+1)
        else: sequence_txt = '0{}'.format(sequence+1)
        
        data_logic = np.logical_and.reduce(np.array((run_eye_data_logic,seq_eye_data_logic)))

        time = (eye_data[data_logic][:,0]-eye_data[data_logic][0,0])/1000 # in sec
        
        angleSteps = np.linspace(0,5*2*np.pi,len(time))
        x_pos = 9.6 * np.cos(dir_sequence[sequence-1]*angleSteps)
        y_pos = 9.6 * -np.sin(dir_sequence[sequence-1]*angleSteps)

        mse_x = np.nanmean(np.square(np.subtract(x_pos,eye_data[data_logic,1])))
        mse_y = np.nanmean(np.square(np.subtract(y_pos,eye_data[data_logic,2])))
        
        time_sac = round(np.sum(np.isnan(eye_data[data_logic,2]))/1000,2)

        fig = plt.figure(figsize=(10,15))
        plt.suptitle('Run {run_txt} - Seq {sequence_txt}'.format(run_txt = run+1, sequence_txt = sequence), fontsize = 14)
        plt.subplot(3,2,(1,2))
        plt.plot(time,eye_data[data_logic,1],color = colmap[col_idx],linewidth = axis_width*2)
        plt.plot(time,x_pos,color = 'k',linewidth = axis_width*1)
        plt.ylim((-20,20))
        plt.xlabel('Time (s)')
        plt.ylabel('Position (deg)')
        plt.title('X axis')
        plt.subplot(3,2,(3,4))
        plt.plot(time,eye_data[data_logic,2],color = colmap[col_idx],linewidth = axis_width*2)
        plt.plot(time,y_pos,color = 'k',linewidth = axis_width*1)
        plt.ylim((-20,20))
        plt.xlabel('Time (s)')
        plt.ylabel('Position (deg)')
        plt.title('Y axis')
        plt.subplot(3,2,5)
        plt.plot(eye_data[data_logic,1],eye_data[data_logic,2],color = colmap[col_idx+1],linewidth = axis_width*1.5)
        plt.plot(x_pos,y_pos,color = 'k',linewidth = axis_width*1)
        plt.xlim((-20,20))
        plt.ylim((-20,20))
        plt.xlabel('Position (deg)')
        plt.ylabel('Position (deg)')
        plt.title('Screen view')
        plt.subplot(3,2,6)
        plt.text(.3, .5, 'MSE - x-axis: {mse_x}\nMSE - y-axis: {mse_y}\nSaccades: {time_sac}s of {time_max}s'.format(mse_x = format(mse_x,'.2f'), mse_y = format(mse_y,'.2f'), time_sac = time_sac, time_max = format(np.max(time),'.2f')))
        plt.axis('off')

        plt.savefig("{file_dir}/add/figures/{task}/run-{run_txt}/{sub}_task-{task}_run-{run_txt}_seq-{seq_txt}_eyetraces.png".format(
                                                    sub = subject,
                                                    task = 'PurVELoc',
                                                    run_txt = run_txt,
                                                    seq_txt = seq_txt,
                                                    file_dir = file_dir),facecolor='w')

        col_idx = col_idx + 2
