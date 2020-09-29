
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
from ANEMO import ANEMO, read_edf
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
occlusion_data = scipy.io.loadmat('{exp_dir}/data/{sub}/add//{sub}_task_occlusion_size.mat'.format(exp_dir = main_dir, sub = subject))
occlusion_data = occlusion_data['occlusion']['occl_color_mat'][0][[1,3]]
# Define colors
# -------------
cmap_steps = 16
col_offset = 0#1/14.0
colmap = np.flipud(np.array(color_palette(palette='Reds',n_colors=cmap_steps)))

#  ANEMO parameters
screen_width_px  = 1920 # px
screen_height_px = 1080 # px
screen_width_cm  = 70   # cm
viewingDistance  = 108.  # cm

tan              = np.arctan((screen_width_cm/2)/viewingDistance)
screen_width_deg = 2. * tan * 180/np.pi
px_per_deg       = screen_width_px / screen_width_deg

param_exp = {# Mandatory :
                # - number of trials per block :
                      'N_trials' : 1,
                # - number of blocks :
                      'N_blocks' : 1,
                # - direction of the target :
                    # list of lists for each block containing the direction of
                    # the target for each trial is to -1 for left 1 for right
                      'dir_target' : [[],[]], # will be defined in the main loop
                # - number of px per degree for the experiment :
                      'px_per_deg' : 1,
              }

sacc_params = {
        'mindur': 1,
        'maxdur': 100,
        'minsep': 2,
        'before_sacc': 5,
        'after_sacc': 10
}
# default:
#     'mindur': 5,
#     'maxdur': 100,
#     'minsep': 30,
#     'before_sacc': 5,
#     'after_sacc': 15


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

    occl_run = occlusion_data[run][:,1,0:30] # frames x dir x TRs
    new_occl_run = []
    for idx,row in enumerate(occl_run.T):
        new_occl_run.extend(row)
    new_occl_run = np.array(new_occl_run)

    
    
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
        
        
        new_len = np.linspace(0,new_occl_run.shape[0]-1, len(time))
        new_len = np.round(new_len).astype(int)
        new_occl_run = new_occl_run[new_len]
        occlusion = new_occl_run == 0

        # creates an ANEMO instance
        A   = ANEMO(param_exp)

        velocity_deg_x = A.velocity_deg(data_x = eye_data[data_logic,1],
                            filt = 'velocity-position', cutoff = 30, sample_rate = 1000)

        velocity_deg_y = A.velocity_deg(data_x = eye_data[data_logic,2],
                            filt = 'velocity-position', cutoff = 30, sample_rate = 1000)

        misac = A.detec_misac(velocity_x = velocity_deg_x,
                                velocity_y = velocity_deg_y,
                                t_0        = eye_data[data_logic,0][0],
                                VFAC       = 5,
                                mindur     = sacc_params['mindur'],
                                maxdur     = sacc_params['maxdur'],
                                minsep     = sacc_params['minsep'])

        new_saccades = [[x[0]-sacc_params['before_sacc'],x[1]-sacc_params['after_sacc']] for x in misac]

        velocity_x_NAN = velocity_deg_x
        velocity_y_NAN = velocity_deg_y
        for sacc in new_saccades:
            idx = np.logical_and(eye_data[data_logic,0] >= sacc[0], eye_data[data_logic,0] < sacc[1])
            velocity_x_NAN[idx] = np.nan
            velocity_y_NAN[idx] = np.nan
            eye_data[data_logic][idx,1:3] = np.nan

        time_sac = round(np.sum(np.isnan(eye_data[data_logic,2]))/1000,2) # sum of saccades' time

        fig = plt.figure(figsize=(25,10))
        plt.suptitle('Run {run_txt} - Seq {sequence_txt}'.format(run_txt = run+1, sequence_txt = sequence), fontsize = 14)
        plt.subplot(2,5,(1,2))
        plt.fill_between(time, -20, 20, where=occlusion, color='gray', alpha=0.1, interpolate=True)
        plt.plot(time,eye_data[data_logic,1],color = colmap[col_idx],linewidth = axis_width*2)
        plt.plot(time,x_pos,color = 'k',linewidth = axis_width*1)
        plt.ylim((-20,20))
        plt.xlabel('Time (s)')
        plt.ylabel('Position (deg)')
        plt.title('X axis')
        plt.subplot(5,2,(6,7))
        plt.fill_between(time, -20, 20, where=occlusion, color='gray', alpha=0.1, interpolate=True)
        plt.plot(time,eye_data[data_logic,2],color = colmap[col_idx],linewidth = axis_width*2)
        plt.plot(time,y_pos,color = 'k',linewidth = axis_width*1)
        plt.ylim((-20,20))
        plt.xlabel('Time (s)')
        plt.ylabel('Position (deg)')
        plt.title('Y axis')
        plt.subplot(5,2,(3,4))
        plt.fill_between(time, -20, 20, where=occlusion, color='gray', alpha=0.1, interpolate=True)
        plt.plot(time,velocity_x_NAN,color = colmap[col_idx],linewidth = axis_width*2)
        plt.ylim((-20,20))
        plt.xlabel('Time (s)')
        plt.ylabel('Velocity')
        plt.title('X axis')
        plt.subplot(5,2,(8,9))
        plt.fill_between(time, -20, 20, where=occlusion, color='gray', alpha=0.1, interpolate=True)
        plt.plot(time,velocity_y_NAN,color = colmap[col_idx],linewidth = axis_width*2)
        plt.ylim((-20,20))
        plt.xlabel('Time (s)')
        plt.ylabel('Velocity')
        plt.title('Y axis')
        plt.subplot(5,2,5)
        plt.plot(eye_data[data_logic,1],eye_data[data_logic,2],color = colmap[col_idx+1],linewidth = axis_width*1.5)
        plt.plot(x_pos,y_pos,color = 'k',linewidth = axis_width*1)
        plt.xlim((-20,20))
        plt.ylim((-20,20))
        plt.xlabel('Position (deg)')
        plt.ylabel('Position (deg)')
        plt.title('Screen view')
        plt.subplot(5,2,10)
        plt.text(.3, .5, 'MSE - x-axis: {mse_x}\nMSE - y-axis: {mse_y}\nSaccades: {time_sac}s of {time_max}s'.format(mse_x = format(mse_x,'.2f'), mse_y = format(mse_y,'.2f'), time_sac = time_sac, time_max = format(np.max(time),'.2f')))
        plt.axis('off')


        plt.savefig("{file_dir}/add/figures/{task}/run-{run_txt}/{sub}_task-{task}_run-{run_txt}_seq-{seq_txt}_eyetraces.png".format(
                                                    sub = subject,
                                                    task = 'PurVELoc',
                                                    run_txt = run_txt,
                                                    seq_txt = seq_txt,
                                                    file_dir = file_dir),facecolor='w')

        col_idx = col_idx + 1
