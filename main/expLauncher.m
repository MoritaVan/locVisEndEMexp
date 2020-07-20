%% General experiment launcher
%  =============================
% By     :  Martin SZINTE & Vanessa C Morita
% Projet :  Localisers experiment
% With   :  Anna MONTAGNINI & Guillaume MASSON
% Version:  1.0

% Version description
% ===================
% Experiment consists of a localiser
% - visually guided vs endogenous (saccades, pursuit, fixation)



% First settings
% --------------
Screen('CloseAll');clear all;clear mex;clear functions;close all;home;AssertOpenGL;

% General settings
% ----------------
const.expName           =   'VElocEMexp';   % experiment name.
const.expStart          =   1;              % Start of a recording exp                          0 = NO  , 1 = YES
const.checkTrial        =   1;              % Print trial conditions (for debugging)            0 = NO  , 1 = YES
const.writeLogTxt       =   1;              % write a log file in addition to eyelink file      0 = NO  , 1 = YES
const.mkVideo           =   0;              % Make a video of a run (on mac not linux)          0 = NO  , 1 = YES

% External controls
% -----------------
const.tracker           =   1;              % run with eye tracker                              0 = NO  , 1 = YES
const.scanner           =   0;              % run in MRI scanner                                0 = NO  , 1 = YES
const.scannerTest       =   1;              % run with T returned at TR time                    0 = NO  , 1 = YES
const.room              =   2;              % run in MRI or eye-tracking room                   1 = MRI , 2 = eye-tracking

% Run order
% ---------
const.cond_run_order = [1;...               % run 01 - Sac      
                        2;...               % run 02 - Pur      
                        1;...               % run 03 - Sac      
                        2];                 % run 04 - Pur     
% Run number per condition
% ------------------------
const.cond_run_num   = [01;01;...
                        02;02];

% Desired screen setting
% ----------------------
const.desiredFD         =   60;            % Desired refresh rate
fprintf(1,'\n\n\tDon''t forget to change before testing\n');
const.desiredRes        =   [1920,1080];    % Desired resolution

% Path
% ----
dir                     =   (which('expLauncher'));
cd(dir(1:end-18));

% Add Matlab path
% ---------------
addpath('config','main','conversion','eyeTracking','instructions','trials','stim','stats');


% Subject configuration
% ---------------------
[const]                 =   sbjConfig(const);
                        
% Main run
% --------
main(const);










