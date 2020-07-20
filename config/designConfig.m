function [expDes]=designConfig(const)
% ----------------------------------------------------------------------
% [expDes]=designConfig(const)
% ----------------------------------------------------------------------
% Goal of the function :
% Define experimental design
% ----------------------------------------------------------------------
% Input(s) :
% const : struct containing constant configurations
% ----------------------------------------------------------------------
% Output(s):
% expDes : struct containg experimental design
% ----------------------------------------------------------------------
% Function created by Martin SZINTE, modified by Vanessa C Morita
% Project :     locVisEndEMexp
% Version :     1.0
% ----------------------------------------------------------------------

%% Experimental random variables

% Cond 1 : task (1 modality)
% =======
expDes.oneC             =   1;
expDes.txt_cond1        =   {'eyemov'};
% 01 = eye movement


% Var 1 : trial types (3 modalities)
% ======
expDes.oneV             =   [1;2;3];
expDes.txt_var1         =   {'sac','pur','fix'};
% 01 = saccade
% 02 = smooth pursuit
% 03 = fixation

% Var 2 : eye movement direction (3 modalities)
% ======
expDes.twoV            =   [01;02;03;];
expDes.txt_var2         =   {'ccw','cw','none'};
% 01 =   counter clock-wise    
% 02 =   clock-wise
% 03 =   none

% Var 3 : eye movement start position (9 modalities)
% ======
expDes.threeV           =   [01, 02;
                             01, 10;
                             01, 02;
                             01, 10;
                             01, 02;
                             01, 10;
                             01, 02;
                             01, 10;
                             17, 17];
expDes.txt_var3=   {  '0 deg',  '45 deg',  '90 deg', '135 deg', '180 deg', '225 deg', '270 deg', '315 deg',...
                    '180 deg', '225 deg', '270 deg', '315 deg',   '0 deg',  '45 deg',  '90 deg', '135 deg',...
                       'none'};
% pursuit           saccade
% 01 =   0.0 deg    01 =  45.0 deg
% 02 =   0.0 deg    02 = 225.0 deg
% 03 =   0.0 deg    03 = 135.0 deg
% 04 =   0.0 deg    04 = 315.0 deg
% 05 =   0.0 deg    05 = 225.0 deg
% 06 =   0.0 deg    06 =  45.0 deg
% 07 =   0.0 deg    07 = 315.0 deg
% 08 =   0.0 deg    08 = 135.0 deg
% 17 = none  

% Var 4 : trial types 2 (3 modalities)
% ======
expDes.fourV            =   [01;02;03;];
expDes.txt_var4         =   {'vis','end','none'};
% 01 =   visually guided    
% 02 =   endogenous 
% 03 =   none

% Var 5 : eye movement amplitude (2 modalities)
% ======
expDes.fiveV             =   [1;2];
expDes.txt_var5         =   {'19.2 dva','none'};
% 01 = 19.2 dva
% 02 = none

% seq order
% ---------
if const.runNum == 1
    % create sequence order
    dir_sequence.eyemov_val = expDes.twoV(randperm(numel(expDes.twoV)-1));
    
    dir_sequence.val                        =     nan(size(const.eyemov_seq));
    dir_sequence.val(const.eyemov_seq==1)   =     numel(expDes.twoV);    
    dir_sequence.val(const.eyemov_seq==2)   =     repmat(dir_sequence.eyemov_val,[length(dir_sequence.val(const.eyemov_seq==2))/length(dir_sequence.eyemov_val),1]);
    
    expDes.dir_sequence   =   dir_sequence.val;
    
    first_task     = const.cond2;
    txt_first_task = expDes.txt_var1{first_task};
 
    expDes.first_task     = first_task;
    expDes.txt_first_task = txt_first_task;
    
    save(const.task_dir_sequence_file, 'dir_sequence', 'first_task', 'txt_first_task');
else
    load(const.task_dir_sequence_file);
    expDes.dir_sequence   =   dir_sequence.val;
    expDes.first_task     =   first_task;
    expDes.txt_first_task =   txt_first_task;
end
%% Experimental configuration :
expDes.nb_cond          =   1;
expDes.nb_var           =   5;
expDes.nb_rand          =   0;
expDes.nb_list          =   0;

%% Experimental loop
rng('default');rng('shuffle');
runT                    =   const.runNum;

t_trial = 0;
for t_seq = 1:size(const.eyemov_seq,2)
    
    cond1 = const.cond1;
    rand_var2 =   expDes.dir_sequence(t_seq);
    
    if rand_var2 == 3
        seq_steps = const.fix_step;
    else
        seq_steps = const.eyemov_step;
    end
    
    for seq_step = 1:seq_steps
        if rand_var2 == 3
            rand_var1 = expDes.oneV(end);
            rand_var3 = expDes.threeV(end,1);
            rand_var4 = expDes.fourV(end);
            rand_var5 = expDes.fiveV(end);
        else
            rand_var1 = expDes.oneV(const.cond2);
            rand_var5 = expDes.fiveV(1);
            
            if seq_step > 8
                idx = mod(seq_step,8)+1;
            else
                idx = seq_step;
            end
            if const.cond2 == 1 % saccade run
                rand_var3 = expDes.threeV(idx,2);
                if mod(seq_step,4) == 0 || mod(seq_step,4) == 3
                    rand_var4 = expDes.fourV(2);
                else
                    rand_var4 = expDes.fourV(1);
                end
            else
                rand_var3 = expDes.threeV(idx,1);
                if mod(seq_step,2) == 0
                    rand_var4 = expDes.fourV(2);
                else
                    rand_var4 = expDes.fourV(1);
                end
            end
        end
    
        t_trial     =   t_trial + 1;
        
        expDes.expMat(t_trial,:)=   [   runT,           t_trial,        cond1,          rand_var1,      rand_var2,      rand_var3,      ...
                                        rand_var4,      rand_var5,      t_seq,          seq_step,       NaN,            NaN];
        % col 01:   Run number
        % col 02:   Trial number
        % col 03:   Task
        % col 04:   Eye mov type
        % col 05:   Eye mov direction
        % col 06:   Eye mov start position
        % col 07:   Type of mov (vis or end)
        % col 08:   Eye mov amplitude
        % col 09:   Sequence number
        % col 10:   Sequence trial (trial number within a sequence of one mov type)
        % col 11:   Trial onset time
        % col 12:   Trial offset time
    end
end
expDes.nb_trials = size(expDes.expMat,1);


end