function [const]=constConfig(scr,const)
% ----------------------------------------------------------------------
% [const]=constConfig(scr,const)
% ----------------------------------------------------------------------
% Goal of the function :
% Define all constant configurations
% ----------------------------------------------------------------------
% Input(s) :
% scr : struct containing screen configurations
% const : struct containing constant configurations
% ----------------------------------------------------------------------
% Output(s):
% const : struct containing constant configurations
% ----------------------------------------------------------------------
% Function created by Martin SZINTE, modified by Vanessa C Morita
% Project :     locVisEndEMexp
% Version :     1.0
% ----------------------------------------------------------------------

% Randomization
rng('default');
rng('shuffle');

%% Colors
const.white             =   [255,255,255];                                                      % white color
const.black             =   [0,0,0];                                                            % black color
const.gray              =   [128,128,128];                                                      % gray color
const.red               =   [200,0,0];                                                          % red
const.background_color  =   const.black;                                                        % background color
const.dot_color         =   const.white;                                                        % define fixation dot color

% Fixation circular aperture
const.fix_out_rim_radVal=   0.3;                                                                % radius of outer circle of fixation bull's eye
const.fix_rim_radVal    =   0.75*const.fix_out_rim_radVal;                                      % radius of intermediate circle of fixation bull's eye in degree
const.fix_radVal        =   0.25*const.fix_out_rim_radVal;                                      % radius of inner circle of fixation bull's eye in degrees
const.fix_out_rim_rad   =   vaDeg2pix(const.fix_out_rim_radVal,scr);                            % radius of outer circle of fixation bull's eye in pixels
const.fix_rim_rad       =   vaDeg2pix(const.fix_rim_radVal,scr);                                % radius of intermediate circle of fixation bull's eye in pixels
const.fix_rad           =   vaDeg2pix(const.fix_radVal,scr);                                    % radius of inner circle of fixation bull's eye in pixels


%% Time parameters
const.TR_dur            =   1.2;                                                                % repetition time
const.TR_num            =   (round(const.TR_dur/scr.frame_duration));                           % repetition time in screen frames

const.eyemov_seq        =   [1,2,1,2,1,2,1,2,1];                                                % 1 = blank/fixation, 2 = eye movement (pursuit or saccade, depending on the run)
const.seq_num           =   numel(const.eyemov_seq);                                            % number of sequences per run

const.eyemov_step       =   32;                                                                 % eye movement steps (in TRs)
const.fix_step          =   16;                                                                 % fixation period step
const.eyemov_step_dur   =   const.TR_dur;                                                       % eye movement steps in seconds
const.eyemov_step_num   =   (round(const.eyemov_step_dur/scr.frame_duration));                  % eye movement step duration in screen frames

const.fix_step_dur      =   const.TR_dur;                                                       % fixation step duration in seconds
const.fix_step_num      =   (round(const.fix_step_dur/scr.frame_duration));                     % fixation step duration in screen frames
const.fix_fixation_dur  =   const.TR_dur;                                                       % visually guided fixation duration in seconds
const.fix_fixation_num  =   (round(const.fix_fixation_dur/scr.frame_duration));                 % visually guided duration in screen frames


const.eyemov_ampVal     =   [19.2];                                                             % eye movement amplitude in visual degrees / diameter of the circunference (to match the paper XXX)
const.eyemov_amp        =   vaDeg2pix(const.eyemov_ampVal,scr);                                 % eye movement amplitude in pixel

const.eyemov_start_step =   90;                                                                 % start position steps in degrees
const.eyemov_direc      =   [1 -1];                                                             % eye movement direction (cw/ccw)
const.purs_steps        =   45;                                                                 % size of each occlusion+pursuit in degrees (4 steps/half circle)
const.occl_start        =   180:const.purs_steps:360-1;                                         % pursuit start position in degrees (for the second half)
const.purs_start        =   const.occl_start + randi([300,375],1,4)/10;                         % pursuit start position in degrees (for the second half)
const.sacc_start        =   45;                                                                 % saccade start position in degrees
const.sacc_positions    =   45:const.eyemov_start_step:360-1;                                   % saccade positions in degrees
% const.purs_occl_start   =   mod(const.purs_start + 270 - const.purs_occl_size/2, 360);          % start position of the occlusion in degrees

const.pursuit_fix_dur   =   0.000;                                                              % first fixation duration in seconds
const.pursuit_fix_num   =   (round(const.pursuit_fix_dur/scr.frame_duration));                  % first fixation duration in screen frames
const.pursuit_dur       =   6*const.TR_dur - const.pursuit_fix_dur;                             % eye movement total duration in seconds (4.8s)
const.pursuit_num       =   (round(const.pursuit_dur/scr.frame_duration));                      % eye movement total duration in screen frames
const.pursuit_end_dur   =   0.000;                                                              % return saccade duration in seconds
const.pursuit_end_num   =   (round(const.pursuit_end_dur/scr.frame_duration));                  % return saccade duration in screen frames
const.pursuit_ang_step  =   360/const.pursuit_num;                                              % eye movement angle step

const.saccade_fix_dur   =   0.400;                                                              % first fixation duration in seconds
const.saccade_fix_num   =   (round(const.saccade_fix_dur/scr.frame_duration));                  % first fixation duration in screen frames
const.saccade_tot_dur   =   0.800;                                                              % eye movement total duration in seconds
const.saccade_tot_num   =   (round(const.saccade_tot_dur/scr.frame_duration));                  % eye movement total duration in screen frames

% define TR for scanner
if const.scanner
    const.TRs = 0;
    for seq = const.eyemov_seq
        if seq == 1
            TR_seq = const.fix_step;
        else
            TR_seq = const.eyemov_step;
        end
        const.TRs = const.TRs + TR_seq;
    end
    const.TRs = const.TRs;
    fprintf(1,'\n\tScanner parameters; %1.0f TRs, %1.2f seconds, %s\n',const.TRs,const.TR_dur,datestr(seconds((const.TRs*const.TR_dur)),'MM:SS'));
end

% compute fixation coordinates
const.fixation_matX(1:const.fix_fixation_num)                    = scr.x_mid;
const.fixation_matY(1:const.fix_fixation_num)                    = scr.y_mid;
const.fixation_matX(const.fix_fixation_num+1:const.fix_step_num) = -scr.x_mid;
const.fixation_matY(const.fix_fixation_num+1:const.fix_step_num) = -scr.y_mid;

% compute pursuit coordinates
% 6 TR = 1 circunference pursuit mov
for eyemov_direc = 1:size(const.eyemov_direc,2)  
    % compute a 6-TR mov
    for purs_diameter = 1:size(const.eyemov_amp,2)
        % fixation
        if const.pursuit_fix_num > 0
            step1 = 1:const.pursuit_fix_num;
            pursuit_matX(step1,eyemov_direc) = scr.x_mid + const.eyemov_amp(purs_diameter)/2 * cosd(0);
            pursuit_matY(step1,eyemov_direc) = scr.y_mid + const.eyemov_amp(purs_diameter)/2 * (-sind(0));
        else 
            step1 = 0;
        end

        % eye movement step
        step2 = (step1(end) + 1):(step1(end) + const.pursuit_num);
        for nbf = step2
            angle = (nbf-1)*const.pursuit_ang_step;

            if (angle >= const.occl_start(1) && angle < const.purs_start(1)) || (angle >= const.occl_start(2) && angle < const.purs_start(2)) || (angle >= const.occl_start(3) && angle < const.purs_start(3)) || (angle >= const.occl_start(4) && angle < const.purs_start(4))% occluded arc
                pursuit_matX(nbf,eyemov_direc) = -scr.x_mid;
                pursuit_matY(nbf,eyemov_direc) = -scr.y_mid;
            else
                pursuit_matX(nbf,eyemov_direc) = scr.x_mid + const.eyemov_amp(purs_diameter)/2 * cosd(360 + const.eyemov_direc(eyemov_direc)*angle);
                pursuit_matY(nbf,eyemov_direc) = scr.y_mid + const.eyemov_amp(purs_diameter)/2 * (-sind(360 + const.eyemov_direc(eyemov_direc)*angle));
            end
        end

        % fixation
        if const.pursuit_end_num > 0
            step3 = (step2(end) + 1):(step2(end) + const.pursuit_end_num);
            pursuit_matX(step3,eyemov_direc) = pursuit_matX(nbf-1,eyemov_direc);
            pursuit_matY(step3,eyemov_direc) = pursuit_matY(nbf-1,eyemov_direc);
        end
    end
end

% repeat movement nRep times / split into 6-TR sequences
nRep = floor(const.eyemov_step/6);
q1 = (1:nRep)*6-5;
q2 = (1:nRep)*6-4;
q3 = (1:nRep)*6-3;
q4 = (1:nRep)*6-2;
q5 = (1:nRep)*6-1;
q6 = (1:nRep)*6-0;
for eyemov_direc = 1:size(const.eyemov_direc,2)
    const.pursuit_matX(:, eyemov_direc, q1) = repmat(pursuit_matX(0*const.eyemov_step_num+1:1*const.eyemov_step_num,eyemov_direc),1,nRep);
    const.pursuit_matX(:, eyemov_direc, q2) = repmat(pursuit_matX(1*const.eyemov_step_num+1:2*const.eyemov_step_num,eyemov_direc),1,nRep);
    const.pursuit_matX(:, eyemov_direc, q3) = repmat(pursuit_matX(2*const.eyemov_step_num+1:3*const.eyemov_step_num,eyemov_direc),1,nRep);
    const.pursuit_matX(:, eyemov_direc, q4) = repmat(pursuit_matX(3*const.eyemov_step_num+1:4*const.eyemov_step_num,eyemov_direc),1,nRep);
    const.pursuit_matX(:, eyemov_direc, q5) = repmat(pursuit_matX(4*const.eyemov_step_num+1:5*const.eyemov_step_num,eyemov_direc),1,nRep);
    const.pursuit_matX(:, eyemov_direc, q6) = repmat(pursuit_matX(5*const.eyemov_step_num+1:6*const.eyemov_step_num,eyemov_direc),1,nRep);

    const.pursuit_matY(:, eyemov_direc, q1) = repmat(pursuit_matY(0*const.eyemov_step_num+1:1*const.eyemov_step_num,eyemov_direc),1,nRep);
    const.pursuit_matY(:, eyemov_direc, q2) = repmat(pursuit_matY(1*const.eyemov_step_num+1:2*const.eyemov_step_num,eyemov_direc),1,nRep);
    const.pursuit_matY(:, eyemov_direc, q3) = repmat(pursuit_matY(2*const.eyemov_step_num+1:3*const.eyemov_step_num,eyemov_direc),1,nRep);
    const.pursuit_matY(:, eyemov_direc, q4) = repmat(pursuit_matY(3*const.eyemov_step_num+1:4*const.eyemov_step_num,eyemov_direc),1,nRep);
    const.pursuit_matY(:, eyemov_direc, q5) = repmat(pursuit_matY(4*const.eyemov_step_num+1:5*const.eyemov_step_num,eyemov_direc),1,nRep);
    const.pursuit_matY(:, eyemov_direc, q6) = repmat(pursuit_matY(5*const.eyemov_step_num+1:6*const.eyemov_step_num,eyemov_direc),1,nRep);
    
    const.pursuit_matX(:, eyemov_direc, 31:32) = scr.x_mid;
    const.pursuit_matY(:, eyemov_direc, 31:32) = scr.y_mid;
end

% compute saccade coordinates
% 4 TR = 4 saccades
step1 = 1:const.saccade_fix_num;                                    % fixation 1
step2 = (step1(end) + 1):(step1(end) + const.saccade_tot_num);      % saccade 1
% step3 = step2(end) + step1;                                         % fixation 2
% step4 = step2(end) + step2;                                         % saccade 2

for eyemov_direc = 1:size(const.eyemov_direc,2)
    for sacc_start = 1:size(const.sacc_start,2)
        idx1 = (sacc_start-1)*8 + 1; % TR 1 - vis
        idx2 = (sacc_start-1)*8 + 2; % TR 2 - vis
        idx3 = (sacc_start-1)*8 + 3; % TR 3 - vis
        idx4 = (sacc_start-1)*8 + 4; % TR 4 - vis
        idx5 = idx4+1:idx4+4;        % TR 5-8 - end
        
        sacc(1) = mod(360 + const.eyemov_direc(eyemov_direc)*const.sacc_positions(sacc_start),360);
        sacc(2) = mod(360 + const.eyemov_direc(eyemov_direc)*const.sacc_positions(sacc_start+1),360);
        sacc(3) = mod(360 + const.eyemov_direc(eyemov_direc)*const.sacc_positions(sacc_start+2),360);
        sacc(4) = mod(360 + const.eyemov_direc(eyemov_direc)*const.sacc_positions(sacc_start+3),360);
        if const.eyemov_direc(eyemov_direc) == -1
            sacc = circshift(sacc,1,2);
        end
        
        % TR 1
        const.saccade_matX(step1,eyemov_direc,idx1) = scr.x_mid + (cosd(sacc(1)) * const.eyemov_amp/2);
        const.saccade_matY(step1,eyemov_direc,idx1) = scr.y_mid + (-sind(sacc(1)) * const.eyemov_amp/2);

        const.saccade_matX(step2,eyemov_direc,idx1) = scr.x_mid + (cosd(sacc(2)) * const.eyemov_amp/2);
        const.saccade_matY(step2,eyemov_direc,idx1) = scr.y_mid + (-sind(sacc(2)) * const.eyemov_amp/2);

        % TR 2
        const.saccade_matX(step1,eyemov_direc,idx2) = scr.x_mid + (cosd(sacc(2)) * const.eyemov_amp/2);
        const.saccade_matY(step1,eyemov_direc,idx2) = scr.y_mid + (-sind(sacc(2)) * const.eyemov_amp/2);

        const.saccade_matX(step2,eyemov_direc,idx2) = scr.x_mid + (cosd(sacc(3)) * const.eyemov_amp/2);
        const.saccade_matY(step2,eyemov_direc,idx2) = scr.y_mid + (-sind(sacc(3)) * const.eyemov_amp/2); 
        
        % TR 3
        const.saccade_matX(step1,eyemov_direc,idx3) = scr.x_mid + (cosd(sacc(3)) * const.eyemov_amp/2);
        const.saccade_matY(step1,eyemov_direc,idx3) = scr.y_mid + (-sind(sacc(3)) * const.eyemov_amp/2);

        const.saccade_matX(step2,eyemov_direc,idx3) = scr.x_mid + (cosd(sacc(4)) * const.eyemov_amp/2);
        const.saccade_matY(step2,eyemov_direc,idx3) = scr.y_mid + (-sind(sacc(4)) * const.eyemov_amp/2);

        % TR 4
        const.saccade_matX(step1,eyemov_direc,idx4) = scr.x_mid + (cosd(sacc(4)) * const.eyemov_amp/2);
        const.saccade_matY(step1,eyemov_direc,idx4) = scr.y_mid + (-sind(sacc(4)) * const.eyemov_amp/2);

        const.saccade_matX(step2,eyemov_direc,idx4) = scr.x_mid + (cosd(sacc(1)) * const.eyemov_amp/2);
        const.saccade_matY(step2,eyemov_direc,idx4) = scr.y_mid + (-sind(sacc(1)) * const.eyemov_amp/2);
        
        % TR 5-8
        const.saccade_matX(:,eyemov_direc,idx5) = -scr.x_mid;
        const.saccade_matY(:,eyemov_direc,idx5) = -scr.y_mid;
    end
    
    step = 8 * size(const.sacc_start,2); % 8 TRs
    for idx = step+1:step:32
        const.saccade_matX(:,eyemov_direc,idx:idx+step-1) = const.saccade_matX(:,eyemov_direc,1:step);
        const.saccade_matY(:,eyemov_direc,idx:idx+step-1) = const.saccade_matY(:,eyemov_direc,1:step);
    end
end


%% Eyelink calibration value
const.ppd               =   vaDeg2pix(1,scr);                                                  % get one pixel per degree
const.maxX              =   scr.scr_sizeX*0.5;                                                 % maximum horizontal amplitude of the screen
const.maxY              =   scr.scr_sizeY*0.5;                                                 % maximum vertical amplitude of the screen
const.calib_maxX     	=   const.maxX/2;
const.calib_maxY        =   const.maxY/2;
const.calib_center      =   [scr.scr_sizeX/2,scr.scr_sizeY/2];

const.calibCoord        =   round([ const.calib_center(1),                     const.calib_center(2),...                       % 01.  center center
                                    const.calib_center(1),                     const.calib_center(2)-const.calib_maxY,...      % 02.  center up
                                    const.calib_center(1),                     const.calib_center(2)+const.calib_maxY,...      % 03.  center down
                                    const.calib_center(1)-const.calib_maxX,    const.calib_center(2),...                       % 04.  left center
                                    const.calib_center(1)+const.calib_maxX,    const.calib_center(2),...                       % 05.  right center
                                    const.calib_center(1)-const.calib_maxX,    const.calib_center(2)-const.calib_maxY,...      % 06.  left up
                                    const.calib_center(1)+const.calib_maxX,    const.calib_center(2)-const.calib_maxY,...      % 07.  right up
                                    const.calib_center(1)-const.calib_maxX,    const.calib_center(2)+const.calib_maxY,...      % 08.  left down
                                    const.calib_center(1)+const.calib_maxX,    const.calib_center(2)+const.calib_maxY,...      % 09.  right down
                                    const.calib_center(1)-const.calib_maxX/2,  const.calib_center(2)-const.calib_maxY/2,...    % 10.  mid left mid up
                                    const.calib_center(1)+const.calib_maxX/2,  const.calib_center(2)-const.calib_maxY/2,...    % 11.  mid right mid up
                                    const.calib_center(1)-const.calib_maxX/2,  const.calib_center(2)+const.calib_maxY/2,...    % 12.  mid left mid down
                                    const.calib_center(1)+const.calib_maxX/2,  const.calib_center(2)+const.calib_maxY/2]);     % 13.  mid right mid down

const.valid_maxX        =   const.calib_maxX * 0.9;
const.valid_maxY        =   const.calib_maxY * 0.9;
const.valid_center      =   const.calib_center;

const.validCoord    	=   round([ const.valid_center(1),                     const.valid_center(2),...                       % 01.  center center
                                    const.valid_center(1),                     const.valid_center(2)-const.valid_maxY,...      % 02.  center up
                                    const.valid_center(1),                     const.valid_center(2)+const.valid_maxY,...      % 03.  center down
                                    const.valid_center(1)-const.valid_maxX,    const.valid_center(2),...                       % 04.  left center
                                    const.valid_center(1)+const.valid_maxX,    const.valid_center(2),...                       % 05.  right center
                                    const.valid_center(1)-const.valid_maxX,    const.valid_center(2)-const.valid_maxY,...      % 06.  left up
                                    const.valid_center(1)+const.valid_maxX,    const.valid_center(2)-const.valid_maxY,...      % 07.  right up
                                    const.valid_center(1)-const.valid_maxX,    const.valid_center(2)+const.valid_maxY,...      % 08.  left down
                                    const.valid_center(1)+const.valid_maxX,    const.valid_center(2)+const.valid_maxY,...      % 09.  right down
                                    const.valid_center(1)-const.valid_maxX/2,  const.valid_center(2)-const.valid_maxY/2,...    % 10.  mid left mid up
                                    const.valid_center(1)+const.valid_maxX/2,  const.valid_center(2)-const.valid_maxY/2,...    % 11.  mid right mid up
                                    const.valid_center(1)-const.valid_maxX/2,  const.valid_center(2)+const.valid_maxY/2,...    % 12.  mid left mid down
                                    const.valid_center(1)+const.valid_maxX/2,  const.valid_center(2)+const.valid_maxY/2]);     % 13.  mid right mid down

end