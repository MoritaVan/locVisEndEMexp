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

const.eyemov_seq        =   [1,2,1,2,1,2,1,2, ...
                             1,2,1,2,1,2,1,2, ...
                             1,2,1,2,1,2,1,2, ...
                             1,2,1,2,1,2,1,2,1];                                                % 1 = blank/fixation, 2 = eye movement (pursuit or saccade, depending on the run)
const.seq_num           =   numel(const.eyemov_seq);                                            % number of sequences per run

% const.eyemov_step       =   32;                                                                 % eye movement steps (possible directions)
% const.fix_step          =   16;                                                                 % fixation period step
const.eyemov_step       =   8;                                                                 % eye movement steps (possible directions)
const.fix_step          =   4;                                                                 % fixation period step
const.eyemov_step_dur   =   const.TR_dur;                                                       % eye movement steps in seconds
const.eyemov_step_num   =   (round(const.eyemov_step_dur/scr.frame_duration));                  % eye movement step duration in screen frames

const.fix_step_dur      =   const.TR_dur;                                                       % fixation step duration in seconds
const.fix_step_num      =   (round(const.fix_step_dur/scr.frame_duration));                     % fixation step duration in screen frames
const.fix_fixation_dur  =   0.050;                                                              % visually guided fixation duration in seconds
const.fix_fixation_num  =   (round(const.fix_fixation_dur/scr.frame_duration));                 % visually guided duration in screen frames


const.eyemov_ampVal     =   [2.5, 5, 7.5, 10];                                                  % eye movement amplitude in visual degrees / diameter of the circunference
const.eyemov_amp        =   vaDeg2pix(const.eyemov_ampVal,scr);                                 % eye movement amplitude in pixel

const.eyemov_start_step =   90;                                                                 % start position steps in degrees
const.purs_start        =   0:const.eyemov_start_step:360-1;                                    % pursuit start position in degrees
const.sacc_start        =   45;                                                                 % saccade start position in degrees
const.sacc_positions    =   45:const.eyemov_start_step:360-1;                                   % saccade positions in degrees
const.purs_occl_size    =   90;                                                                 % size of the occlusion in degrees
const.purs_occl_start   =   mod(const.purs_start + 270 - const.purs_occl_size/2, 360);          % start position of the occlusion in degrees

const.pursuit_fix_dur   =   0.000;                                                              % first fixation duration in seconds
const.pursuit_fix_num   =   (round(const.pursuit_fix_dur/scr.frame_duration));                  % first fixation duration in screen frames
const.pursuit_dur       =   2.400 - const.pursuit_fix_dur;                                      % eye movement total duration in seconds
const.pursuit_num       =   (round(const.pursuit_dur/scr.frame_duration));                      % eye movement total duration in screen frames
const.pursuit_end_dur   =   0.000;                                                              % return saccade duration in seconds
const.pursuit_end_num   =   (round(const.pursuit_end_dur/scr.frame_duration));                  % return saccade duration in screen frames
const.pursuit_ang_step  =   360/const.pursuit_num;                                              % eye movement angle step

const.saccade_fix_dur   =   0.200;                                                              % first fixation duration in seconds
const.saccade_fix_num   =   (round(const.saccade_fix_dur/scr.frame_duration));                  % first fixation duration in screen frames
const.saccade_tot_dur   =   0.400;                                                              % eye movement total duration in seconds
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
% 2 TR = 1 circunference pursuit mov
for purs_diameter = 1:size(const.eyemov_amp,2)
    % compute a 2-TR mov
    for purs_start = 1:size(const.purs_start,2)  
        % fixation
        if const.pursuit_fix_num > 0
            step1 = 1:const.pursuit_fix_num;
            pursuit_matX(step1,purs_diameter,purs_start) = scr.x_mid + const.eyemov_amp(purs_diameter)/2 * cosd(const.purs_start(purs_start));
            pursuit_matY(step1,purs_diameter,purs_start) = scr.y_mid + const.eyemov_amp(purs_diameter)/2 * (-sind(const.purs_start(purs_start)));
        else 
            step1 = 0;
        end
        
        % eye movement step
        step2 = (step1(end) + 1):(step1(end) + const.pursuit_num);
        for nbf = step2
            angle = mod(const.purs_start(purs_start) + (nbf-1)*const.pursuit_ang_step, 360);
            
            if angle >= const.purs_occl_start(purs_start) && angle < const.purs_occl_start(purs_start) + const.purs_occl_size % occluded arc
                pursuit_matX(nbf,purs_diameter,purs_start) = -scr.x_mid;
                pursuit_matY(nbf,purs_diameter,purs_start) = -scr.y_mid;
            else
                pursuit_matX(nbf,purs_diameter,purs_start) = scr.x_mid + const.eyemov_amp(purs_diameter)/2 * cosd(angle);
                pursuit_matY(nbf,purs_diameter,purs_start) = scr.y_mid + const.eyemov_amp(purs_diameter)/2 * (-sind(angle));
            end
        end
        
        % fixation
        if const.pursuit_end_num > 0
            step3 = (step2(end) + 1):(step2(end) + const.pursuit_end_num);
            pursuit_matX(step3,purs_diameter,purs_start) = pursuit_matX(nbf-1,purs_diameter,purs_start);
            pursuit_matY(step3,purs_diameter,purs_start) = pursuit_matY(nbf-1,purs_diameter,purs_start);
        end

    end
    
    % repeat movement nRep times and split into 1-TR sequences
    nRep = (const.eyemov_step/2) / length(const.purs_start);
    half1 = 1:2:nRep;
    half2 = (1:2:nRep) + 1;
    for start = 1:size(const.purs_start,2)
        const.pursuit_matX(:, purs_diameter, start, half1) = repmat(pursuit_matX(1:const.eyemov_step_num,purs_diameter,start),1,nRep);
        const.pursuit_matX(:, purs_diameter, start, half2) = repmat(pursuit_matX(const.eyemov_step_num+1:end,purs_diameter,start),1,nRep);
        const.pursuit_matY(:, purs_diameter, start, half1) = repmat(pursuit_matY(1:const.eyemov_step_num,purs_diameter,start),1,nRep);
        const.pursuit_matY(:, purs_diameter, start, half2) = repmat(pursuit_matY(const.eyemov_step_num+1:end,purs_diameter,start),1,nRep);
    end

end

% compute saccade coordinates
% 2 TR = 4 saccades
step1 = 1:const.saccade_fix_num;                                    % fixation 1
step2 = (step1(end) + 1):(step1(end) + const.saccade_tot_num);      % saccade 1
step3 = step2(end) + step1;                                         % fixation 2
step4 = step2(end) + step2;                                         % saccade 2

for sacc_amp = 1:size(const.eyemov_amp,2)
    for sacc_start = 1:size(const.sacc_start,2)
        idx1 = (sacc_start-1)*4 + 1; % TR 1 - vis
        idx2 = (sacc_start-1)*4 + 2; % TR 2 - vis
        idx3 = (sacc_start-1)*4 + 3; % TR 3 - end
        idx4 = (sacc_start-1)*4 + 4; % TR 4 - end
        
        sacc1 = const.sacc_positions(sacc_start);
        if sacc_start+1 <= 4, sacc2 = const.sacc_positions(sacc_start+1); else sacc2 = const.sacc_positions(mod(sacc_start+1,4)); end
        if sacc_start+2 <= 4, sacc3 = const.sacc_positions(sacc_start+2); else sacc3 = const.sacc_positions(mod(sacc_start+2,4)); end
        if sacc_start+3 <= 4, sacc4 = const.sacc_positions(sacc_start+3); else sacc4 = const.sacc_positions(mod(sacc_start+3,4)); end
        
        % TR 1
        const.saccade_matX(step1,sacc_amp,idx1) = scr.x_mid + (cosd(sacc1) * const.eyemov_amp(sacc_amp)/2);
        const.saccade_matY(step1,sacc_amp,idx1) = scr.y_mid + (-sind(sacc1) * const.eyemov_amp(sacc_amp)/2);

        const.saccade_matX(step2,sacc_amp,idx1) = scr.x_mid + (cosd(sacc2) * const.eyemov_amp(sacc_amp)/2);
        const.saccade_matY(step2,sacc_amp,idx1) = scr.y_mid + (-sind(sacc2) * const.eyemov_amp(sacc_amp)/2);

        const.saccade_matX(step3,sacc_amp,idx1) = scr.x_mid + (cosd(sacc2) * const.eyemov_amp(sacc_amp)/2);
        const.saccade_matY(step3,sacc_amp,idx1) = scr.y_mid + (-sind(sacc2) * const.eyemov_amp(sacc_amp)/2);

        const.saccade_matX(step4,sacc_amp,idx1) = scr.x_mid + (cosd(sacc3) * const.eyemov_amp(sacc_amp)/2);
        const.saccade_matY(step4,sacc_amp,idx1) = scr.y_mid + (-sind(sacc3) * const.eyemov_amp(sacc_amp)/2); 
        
        % TR 2
        const.saccade_matX(step1,sacc_amp,idx2) = scr.x_mid + (cosd(sacc3) * const.eyemov_amp(sacc_amp)/2);
        const.saccade_matY(step1,sacc_amp,idx2) = scr.y_mid + (-sind(sacc3) * const.eyemov_amp(sacc_amp)/2);

        const.saccade_matX(step2,sacc_amp,idx2) = scr.x_mid + (cosd(sacc4) * const.eyemov_amp(sacc_amp)/2);
        const.saccade_matY(step2,sacc_amp,idx2) = scr.y_mid + (-sind(sacc4) * const.eyemov_amp(sacc_amp)/2);

        const.saccade_matX(step3,sacc_amp,idx2) = scr.x_mid + (cosd(sacc4) * const.eyemov_amp(sacc_amp)/2);
        const.saccade_matY(step3,sacc_amp,idx2) = scr.y_mid + (-sind(sacc4) * const.eyemov_amp(sacc_amp)/2);

        const.saccade_matX(step4,sacc_amp,idx2) = scr.x_mid + (cosd(sacc1) * const.eyemov_amp(sacc_amp)/2);
        const.saccade_matY(step4,sacc_amp,idx2) = scr.y_mid + (-sind(sacc1) * const.eyemov_amp(sacc_amp)/2);
        
        % TR 3
        const.saccade_matX(:,sacc_amp,idx3) = -scr.x_mid;
        const.saccade_matY(:,sacc_amp,idx3) = -scr.y_mid;
        % TR 4
        const.saccade_matX(:,sacc_amp,idx4) = -scr.x_mid;
        const.saccade_matY(:,sacc_amp,idx4) = -scr.y_mid;
    end
    
    step = 4 * size(const.sacc_start,2); % 4 TRs
    for idx = step+1:step:32
        const.saccade_matX(:,sacc_amp,idx:idx+step-1) = const.saccade_matX(:,sacc_amp,1:step);
        const.saccade_matY(:,sacc_amp,idx:idx+step-1) = const.saccade_matY(:,sacc_amp,1:step);
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