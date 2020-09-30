function drawEmptyTarget(scr,const,targetX,targetY,color)
% ----------------------------------------------------------------------
% drawEmptyTarget(scr,const,targetX,targetY,color)
% ----------------------------------------------------------------------
% Goal of the function :
% Draw bull's eye target with empty center
% ----------------------------------------------------------------------
% Input(s) :
% scr : struct containing screen configurations
% const : struct containing constant configurations
% targetX: target coordinate X
% targetY: target coordinate Y
% color: bull's eye color
% ----------------------------------------------------------------------
% Output(s):
% none
% ----------------------------------------------------------------------
% Function created by Martin SZINTE (martin.szinte@gmail.com)
% Last update : 07 / 08 / 2020
% Project :     pMFexp
% Version :     1.0
% ----------------------------------------------------------------------


Screen('DrawDots',scr.main,[targetX,targetY],const.fix_out_rim_rad*2, color , [], 2);
Screen('DrawDots',scr.main,[targetX,targetY],const.fix_rim_rad*2, const.background_color, [], 2);
    

end