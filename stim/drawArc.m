function drawArc(scr,const,ang,siz)
% ----------------------------------------------------------------------
% drawTarget(scr,const,targetX,targetY)
% ----------------------------------------------------------------------
% Goal of the function :
% Draw occlusion arcs
% ----------------------------------------------------------------------
% Input(s) :
% scr : struct containing screen configurations
% const : struct containing constant configurations
% ang: start angle
% siz: size of the arc
% ----------------------------------------------------------------------
% Output(s):
% none

rect = [
        floor(scr.x_mid - (const.eyemov_amp/2 + 2*const.fix_out_rim_rad)) 
        floor(scr.y_mid - (const.eyemov_amp/2 + 2*const.fix_out_rim_rad))
        ceil(scr.x_mid + (const.eyemov_amp/2 + 2*const.fix_out_rim_rad))
        ceil(scr.y_mid + (const.eyemov_amp/2 + 2*const.fix_out_rim_rad))
       ];

Screen('FillArc',scr.main,[0],rect,ang,siz);

end