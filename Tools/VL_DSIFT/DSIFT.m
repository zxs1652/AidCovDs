function [DSIFT] = DSIFT(Img,Size,Step)
%% Initialization
magnif = 3 ;
%% Operation
Img = vl_imsmooth(single(Img),sqrt((Size/magnif)^2-.25));
[~,DSIFT] = vl_dsift(single(Img),'size',Size,'step',Step);%,'FloatDescriptors'); 
%% Finalization
end

