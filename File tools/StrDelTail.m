%% Delete string tail after the '\0' flag.
% * Author: WANG Lei,USTB
%
% * Link: <https://github.com/shidafu/ViewConeCalibration.git>
%
% * Date:2016/3/8
%
% 
% * Inputs:
%
%     srcStr ---- source string
% 
% * Outputs:
%
%     dstStr---- dest string
%
function dstStr= StrDelTail(srcStr)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    allStrLen=length(srcStr);
    srcStrLen=allStrLen;
    for i=allStrLen:-1:1
        if srcStr(1,srcStrLen)==' '
            srcStrLen=srcStrLen-1;
        else
            break;
        end
    end
    if srcStrLen<0
        srcStrLen=0;
    end
    dstStr=srcStr(1,1:srcStrLen);
end

