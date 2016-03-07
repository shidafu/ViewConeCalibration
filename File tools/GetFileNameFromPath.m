%% Get File Name From Path
% * Author: WANG Lei,USTB
%
% * Link: <https://github.com/shidafu/ViewConeCalibration.git>
%
% * Date:2016/3/8
%
% 
% * Inputs:
%
%     filePath ---- File path
% 
% * Outputs:
%
%     fileName----  File name
%
function fileName= GetFileNameFromPath(filePath)
%Get File Name From Path
    filePath=StrDelTail(filePath);
    srcStrLen=length(filePath);
    fileName=filePath;
    for i=srcStrLen:-1:1
        if strcmp(filePath(1,i),'\') ||strcmp(filePath(1,i),'/')
            fileName=filePath(1,i+1:srcStrLen);
            break;
        end
    end
    filePath=fileName;
    srcStrLen=length(filePath);
    for i=srcStrLen:-1:1
        if strcmp(filePath(1,i),'.')
            fileName=filePath(1,1:i-1);
            break;
        end
    end
end

