%% Get File ext From Path
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
%     fileExt----  File ext
%
function fileExt= GetFileExtFromPath(filePath)
%Get File ext From Path
    filePath=StrDelTail(filePath);
    srcStrLen=length(filePath);
    fileExt=filePath;
    for i=srcStrLen:-1:1
        if strcmp(filePath(1,i),'\') ||strcmp(filePath(1,i),'/')
            fileExt=filePath(1,i+1:srcStrLen);
            break;
        end
    end
    filePath=fileExt;
    srcStrLen=length(filePath);
    for i=srcStrLen:-1:1
        if strcmp(filePath(1,i),'.')
            fileExt=filePath(i,i+1:srcStrLen);
            break;
        end
    end
end

