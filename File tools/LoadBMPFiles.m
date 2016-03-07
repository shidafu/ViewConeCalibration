%% Load BMP Files from folder
% * Author: WANG Lei,USTB
%
% * Link: <https://github.com/shidafu/ViewConeCalibration.git>
%
% * Date:2016/3/8
%
% 
% * Inputs:
%
%     pathName ---- folder path
% 
% * Outputs:
%
%     fileNum---- Inner and outer camera para matrix
%     filePathArray---- axis vector and angle
%     fileNameArray ---- Object focus length
%
function [fileNum filePathArray fileNameArray] =LoadBMPFiles(pathName)
fileNameArray = ls(strcat(pathName,'\*.bmp'));
lengthStr=length(pathName);
if ~strcmp(pathName(1,lengthStr),'\')
    pathName=[pathName '\'];
end
filePathArray = strcat(pathName,fileNameArray); % Get file path
%filePathSet=cell2mat(filePathSet);
%fileNum = length(size(filePathArray)); % Get file num
fileNum = size(filePathArray,1); % Get file num
clc;