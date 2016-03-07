%% Load Files from folder of required type
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
%     filetype ---- 'bmp' e.g.
% 
% * Outputs:
%
%     fileNum---- Inner and outer camera para matrix
%     filePathArray---- axis vector and angle
%     fileNameArray ---- Object focus length
%
function [fileNum filePathArray fileNameArray] =LoadFiles(pathName,filetype)
fileNameArray = ls(strcat(pathName,'\*.',filetype));
filePathArray = strcat(pathName,fileNameArray); % Get file path
%filePathSet=cell2mat(filePathSet);
%fileNum = length(size(filePathArray)); % Get file num
fileNum = size(filePathArray,1); % Get file filenum
clc;