%% Get UV from XYZ according to 3D Homography of a projective geometry transform
% Author:WANG Lei,USTB,email£ºtowanglei@163.com
%
% Date:2016/3/3
%
% Algorithom:
%
% $$\mathit{s} \left(\begin{array}{c} \mathit{u}\\ \mathit{v}\\ 1 \end{array}\right)=
% \left(\begin{array}{cccc} 
% \mathit{h}_{11} & \mathit{h}_{12} & \mathit{h}_{13} & \mathit{h}_{14}\\
% \mathit{h}_{21} & \mathit{h}_{22} & \mathit{h}_{23} & \mathit{h}_{24}\\
% \mathit{h}_{31} & \mathit{h}_{32}  & \mathit{h}_{33} & 1 \end{array}\right)
% \left(\begin{array}{c} \mathit{x}\\ \mathit{y}\\ \mathit{z}\\ 1 \end{array}\right)$$
%
% 

function UV = GetUVFromXYZ(H,XYZ)
% Get UV from XYZ according to 3D Homography of a projective geometry transform
%
% Inputs:
%     H----3 by 3 matrix: [h11 h12 h13 h14;
%                          h21 h22 h23 h24;
%                          h31 h32 h33 1]
%     XY----pointNum by cordNum matrix,
%                  cordNum==3,pointNum>=4,
%                  [x1,y1,z1;
%                   x2,y2,z1;
%                     :]
%
% Outputs:
%     UV----pointNum by cordNum matrix,
%                  cordNum==2,pointNum>=4,
%                  [u1,v1;
%                   u2,v2;
%                     :]
%% Initial
[pointNum, cordNum]=size(XYZ);
if ~(cordNum==3 || cordNum==4)
    error('Input matrix size error!');
end
[hH, wH]=size(H);
if hH~=3 || wH~=4
    error('Input matrix size error!');
end
XYZ1=ones(pointNum,4,'double');
XYZ1(:,1:3)=XYZ(:,1:3);
UV=zeros(pointNum,2,'double');
%% Algorithm.
UVu=H(1,:)*XYZ1'./H(3,:)*XYZ1';
UV(:,1)=UVu';
UVv=H(2,:)*XYZ1'./H(3,:)*XYZ1';
UV(:,2)=UVv';