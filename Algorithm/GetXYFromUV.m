%% Get XY from UV according to 2D Homography of a projective geometry transform
% Author:WANG Lei,USTB,email£ºtowanglei@163.com
%
% Date:2016/3/3
%
% Algorithom:
%
% $$\mathit{s} \left(\begin{array}{c} \mathit{x}\\ \mathit{y} \end{array}\right)=
% \left(\begin{array}{ccc} 
% \mathit{h}_{11} & \mathit{h}_{12} & \mathit{h}_{13}\\
% \mathit{h}_{21} & \mathit{h}_{22} & \mathit{h}_{23}\\
% \mathit{h}_{31} & \mathit{h}_{32} & 1 \end{array}\right)
% \left(\begin{array}{c} \mathit{x}\\ \mathit{y}\\ 1 \end{array}\right)$$
%

function XY = GetXYFromUV(H,UV)
% Get XY from UV according to 2D Homography of a projective geometry transform
%
% Inputs:
%     H----3 by 3 matrix: [h11 h12 h13;
%                          h21 h22 h23;
%                          h31 h32  1]
%     UV----pointNum by cordNum matrix,
%                  cordNum==2,pointNum>=4,
%                  [u1,v1;
%                   u2,v2;
%                     :]
%
% Outputs:
%     XY----pointNum by cordNum matrix,
%                  cordNum==2,pointNum>=4,
%                  [x1,y1;
%                   x2,y2;
%                     :]

%% Initial
[pointNum, cordNum]=size(XY);
if ~(cordNum==2 || cordNum==3)
    error('Input matrix size error!');
end
[hH, wH]=size(H);
if hH~=3 ||~(wH==3 || wH==4)
    error('Input matrix size error!');
end
if wH==4
    H=[H(:,1:2) H(:,4)];
end
XY1=ones(pointNum,3,'double');
XY1(:,1:2)=XY(:,1:2);
UV=zeros(pointNum,2,'double');
%% Algorithm.
UVu=H(1,:)*XY1'./H(3,:)*XY1';
UV(:,1)=UVu';
UVv=H(2,:)*XY1'./H(3,:)*XY1';
UV(:,2)=UVv';