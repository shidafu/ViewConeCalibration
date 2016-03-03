%% Get UV from XY according to 2D Homography of a projective geometry transform
% * Author: WANG Lei,USTB
%
% * Link: <https://github.com/shidafu/ViewConeCalibration.git>
%
% * Date:2016/3/3
%
% * Algorithom:
%
% Get [ _UV_ ] By solving:
%
% $$\mathit{s} \left[\begin{array}{c} \mathit{u}\\ \mathit{v}\\ 1 \end{array}\right]=
% \left[\begin{array}{ccc} 
% \mathit{h}_{11} & \mathit{h}_{12} & \mathit{h}_{13}\\
% \mathit{h}_{21} & \mathit{h}_{22} & \mathit{h}_{23}\\
% \mathit{h}_{31} & \mathit{h}_{32} & 1 \end{array}\right]
% \left[\begin{array}{c} \mathit{x}\\ \mathit{y}\\ 1 \end{array}\right]$$
%
% * Inputs:
%
%     H----3 by 3 matrix: [h11 h12 h13;
%                          h21 h22 h23;
%                          h31 h32  1]
%     XY----cordNum by pointNum matrix,
%                  cordNum==2,pointNum>=4,
%                  [x1,x2,...;
%                   y1,y2,...]
%
% * Outputs:
%
%     UV----cordNum by pointNum matrix,
%                  cordNum==2,pointNum>=4,
%                  [u1,u2,...;
%                   v1,v2,...]
function UV = GetUVFromXY(H,XY)
% Initial
[cordNum, pointNum]=size(XY);
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
XY1=ones(3,pointNum,'double');
XY1(1:2,:)=XY(1:2,:);
UV=zeros(2,pointNum,'double');
% Algorithm
UVu=H(1,:)*XY1./H(3,:)*XY1;
UV(:,1)=UVu;
UVv=H(2,:)*XY1./H(3,:)*XY1;
UV(:,2)=UVv;