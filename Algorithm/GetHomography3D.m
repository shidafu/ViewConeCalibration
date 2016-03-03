%% Get 3D Homography of a projective geometry transform
% Author:WANG Lei,USTB,email£ºtowanglei@163.com
%
% Date:2016/3/3
%
% Find _s_ and  [ _H_ ]:
%
% $$\mathit{s} \left(\begin{array}{c} \mathit{u}\\ \mathit{v}\\ 1 \end{array}\right)=
% \left(\begin{array}{cccc} 
% \mathit{h}_{11} & \mathit{h}_{12} & \mathit{h}_{13} & \mathit{h}_{14}\\
% \mathit{h}_{21} & \mathit{h}_{22} & \mathit{h}_{23} & \mathit{h}_{24}\\
% \mathit{h}_{31} & \mathit{h}_{32}  & \mathit{h}_{33} & 1 \end{array}\right)
% \left(\begin{array}{c} \mathit{x}\\ \mathit{y}\\ \mathit{z}\\ 1 \end{array}\right)$$
%
% Algorithom:
%
% $$ \left(\begin{array}{ccccccccccc}
% \mathit{x}_{1} & \mathit{y}_{1} & \mathit{z}_{1} &1 &0 &0 &0 &0 & \mathit{-u}_{1}\mathit{x}_{1}& \mathit{-u}_{1}\mathit{y}_{1}& \mathit{-u}_{1}\mathit{z}_{1}\\
% 0 &0 &0 &0 & \mathit{x}_{1} & \mathit{y}_{1} & \mathit{z}_{1} &1 & \mathit{-v}_{1}\mathit{x}_{1}& \mathit{-v}_{1}\mathit{y}_{1}& \mathit{-v}_{1}\mathit{z}_{1}\\
%                               &   &   &   &   &   & \vdots\\
% \mathit{x_{n}} & \mathit{y_{n}} & \mathit{z_{n}} &1 &0 &0 &0 &0 & \mathit{-u_{n}}\mathit{x_{n}}& \mathit{-u_{n}}\mathit{y_{n}}& \mathit{-u_{n}}\mathit{z_{n}}\\
% 0 &0 &0 &0 & \mathit{x_{n}} & \mathit{y_{n}} & \mathit{z_{n}} &1 & \mathit{-v_{n}}\mathit{x_{n}}& \mathit{-v_{n}}\mathit{y_{n}}& \mathit{-v_{n}}\mathit{z_{n}}\\
% \end{array}\right)
% \left(\begin{array}{ccccc} 
% \mathit{h}_{11} \\ \mathit{h}_{12} \\ \mathit{h}_{13}\\ \mathit{h}_{14}\\
% \mathit{h}_{21} \\ \mathit{h}_{22} \\ \vdots\\
% \mathit{h}_{31} \\ \mathit{h}_{32}
% \end{array}\right)=
% \left(\begin{array}{c} 
% \mathit{u}_{1} \\ \mathit{v}_{1} \\ \vdots \\ \mathit{u_{n}} \\ \mathit{v_{n}}
% \end{array}\right)$$
%
% $$[\mathit{H_{list}}]=\left([\mathit{XYZUV_{list}}]^{T}\cdot [\mathit{XYZUV_{list}}]\right)^{-1}
%                        \cdot[\mathit{XYZUV_{list}}]^{T}\cdot [\mathit{UV_{list}}]$$
%
% $$\mathit{s}=[\mathit{H}]\cdot[\mathit{XYZ1}]^{T}\cdot[\mathit{UV1}]
%             \cdot\left( [\mathit{UV1}^{T} \cdot [\mathit{UV1} \right)$$
% 

function [H,s] = GetHomography3D(XYZ,UV)
% Get 2D Homography of a projective geometry transform
%
% Inputs:
%     UV----pointNum by cordNum matrix,
%                  cordNum==2,pointNum>=4,
%                  [u1,v1;
%                   u2,v2;
%                     :]
%     XY----pointNum by cordNum matrix,
%                  cordNum==3,pointNum>=4,
%                  [x1,y1,z1;
%                   x2,y2,z1;
%                     :]
%
% Outputs:
%     H----3 by 3 matrix: [h11 h12 h13 h14;
%                          h21 h22 h23 h24;
%                          h31 h32 h33 1]
%     s----projective para
%% Initial
[pointNum, cordNum]=size(XYZ);
if ~(cordNum==3 || cordNum==4)
    error('Input matrix size error!');
end
XYZ1=ones(pointNum,4,'double');
XYZ1(:,1:3)=XYZ(:,1:3);
UV1=ones(pointNum,3,'double');
UV1(:,1:2)=XYZ(:,1:2);
HList=zeros(11,1,'double');
XYZUVList=zeros(2*pointNum,11,'double');
UVList=zeros(2*pointNum,1,'double');
for i=1:pointNum
    XYZUVList(i*2-1,:)=[XYZ(i,1), XYZ(i,2), XYZ(i,3), 1, 0, 0, 0, 0, -1*UV(i,1)*XYZ(i,1), -1*UV(i,1)*XYZ(i,2), -1*UV(i,1)*XYZ(i,3)];    
    XYZUVList(i*2,:)=[0, 0, 0, 0, XYZ(i,1), XYZ(i,2), XYZ(i,3), 1, -1*UV(i,2)*XYZ(i,1), -1*UV(i,2)*XYZ(i,2), -1*UV(i,2)*XYZ(i,3)];
    UVList(i*2-1,:)=UV(1,i);
    UVList(i*2,:)=UV(2,i);
end
%% Algorithm.
HList=inv(XYZUVList'*XYZUVList)*XYZUVList'*UVList;
%% Set outputs.
H(1,:)=HList(1:4,1);
H(2,:)=HList(5:8,1);
H(3,1:3)=HList(9:11,1);
H(3,4)=1;
s=H*XYZ1'*UV1*inv(UV1'*UV1);