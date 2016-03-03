%% Get 2D Homography of a projective geometry transform
% Author:WANG Lei,USTB,email£ºtowanglei@163.com
%
% Date:2016/3/3
%
% Find _s_ and  [ _H_ ]:
%
% $$\mathit{s} \left(\begin{array}{c} \mathit{u}\\ \mathit{v}\\ 1 \end{array}\right)=
% \left(\begin{array}{ccc} 
% \mathit{h}_{11} & \mathit{h}_{12} & \mathit{h}_{13}\\
% \mathit{h}_{21} & \mathit{h}_{22} & \mathit{h}_{23}\\
% \mathit{h}_{31} & \mathit{h}_{32} & 1 \end{array}\right)
% \left(\begin{array}{c} \mathit{x}\\ \mathit{y}\\ 1 \end{array}\right)$$
%
% Algorithom:
%
% $$ \left(\begin{array}{cccccccc}
% \mathit{x}_{1} & \mathit{y}_{1} &1 &0 &0 &0 & \mathit{-u}_{1}\mathit{x}_{1}& \mathit{-u}_{1}\mathit{y}_{1}\\
% 0 &0 &0 & \mathit{x}_{1} & \mathit{y}_{1} &1 & \mathit{-v}_{1}\mathit{x}_{1}& \mathit{-v}_{1}\mathit{y}_{1}\\
%                               &   &   &   & \vdots\\
% \mathit{x_{n}} & \mathit{y_{n}} &1 &0 &0 &0 & \mathit{-u_{n}}\mathit{x_{n}}& \mathit{-u_{n}}\mathit{y_{n}}\\
% 0 &0 &0 & \mathit{x_{n}} & \mathit{y_{n}} &1 & \mathit{-v_{n}}\mathit{x_{n}}& \mathit{-v_{n}}\mathit{y_{n}}\\
% \end{array}\right)
% \left(\begin{array}{c} 
% \mathit{h}_{11} \\ \mathit{h}_{12} \\ \mathit{h}_{13}\\
% \mathit{h}_{21} \\ \mathit{h}_{22} \\ \mathit{h}_{23}\\
% \mathit{h}_{31} \\ \mathit{h}_{32}
% \end{array}\right)=
% \left(\begin{array}{c} 
% \mathit{u}_{1} \\ \mathit{v}_{1} \\ \vdots \\ \mathit{u_{n}} \\ \mathit{v_{n}}
% \end{array}\right)$$
%
% $$[\mathit{H_{list}}]=\left([\mathit{XYUV_{list}}]^{T}\cdot [\mathit{XYUV_{list}}]\right)^{-1}
%                        \cdot[\mathit{XYUV_{list}}]^{T}\cdot [\mathit{UV_{list}}]$$
%
% $$\mathit{s}=[\mathit{H}]\cdot[\mathit{XY1}]^{T}\cdot[\mathit{UV1}]
%             \cdot\left( [\mathit{UV1}^{T} \cdot [\mathit{UV1} \right)$$
% 

function [H,s] = GetHomography2D(XY,UV)
% Get 2D Homography of a projective geometry transform
%
% Inputs:
%     UV----pointNum by cordNum matrix,
%                  cordNum==2,pointNum>=4,
%                  [u1,v1;
%                   u2,v2;
%                     :]
%     XY----pointNum by cordNum matrix,
%                  cordNum==2,pointNum>=4,
%                  [x1,y1;
%                   x2,y2;
%                     :]
%
% Outputs:
%     H----3 by 3 matrix: [h11 h12 h13;
%                          h21 h22 h23;
%                          h31 h32  1]
%     s----projective para
%% Initial
[pointNum, cordNum]=size(XY);
if ~(cordNum==2 || cordNum==3)
    error('Input matrix size error!');
end
XY1=ones(pointNum,3,'double');
XY1(:,1:2)=XY(:,1:2);
UV1=ones(pointNum,3,'double');
UV1(:,1:2)=XY(:,1:2);
HList=zeros(8,1,'double');
XYUVList=zeros(2*pointNum,8,'double');
UVList=zeros(2*pointNum,1,'double');
for i=1:pointNum
    XYUVList(i*2-1,:)=[XY(i,1), XY(i,2), 1, 0, 0, 0, -1*UV(i,1)*XY(i,1), -1*UV(i,1)*XY(i,2)];    
    XYUVList(i*2,:)=[0, 0, 0, XY(i,1), XY(i,2), 1, -1*UV(i,2)*XY(i,1), -1*UV(i,2)*XY(i,2)];
    UVList(i*2-1,:)=UV(1,i);
    UVList(i*2,:)=UV(2,i);
end
%% Algorithm.
HList=inv(XYUVList'*XYUVList)*XYUVList'*UVList;
%% Set outputs.
H(1,:)=HList(1:3,1);
H(2,:)=HList(4:6,1);
H(3,1:2)=HList(7:8,1);
H(3,3)=1;
s=H*XY1'*UV1*inv(UV1'*UV1);