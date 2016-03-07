%% Image calibrate based on view cone
% * Author: WANG Lei,USTB
%
% * Link: <https://github.com/shidafu/ViewConeCalibration.git>
%
% * Date:2016/3/3
%
% * Algorithom:  
%
% <<..\png\ViewConeModel.png>>
% 
% # *Get 2D Homography
% # Get C' ¦£c' ¦£u' ¦£v'
% # Get Zc from ellipse ¦£c'
% # Get Xc Yc from ellipse ¦£c'
% # *Or Get Zc from ellipse ¦£c'and focus
% 
% * Inputs:
%
%     Homography2D ---- UV must be center based but not lefttop based
%     imageHeight,imageWidth
%     %radius ---- radius of ¦£c < min(imageWidth,imageHeight)
%     'F',F ---- Object focus length if known
%     'F_Pixlength',F_Pixlength ---- Focus/Pixlength if known from other image
% 
% * Outputs:
%
%     A,Rd---- Inner and outer camera para matrix
%     Xc,Yc,Zc,C,O,thetaZc,phiZc---- axis vector and angle
%     F ---- Object focus length
%       Pixlength ---- center pixel real size in vertical object surface
%     Circle ---- Circle matrix in image plane
%                 [UV]'[Circle][UV]=0
%       ConePointsinUV ---- [CenterU,ConePointsinU1,ConePointsinU2,...;
%                            CenterV,ConePointsinV1,ConePointsinV2,...]
%       r ---- Radius of circle in UV plane
%     Ellipse ---- Ellipse matrix in XY plane
%                 [XY]'[Ellipse][XY]=0
%       ConePointsinXY ---- [CenterX,ConePointsinX1,ConePointsinX2,...;
%                            CenterY,ConePointsinY1,ConePointsiY2,...]
%       HEllipse,Rotate,Shift,a,b,c,e,phi,Focus,Peaks1,Peaks2
%                      ---- Ellipse paras in XY plane
%
function [A,Rd,...
          Xc,Yc,Zc,C,O,thetaZc,phiZc,...
          F,Pixlength,...
          Circle,ConePointsinUV,r,...
          Ellipse,ConePointsinXY,HEllipse,Rotate,Shift,a,b,c,e,phi,Focus,Peaks1,Peaks2]...
       = ViewConeCalibrate(varargin)
% function [A,Rd,...
%           Xc,Yc,Zc,C,O,thetaZc,phiZc,...
%           F,Pixlength,...
%           Circle,ConePointsinUV,...
%           Ellipse,ConePointsinXY,HEllipse,Rotate,Shift,a,b,c,e,phi,Focus,Peaks1,Peaks2]...
%     = ViewConeCalibriteBy2DCorners(Homography2D,imageHeight,imageWidth,radius,'F',F,'F_Pixarea',F_Pixarea)
% Get inputs
F=NaN;
F_Pixlength=NaN;
Homography2D=varargin{1};
imageHeight=varargin{2};
imageWidth=varargin{3};
% coneRadius=varargin{4};
% r=coneRadius;
if nargin>=5
    paraIndex=4;
    while paraIndex<nargin
        paraIndex=paraIndex+1;
        if(ischar(varargin{paraIndex}))
            if(strcmp(varargin{paraIndex},'F'))
                  if(nargin>=paraIndex+1)
                      if(~ischar(varargin{paraIndex+1})) 
                          F=varargin{paraIndex+1};
                      else
                          error('Input parameter area error');
                      end
                      paraIndex=paraIndex+1;
                      continue;
                  end
             elseif(strcmp(varargin{paraIndex},'F_Pixarea'))
                  if(nargin>=paraIndex+1)
                      if(~ischar(varargin{paraIndex+1})) 
                          F_Pixlength=varargin{paraIndex+1};
                      else
                          error('Input parameter area error');
                      end
                      paraIndex=paraIndex+1;
                      continue;
                  end
             end
        end
    end
end
% 1. *Get 2D Homography
Homography2D;

% 2. Get C' ConePointsinXY  ¦£c' ¦£u' ¦£v' 
% C:
CenterinUV=[0;0];
% C':
CenterinXY = GetXYFromUV(Homography2D,CenterinUV);
% ConePointsinUV
angleDiv=4;
angleNum=angleDiv*4;
angleAdd=2*pi/angleNum;
coneRadius=3*min(imageHeight,imageWidth)/8;%,min(height,width)/4
r=coneRadius;
conePointNum=angleNum+1;
ConePointsinUV=nan(2,conePointNum,'double');
for i=1:conePointNum
    ConePointsinUV(:,i)=double(CenterinUV);
end
for i=1:angleNum
    index=i+1;
    ConePointsinUV(1,index)=ConePointsinUV(1,index)+coneRadius*cos(angleAdd*(i-1));
    ConePointsinUV(2,index)=ConePointsinUV(2,index)+coneRadius*sin(angleAdd*(i-1));
end
% ConePointsinXY
% ConePointsinXY=nan(2,conePointNum,'double');
ConePointsinXY = GetXYFromUV(Homography2D,ConePointsinUV);

% ¦£c:[UV1]'[Circle][UV1]=0
Circle=[1,0,0;...
        0,1,0;...
        0,0,-coneRadius*coneRadius];
% ¦£c':[XY1]'[Ellipse][XY1]=0
Ellipse=Homography2D'*Circle*Homography2D;
[HEllipse,Rotate,Shift,a,b,c,e,phi,CenterinXY2,Focus,Peaks1,Peaks2]...
     = GetHorizontalEllipse(Ellipse);
[Focus,Peaks1] = FixDirection(ConePointsinXY,Focus,Peaks1);
CenterinXY-CenterinXY2;
% ¦£u:[LineU][UV1]=0
LineU=[1 0 0];
% ¦£u':[LineUinXY][XY1]=0
LineUinXY=LineU*Homography2D;
% ¦£v:LineV][UV1]=0
LineV=[0 1 0];
% ¦£v':[LineVinXY][XY1]=0
LineVinXY=LineV*Homography2D;

% 3. Get Zc from ellipse ¦£c'
%Vector2D in XY of camera/viewcone axis pointout from camera
%Angle of Zc axis with XY plane
ConePointsinXY1Fix=inv([Peaks1(1,2)-Peaks1(1,1) Peaks1(2,2)-Peaks1(2,1);...
               Peaks1(2,1)-Peaks1(2,2) Peaks1(1,2)-Peaks1(1,1)])*...
              ([ConePointsinXY(1,1) ConePointsinXY(2,1);Peaks1(2,1) -Peaks1(1,1)]*(Peaks1(:,2)-Peaks1(:,1)));
%plot(ConePointsinXY1Fix(1,:),ConePointsinXY1Fix(2,:), 'b+');
Zcxy=(Peaks1(:,1)-Peaks1(:,2))/norm(Peaks1(:,1)-Peaks1(:,2));
phiZc=acossin([Zcxy(1,1),Zcxy(2,1)]);
%Angle of Zc axis with XY plane !!!!!
D1D3=norm(Peaks2(:,1)-Peaks2(:,2));
D2C=norm(Peaks1(:,1)-ConePointsinXY1Fix);
D4C=norm(Peaks1(:,2)-ConePointsinXY1Fix);
if isnan(F)
    thetaZc=asin(D1D3*((1/D2C)+(1/D4C))/4);
    F=D1D3*(D2C+D4C)/(2*tan(thetaZc)*abs(D4C-D2C));
else
    thetaZc=acos(0.5*F*abs((1/D4C)-(1/D2C)));
end

% 4. Get Xc Yc from ellipse ¦£c'
%Vector3D in XYZ space
Zcxy=Zcxy*cos(thetaZc);
Zc=[Zcxy;-sin(thetaZc)];
normZc=norm(Zc);
Zc=Zc/normZc;
Xcxy=NaN;
dist3=norm(ConePointsinXY(:,2)-ConePointsinXY(:,1));
dist31=norm(ConePointsinXY(:,2+angleNum/2)-ConePointsinXY(:,1));
Xcxy=(ConePointsinXY(:,2)-ConePointsinXY(:,1))/dist3...
    +(ConePointsinXY(:,1)-ConePointsinXY(:,2+angleNum/2))/dist31;
Xcxy=Xcxy/2;
Xcxy=[Xcxy;0];
Xcv=cross(Zc,Xcxy);% vertical inter vector3D
Xc= cross(Xcv,Zc);
Xc=Xc/(norm(Xc));
Ycxy=NaN;
dist4=norm(ConePointsinXY(:,2+angleNum/4)-ConePointsinXY(:,1));
dist41=norm(ConePointsinXY(:,2+3*angleNum/4)-ConePointsinXY(:,1));
Ycxy=(ConePointsinXY(:,2+angleNum/4)-ConePointsinXY(:,1))/dist4...
    +(ConePointsinXY(:,1)-ConePointsinXY(:,2+3*angleNum/4))/dist41;
Ycxy=Ycxy/2;
Ycxy=[Ycxy;0];
Ycv=cross(Zc,Ycxy);% vertical inter vector3D
Yc= cross(Ycv,Zc);
Yc=Yc/(norm(Yc));

% 5. *Or Get Zc from ellipse ¦£c'and focus
%  Take the chord vertical with major axis of ellipse and through
%         circle center's mapping point inside ellipse,the chord's length
%         is in direct proportion to object Distance;
chordVector=Peaks2(:,1)-Peaks2(:,2);
syms x y real;
[solutions_x, solutions_y]=solve([x y 1]*M*[x;y;1]==0,...
    (x-ConePointsinXY(1,1))*chordVector(2,1)==(y-ConePointsinXY(2,1))*chordVector(1,1),x,y);
chordVector=[double(solutions_x(1))-double(solutions_x(2)),...
    double(solutions_y(1))-double(solutions_y(2))];
chordLength=sqrt(chordVector*chordVector');
Pixlength=chordLength/(2*coneRadius);
defaultObjectDistance=100;
if ~isnan(F)
else
    if ~isnan(F_Pixlength)
        F=Pixlength*F_Pixlength;
    else
        F=defaultObjectDistance;
    end
end
xc=ConePointsinXY(1,1);
yc=ConePointsinXY(2,1);
zc=0;
O=[xc;yc;zc];
xc=xc-F*Zc(1,1);
yc=yc-F*Zc(2,1);
zc=zc-F*Zc(3,1);
C=[xc;yc;zc];
end

%% Fix the direction of view cone ellipse
% * Author: WANG Lei,USTB
%
% * Link: <https://github.com/shidafu/ViewConeCalibration.git>
%
% * Date:2016/3/3
%
% * Inputs:
%
%     ConePointsinUV ---- UV must be center based but not lefttop based
%     Focus  ---- 2 by 1 array,Focus point of the source ellipse;
%     Peaks1 ---- 2 by 1 array,Peak points in long axis of the source ellipse;
% 
% * Outputs:
%
%     Focus  ---- 2 by 1 array,Focus point of the source ellipse;
%     Peaks1 ---- 2 by 1 array,Peak points in long axis of the source ellipse;
%
function [Focus,Peaks1] = FixDirection(ConePointsinXY,Focus,Peaks1)
[hXY,wXY]=size(ConePointsinXY);
if hXY==2
    pointNum=wXY;
elseif wXY==2
    ConePointsinXY=ConePointsinXY';
    pointNum=hXY;
else
    error('Error size of input matrix!');
end
%sum point left or right
leftNum=0;
rightNum=0;
for i=1:pointNum
    if norm(Peaks1(:,1)-ConePointsinXY(:,i))<norm(Peaks1(:,2)-ConePointsinXY(:,i));
        leftNum=leftNum+1;
    else
        rightNum=rightNum+1;
    end
end
if leftNum>rightNum
    tmp=Peaks1;
    Peaks1(:,1)=tmp(:,2);
    Peaks1(:,2)=tmp(:,1);
end
if norm(Peaks1(:,1)-Focus(:,1))>norm(Peaks1(:,2)-Focus(:,1));
    tmp=Focus;
    Focus(:,1)=tmp(:,2);
    Focus(:,2)=tmp(:,1);
end
end