clc;
clear;
load('Calib3.mat')

Gaze(:,1) = Data.gazeL(:,7);
Gaze(:,3) = Data.gazeL(:,8);

Gaze(:,2) = Data.gazeR(:,7);
Gaze(:,4) = Data.gazeR(:,8);


Gaze(Gaze(:,1) < 0, :)=[];
Gaze(Gaze(:,2) < 0, :)=[];
Gaze(Gaze(:,3) < 0, :)=[];
Gaze(Gaze(:,4) < 0, :)=[];


Gaze(:,1:2) = round(Gaze(:,1:2)*1920);
Gaze(:,3:4) = round(Gaze(:,3:4)*1080);

Gaze(Gaze(:,1) > 1920, :)=[];
Gaze(Gaze(:,2) > 1920, :)=[];
Gaze(Gaze(:,3) > 1080, :)=[];
Gaze(Gaze(:,4) > 1080, :)=[];

VisMatrix = zeros(1920,1080);
Matrix=zeros(1920,1080);
for i = 1:length(Gaze)
    xl = Gaze(i,1);
    yl = Gaze(i,3);
    xr = Gaze(i,2);
    yr = Gaze(i,4);
    Matrix(xl,yl) = Matrix(xl,yl) + 1;
    Matrix(xr,yr) = Matrix(xr,yr) + 1;
    if ((xl + 3) && (xr + 3) <= 1920) && ((yl + 3)&&(yr + 3)<=1080) && ((xl - 3) && (xr - 3) > 0) && ((yl - 3)&&(yr - 3)>0)
        VisMatrix((xl-3):(xl+3),(yl-3):(yl+3)) = VisMatrix((xl-3):(xl+3),(yl-3):(yl+3)) + 1;
        VisMatrix((xr-3):(xr+3),(yr-3):(yr+3)) = VisMatrix((xr-3):(xr+3),(yr-3):(yr+3)) + 1;
    else
    VisMatrix(xl,yl) = VisMatrix(xl,yl) + 1;
    VisMatrix(xr,yr) = VisMatrix(xr,yr) + 1;
    end
    end
Vind = find(Matrix >= 1);

figure(3)
axis ([0 1920 0 1080])
imagesc(VisMatrix)

axis ij
hold on




x=80;
C = Ellipse(x,x*pi);
Csize = size(C);
CirMatrix = zeros(1920,1080);
% find coordinates on Matrix to place circle at
x1 = 1920/2-floor(Csize(1)/2)+1;
x2 = 1920/2+ceil(Csize(1)/2);
y1 = 1080/2-floor(Csize(2)/2)+1;
y2 = 1080/2+ceil(Csize(2)/2);
% Place circle into Matrix
VisMatrix(x1:x2,y1:y2) = C;
CirMatrix(x1:x2,y1:y2) = C;
Cind = find(CirMatrix == 1);

% ValidCount = size(intersect(Vind,Cind));
% if ValidCount/2 > length(Gaze)/2
% else
% end

imagesc(VisMatrix)





% intersect(Vind, Cind)

% 
% 
% figure(1)
% scatter(Gaze(:,1)*1920,Gaze(:,2)*1080,'r')
% 
% hold on
% scatter(Gaze(:,3)*1920,Gaze(:,4)*1080,'b')
% 
% C = Circle(100)
% imagesc(C==1)
% %viscircles([1920/2,1080/2],100,'Color','k')
% 
% xlim([0 1920])
% ylim([0 1080])
% hold off
% 
% figure(2)
% C = viscircles([1920/2,1080/2],100,'Color','k')
% 
% 
% x = 0:1920;
% y = 0:1080;
% [X Y] = meshgrid(x,y);
% u = zeros(size(X));
% u((X.^2+Y.^2)<100^2)=1;   % radius 100, center at the origin
% % hard boundary
% figure(1)
% imagesc(u)
% % weight the points: point itself; average of nearest neighbors;
% % averaged of diagonal neighbors.  These must add up to 1.
% wp = .4;  wn = .4;  wd = .2;
% ind = 2:length(x)-1;
% u(ind,ind) = wp*u(ind,ind) ...
%   + (wn/4)*(u(ind-1,ind  ) + u(ind+1,ind  ) + u(ind  ,ind-1) + u(ind  ,ind+1) ) ...
%   + (wd/4)*(u(ind-1,ind-1) + u(ind-1,ind+1) + u(ind+1,ind-1) + u(ind+1,ind+1) );
% % extended boundary
% fig(2)
% imagesc(u)