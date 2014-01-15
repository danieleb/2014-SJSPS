function drawplane(Y,c)
if ~nargin, unittest; return; end
if nargin<2, c='k'; end

a = Y(:,1:3); % select 3 points in Y
% find normal vector
a2ma1 = a(:,2) - a(:,1);
a3ma1 = a(:,3) - a(:,1);
n = cross(a2ma1,a3ma1);
% plot the plane
d = n'*a(:,1);
alpha = 1.5;
xmin = alpha*min(Y(1,:));
xmax = alpha*max(Y(1,:));
ymin = alpha*min(Y(2,:));
ymax = alpha*max(Y(2,:));
x = linspace(xmin,xmax,10);
y = linspace(ymin,ymax,10);
[X,Y] = meshgrid(x,y);
Z = (d - (n(1)*X) - (n(2)*Y))/n(3);
h = surf(X,Y,Z);
set(h,'FaceColor',c,'FaceAlpha',0.3,'LineStyle','none');

function unittest
close all
X = randn(3,100);
D = pca(X');
Y = D(:,1:2)*pinv(D(:,1:2))*X;
scatter3(Y(1,:),Y(2,:),Y(3,:));
hold on
drawplane(Y);