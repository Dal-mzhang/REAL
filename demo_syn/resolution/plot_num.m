%plot ps distribution in 3d
res = 'resolution.txt'
a = load(res);
lat = a(:,2);
lon = a(:,3);
dep = a(:,4);
ps = a(:,7);

nx = floor((max(lon) - min(lon))/0.01);
ny = floor((max(lat) - min(lat))/0.01);
nz = floor((max(dep) - min(dep))/1.0);

nx = 61;ny = 61; nz = 41;

xlin = linspace(min(lon),max(lon),nx);
ylin = linspace(min(lat),max(lat),ny);
zlin = linspace(min(dep),max(dep),nz);
[X,Y,Z] = meshgrid(xlin,ylin,zlin);
R = griddata(lon,lat,dep,ps,X,Y,Z);

xslice =  13.25; yslice = 42.75;  zslice = 10;

slice(X,Y,Z,R,xslice,yslice,zslice)
xlim([13.15,13.35]);
ylim([42.6,42.9]);
%xlabel('East-West');ylabel('North-South'); zlabel('Depth (km)');
set(gca,'FontSize',15);
xlabel('East-West','FontSize',15);ylabel('North-South','Fontsize',15); zlabel('Depth (km)','Fontsize',15);
colormap jet; 
set(gca,'zDir','reverse');
colorbar('FontSize',12,'Position',[0.92 0.27 0.025 0.44])
saveas(gcf,'num.pdf')

