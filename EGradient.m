clear 
clc

R = 1.65;
v = 1e-6;
g = 9.81;

sInput = 1.6e-4;
hInput = 5;
bInput = 200/2;
uInput = 1.46;
Cf = 5e-3;

rng('default')
year = 365*24*60*60;
noNode = 270;%number of modeling node
noYear = 1000;
yearSkip = 50;
tf = 1;
noTime = noYear/tf;%modeling time steps
T1 = nan(1,noTime);%first cut off time
xStep = bInput*3;
xlimit = noNode*xStep;
x = xStep:xStep:xlimit;
y = zeros(1,noNode) + rand(1,noNode)*bInput/10;
xNor = x/x(end)*25*pi;
xi = x;

X = nan(noYear,1000);
Y = nan(noYear,1000);
SINU = nan(noYear,1000);
XOver = nan(noYear,1000);
YOver = nan(noYear,1000);
AccuSOver = nan(noYear,1000);
Etotal = nan(noYear,1000);

us0 = uInput; %normal flow velocity
h0 = hInput; %normal flow depth
b = bInput;%half channel width
g = 9.81;%gravitational acceleration
A = 10; %need to verify
alfa = 0.077;
Chi1 = alfa/sqrt(Cf);
Chi = Chi1 - 1/3;
As = 181 * (h0/b)^2 * 1/Chi1 * (2*Chi^2 + 0.8*Chi + 1/15);
A = A + As -1;
k = 80; %searching index range for cutoff point
E = fliplr(linspace(2,8,noNode)*1e-7);
Ei = E;

%allocate for parameters
curvature = nan(1,noNode);
dCurve = nan(1,noNode);
s = nan(1,noNode);
accuS = nan(1,noNode);
ds = nan(1,noNode);
dx1 = nan(1,noNode);
dx2 = nan(1,noNode);
dy1 = nan(1,noNode);
dy2 = nan(1,noNode);
usb = zeros(1,noNode);
migrRate = zeros(1,noNode);
maxSinu = zeros(1,noYear);
aveSinu = zeros(1,noYear);
maxSinu2 = zeros(1,noYear);
aveSinu2 = zeros(1,noYear);
aveMigRate = zeros(1,noYear);
maxMigRate = zeros(1,noYear);
aveCurve = zeros(1,noYear);
maxCurve = zeros(1,noYear);
sinu25 = zeros(1,noYear);
sinu75 = zeros(1,noYear);
coXY = [];

figure(1);hold on
% subplot(2,1,1);hold on
% h = plot(x/1e3,y/1e3,'-k','linewidth',1.5);
% ylim([-20000 20000]/1e3);
% xlim([0 110000]/1e3);
% xlim([0 3])
set(gca,'DataAspectRatio', [1 1 1]);
set(gcf,'Unit','in','position',[.5,.5,12,6]);

for t=1:noTime
    %calculate alongstream distance, curvature and velocity increment
    ds(2:end) = sqrt((x(2:end)-x(1:end-1)).^2 + (y(2:end)-y(1:end-1)).^2);
    ds(1) = ds(2);
    dx1(2:end) =  (x(2:end)-x(1:end-1))./ds(2:end);
    dx1(1) = dx1(2);
    dx2(2:end) = (dx1(2:end)-dx1(1:end-1))./ds(2:end);
    dx2(1) = dx2(2);
    dy1(2:end) = (y(2:end)-y(1:end-1))./ds(2:end);
    dy1(1) = dy1(2);
    dy2(2:end) = (dy1(2:end)-dy1(1:end-1))./ds(2:end);
    dy2(1) = dy2(2);
    curvature = (dx1.*dy2-dy1.*dx2)./(dx1.^2+dy1.^2).^(3/2);
    dCurve(2:end) = (curvature(2:end)-curvature(1:end-1))./ds(2:end);
    dCurve(1) = dCurve(2);
    dnx = (-dy1)./sqrt(dx1.^2+dy1.^2);
    dny = dx1./sqrt(dx1.^2+dy1.^2);
    usb(2:end) = b./(us0./ds(2:end) + 2*us0./h0*Cf) .*...
        (-us0.^2.*dCurve(2:end) + ...
        Cf*curvature(2:end) .* ((us0.^4/g./h0.^2) + A*us0.^2./h0)...
        + us0./ds(2:end) .* usb(1:end-1)/b);
    usb(1) = 0;
    migrRate = E.*usb*year*tf;
    L = sqrt(dnx.^2+dny.^2);
    dx = -dnx.*migrRate;
    dy = -dny.*migrRate;
    dx(end) = 0;
    dy(end) = 0;
    x = x + dx;
    y = y + dy;

    %save x, curvature and s from previous time step,
    x_t = x;
    y_t = y;
    curvature_t = curvature;

    accuS = streamLineDistance(ds);
    accuS_t = accuS;%
    maxCurve(t) = max(abs(curvature_t));   
    aveCurve(t) = mean(abs(curvature_t));
    aveMigRate(t) = mean(abs(migrRate));
    maxMigRate(t) = max(abs(migrRate));

    %calculate cutoff
    [x,y,T1,noNode,coxy]=cutOff(noNode,k,x,y,b,t,T1);
    coXY = [coXY,coxy];
    %calculate cutoff for the most downstream section
    [x,y,T1,noNode]=cutOff2(noNode,k,x,y,b,t,T1);
    
    %new channel path
    [~,noNode] = size(x);
    ds = nan(1,noNode);
    ds(2:end) = sqrt((x(2:end)-x(1:end-1)).^2 + (y(2:end)-y(1:end-1)).^2);
    ds(1) = 0;
    
    accuS = streamLineDistance(ds);
    s_reMesh = accuS(1):xStep:accuS(end);
    x_reMesh = interp1(accuS,x,s_reMesh,'spline');%'makima'
    y_reMesh = interp1(accuS,y,s_reMesh,'spline');

    x_reMesh = [x_reMesh, xlimit];%add the end point, so model is stable
    y_reMesh = [y_reMesh, 0];
    x = x_reMesh;
    y = y_reMesh;

    y(1) = y(1) + (rand(1,1)-0.5)*1;% small perturbation at the upstream, without this, upstream become straight    
    y(end) = y(end-1) + (rand(1,1)-0.5)*1;

    [~,noNode] = size(x);
    E = nan(1,noNode);%resize E, bucause there are cutoff elements.
    [x_sort,idx_sort] = sort(x);
    E(idx_sort) = interp1(xi,Ei,x_sort);

    us0_t = us0;
    usb_t = usb;
    migrRate_t = migrRate;
    curvature = nan(1,noNode);
    dCurve = nan(1,noNode);
    s = nan(1,noNode);
    accuS = nan(1,noNode);
    ds = nan(1,noNode);
    dx1 = nan(1,noNode);
    dx2 = nan(1,noNode);
    dy1 = nan(1,noNode);
    dy2 = nan(1,noNode);
    usb = zeros(1,noNode);
    migrRate = zeros(1,noNode);

    %calculate intrinsic sinuosity
    output = calculateSinu2(x_t,y_t,accuS_t,b);
    sinu2 = output.sinu2;
    xOver = output.xOver;
    yOver = output.yOver;
    accuSOver = output.accuSOver;

    %store sinuosity through time
    [~,xSize] = size(x_t);
    [~,sinuSize] = size(sinu2);
    X(t,1:xSize) = x_t;
    Y(t,1:xSize) = y_t;
    SINU(t,1:sinuSize) = sinu2;
    XOver(t,1:sinuSize) = xOver(2:end);%there're n xy coordinates and n-1 sinu meansurements
    YOver(t,1:sinuSize) = yOver(2:end);
    AccuSOver(t,1:sinuSize) = accuSOver(2:end);
    
    % assign erosion coefficient
    E_t = nan(1,xSize);
    [x_sort,idx_sort] = sort(x_t);
    E_t(idx_sort) = interp1(xi,Ei,x_sort);
    Etotal(t,1:xSize) = E_t;
    
    %plotting
    movingWindows = linspace(0,noYear,noYear/yearSkip);
    cmap = parula(length(movingWindows));
    titleText = string(sprintf('t = %d year',t));

    if t/yearSkip-round(t/yearSkip)==0
        figure(1);hold on
        title(titleText);
        h = plot(x/1e3,y/1e3,'-k','linewidth',1.5);
        pause(0.00000001);
        delete(h)
        set(gca,'DataAspectRatio', [1 1 1]);            
    end

end

close all