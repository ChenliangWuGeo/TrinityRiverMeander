% Moving Boundary method from Parker et al., 2008b
clear
clc

tic
seaLevel=0;

%Grid setup
%Distance=500000; % Total modeling distance
noNodes = 100; % number of modeling nodes
dx = 1000; % x increment, unit is meter
dxp = dx; % previous dx for sed partition calculation
dxn = 1/noNodes; % distretized interval in moving boundary coordinates.
xn=0:dxn:1; %dimensionless grid

%simulation time-stepping
tf = 1; % ratio to adjust number of time steps and sea level rate
noTimeSteps = 1000/tf; % modeling time period
ti = 500/tf;  % sample interval for ploting figure
dt = 60*60*24*365*tf; % time increment of 1 year X tf

%sampling simulation results
noSkipSteps = 10/tf; %grab system information every noSkipSteps*dt time interval

%system thresholds (for autoretreat calculation)
rS_threshold = 0.3;
rS = ones(1,noNodes+1); %initial ratio of Sfrix/Sx

%channel system slope
Sfi = 1.6e-4; % initial slope of fluvial reach%system slopes
Sbb = 4.6e-4; % Slope of subaerial bedrock reach
Ssb = 1e-5; % slope of subaqueous basement
Sfore = 1.6e-3; % slope of foreset.
Sx = zeros(1,noNodes+1); % channel slope over x space

% fitted initial channel bed elevation and x node location
x = 0:dx:dx*noNodes;
x = x./1000;
xp = x; % x coordinates for previous time step
be = zeros(1,noNodes+1); %preallocate bed elevation
% be = 2.241e-5 * x.^2 + 0.03048 * x -23.86; % bed elevation for alluvial domain.
% % be = 0.0001151 * x.^2 - 0.008928 * x -20.85; % fitted bed elevation for just backwater reach
% be = fliplr(be);
be(end) = -5; %bed elevation
be(1) = dx*noNodes*Sfi+be(end);%bed elevation
be = be(1):((be(end)-be(1))/noNodes):be(end);%bed elevation

bep = be; % bed elevation for previous time step
aggradationRate = be; % bed elevation increment/aggradation rate per year

%inital basin configuration
bebasei = -8; % initial bed elevation at foreset-basement break
sbai = -dx*noNodes; % initial x coordinate at bed rock/alluvial contact
stfi = 0; % initial x coordinate at shoreline break
ssbi = stfi+(be(end)-bebasei)/Sfore; % initial x coodinate of foreset bottom set boundary
sba = sbai; % x coordinate at bed rock/alluvial contact
stf = stfi; % x coordinate at shoreline break
ssb = ssbi; % x coordinate at foreset/bottom set contact.
bebase = be(end)-Sfore*(ssb-stf); % elevation of foreset/basement break, the value is updated each step, the value here is not important.
bebasep = bebase; %elevation of basement break for previous step, foreset increment calculation
autosb = noNodes+1; % index of shoreline break during shoreline break
stfauto = []; % location of shoreline location during autoretreat

%basin rate variables
rsba = 0; % rate of migration of bedrock/alluvial contact 
rstf = 0; % rate of migration of shoreline break
rssb = 0; % rate of migration of foreset/bottom set contact

%storage variables
STF = zeros(1,noTimeSteps); % x coordinate at shoreline break throughout the modeling time frame
RSFT = zeros(1,noTimeSteps);% rate of migration of shoreline break for all the timed interval
STFauto = zeros(1,noTimeSteps); % collection of shoreline location during autoretreat
migRate = zeros(1,noTimeSteps);

%preallocation of derivatives
dndt = zeros(1,noNodes+1); %elevation time der.
dqdxn = zeros(1,noNodes+1);  %sed flux time der.
dFdxn = zeros(11,noNodes+1); %grainsize classes time derv.
dLadxn = zeros(1,noNodes+1); %Active thickness of active layer time der
dDgdx = zeros(1,noNodes+1); %spatial change in median grainsize

x = xn.*(stf-sba)+sba;
figure(1)
hold on
plot(x./1e3,be,'k')

Sx = zeros(1,noNodes+1); % slope over x space
Sfrix = zeros(1,noNodes+1); % friction slope over x space
dSdx = zeros(1,noNodes+1); % slope change over x space
dS = zeros(1,noNodes+1); % slope change over x space

% hydraulics
un = 1.46;% uniform flow velocity m/s during flood
seaLevelRate = 1e-3*tf; % rate of sea level rise m/yr
If = 0.1;
Cf = 0.0036; % friction coefficient
kc = 75e-3; % Roughness height
ar = 8.1; % Coefficient in Manning-Strickler resistance relation
g = 9.81; % gravitational acceleration
qw = (kc/ar^6*un^10/g^3/Sfi^3)^(1/4); % calculated qw
% qw = 700/200; % discharge/width
bei = be; % set the initial bed elevation for comparison
sp = zeros(1,noNodes); % Flow surface profile
Hn = qw/un; % normal flow depth
qs = zeros(1,noNodes); % sediment flux
D50 = 250e-6; % medium grain size
R = 1.65; % Submerged specific gravity
au = 0.7; % explicit vs. implicit
LAMBDA = 1; % fraction of mud deposited in flood plain
Sinu = 1.7; % Sinuosity
rB = 36; % ratio between channel and flood plain
lambda = 0.4; % bed porosity

% flow depth calculation preallocation
H = zeros(1,noNodes+1); % flow depth.
H(end) = seaLevel-be(end); % flow depth this sets downstream boundary value as well
Hn = qw/un; % normal flow depth
H(1:noNodes) = Hn;

% % grain size group and initial fractions for each group
% GSGr = [-5 -4.5 -4 -3.5 -3 -2.5 -2 -1.5 -1 -0.5 0]; % grain size group
% noGrainSizeGroups = length(GSGr);
% 
% %grain size calculation preallocation 
% D = zeros(noGrainSizeGroups,noNodes+1); % grain size groups
% DgPhi = zeros(1,noNodes+1); % geometric mean in phi scale
% DgFsPhi = zeros(1,noNodes+1); % geometric mean in phi scale in substrate
% Dv = zeros(1,noNodes+1); % Deviation of grain size
% DvFs = zeros(1,noNodes+1); % Deviation of grain size in substrate
% SDv = zeros(1,noNodes+1); % standard deviation
% SDvFs = zeros(1,noNodes+1); % standard deviation
% La = zeros(1,noNodes+1); % thickness of active layer
% Lai = zeros(1,noNodes+1); % thickness of active layer for previous step
% A = zeros(noGrainSizeGroups,noNodes+1); %
% B = zeros(noGrainSizeGroups,noNodes+1); %
% qit = zeros(1,noNodes+1); %
% Chi = 1; % coefficient for calculating bedload layer sediment fraction
% taog = ones(1,noNodes+1);
% dSfdx = ones(1,noNodes+1);
% dSxdx = ones(1,noNodes+1);
% 
% GSGr = [-5 -4.5 -4 -3.5 -3 -2.5 -2 -1.5 -1 -0.5 0]; % grain size group
% GSGf = [0.6 1.7 4.4 9.2 15 38.2 15 9.2 4.4 1.7 0.6]/100; % grain size group fraction
% % GSGf = [1/11 1/11 1/11 1/11 1/11 1/11 1/11 1/11 1/11 1/11 1/11]; % grain size group fraction
% Phi = zeros(11,noNodes+1); % Grain size group
% F = zeros(11,noNodes+1); % fraction for each grain size group in active layer
% f = zeros(11,noNodes+1); % fraction for each grain size group in bedload layer
% for  ii=1:11
% Phi(ii,:) = GSGr(ii);
% F(ii,:) = GSGf(ii);
% f(ii,:) = GSGf(ii);
% end
% Fs = F; % fraction each grain size group in substrate
% Phia = sum(Phi(:,:).*F(:,:));
% PhiaFs = sum(Phi(:,:).*Fs(:,:));
% D(:,:) = 2.^Phi(:,:); % calculate mean grain size for each group in mm.
% Dg = 2.^Phia/1000; % geometric mean grain size 
% DgFs = 2.^PhiaFs/1000; % geometric mean grain size

% backwater calculation
for tspan=1:noTimeSteps
    M=noNodes+1;
    for i=1:(noNodes+1) 

    if M<noNodes+1
        kh1 = ((-be(M+1)+be(M))/dxn/(stf-sba)...
        -(kc^(1/3)*qw^2/ar^2/g/H(M+1)^(10/3)))/(1-qw^2/g/H(M+1)^3)*dxn*(stf-sba);      
        kh2 = ((-be(M+1)+be(M))/dxn/(stf-sba)...
        -(kc^(1/3)*qw^2/ar^2/g/(H(M+1)-kh1/2)^(10/3)))/(1-qw^2/g/(H(M+1)-kh1/2)^3)*dxn*(stf-sba);        
        kh3 = ((-be(M+1)+be(M))/dxn/(stf-sba)...
        -(kc^(1/3)*qw^2/ar^2/g/(H(M+1)-kh2/2)^(10/3)))/(1-qw^2/g/(H(M+1)-kh2/2)^3)*dxn*(stf-sba);  
        kh4 = ((-be(M+1)+be(M))/dxn/(stf-sba)...
        -(kc^(1/3)*qw^2/ar^2/g/(H(M+1)-kh3)^(10/3)))/(1-qw^2/g/(H(M+1)-kh3)^3)*dxn*(stf-sba);  
        H(M) = H(M+1)-kh1/6-kh2/3-kh3/3-kh4/6;         
        if H(M)<=Hn            
           H(M)=Hn;
        end
    end
    qs(M)=(R*g*D50^3)^0.5/(R*g*D50)^3*0.0355/ar^4*(kc/H(M))^(2/3)*(qw/H(M))^6;

%     % grain size calculation
%     taog(M) = 1/(ar^2)*(kc/H(M))^(1/3)*qw^2/R/g/Dg(M)/H(M)^2;
%     % calculate coefficient A and B
%     A(:,M) = 0.455*(D(:,M)/Dg(M)).^(-0.839);
%     B(:,M) = 0.353*(D(:,M)/Dg(M)).^(-1.157);
%     qit = sum(F(:,M).*A(:,M).*(taog(M)*Dg(M)./D(:,M)).^B(:,M));
% 
%     % calculate fraction of grain size in bedload
%     f(:,M) = F(:,M).*A(:,M).*(taog(M)*Dg(M)./D(:,M)).^B(:,M)/qit;    
% 
%     % calculate deviation
%     Dv(M) = Dv(M)+sum((Phi(:,M)-Phia(M)).^2 .* F(:,M));
%     SDv(M) = sqrt(Dv(M)); % standard deviation
%     La(M) = 2*8*H(M)*(Dg(M)/H(M))^0.3; %Dg(M) * 2^(1.28*SDv(M)); % thickness of active layer
% 
%     % calculate fraction in bedload layer
%     Fs(:,M) = Chi * F(:,M) + (1-Chi) * f(:,M);
%     Phia(M) = Phia(M)+ sum(Phi(:,M).*F(:,M));
%     PhiaFs(M) = PhiaFs(M)+ sum(Phi(:,M).*Fs(:,M));
%     D(:,M) = 2.^Phi(:,M); % calculate mean grain size for each group in mm.
%     Dg(M) = 2^Phia(M)/1000; % geometric mean grain size 
%     DgFs(M) = 2^PhiaFs(M)/1000; % geometric mean grain size

    % end of grain size calculation        
    M=M-1;
    end

    if tspan==1
    dqdxn(end) = (qs(end)-qs(end-1))/dxn;
    dqdxn(1) = (qs(2)-qs(1))/dxn;
    dndt(end) = 1e-5*(-Sfi)-dqdxn(end)*If*(1+LAMBDA)*Sinu/rB/(1-lambda)/(stf-sba);
    dndt(1) = 1e-5*(-Sfi)-dqdxn(1)*If*(1+LAMBDA)*Sinu/rB/(1-lambda)/(stf-sba);
    end
    rstf = 1/Sfore*(If*Sinu*(1+LAMBDA)*qs(end)/rB/(1-lambda)/(ssb-stf)-dndt(end)); % rate of migration of shoreline break
    rsba = -1/Sbb*dndt(1); % rate of migration of bedrock/alluvial contact
    rssb = 1/(Sfore-Ssb)*(Sfore*rstf+dndt(end)); % rate of rate of migration of foreset/bottom set contact
%     rstf = 0;
%     rsba = 0;
%     rssb = 0;
    
    Sxp = Sx(end); % Sx at downstream end from previous iteration, for sediment partition calculation 
    bet = be(end); % temporary be(end) for mass balance calculation

    % downstream point
    dqdxn(noNodes+1) = (qs(noNodes+1)-qs(noNodes))/dxn;
    dqdx(noNodes+1) = (qs(noNodes+1)-qs(noNodes))/dx;
    dndt(noNodes+1) = (xn(noNodes+1)*rstf+(1-xn(noNodes+1))*rsba)/(stf-sba)*(be(noNodes+1)-be(noNodes))/dxn...
    -dqdxn(noNodes+1)*If*(1+LAMBDA)*Sinu/rB/(1-lambda)/(stf-sba);
%     dFdxn(:,noNodes+1) = (F(:,noNodes+1)-F(:,noNodes))/dxn;
%     dLadxn(noNodes+1) = (La(noNodes+1)-La(noNodes))/dxn;
    Sx(noNodes+1) = -(be(noNodes+1)-be(noNodes))/(x(noNodes+1)-x(noNodes));
    dHdx(noNodes+1) = (H(noNodes+1)-H(noNodes))/(x(noNodes+1)-x(noNodes));

    % upstream point
    dqdxn(1) = (qs(2)-qs(1))/dxn;
    dqdx(1) = (qs(2)-qs(1))/dx;
    dndt(1) = (xn(1)*rstf+(1-xn(1))*rsba)/(stf-sba)*(be(2)-be(1))/dxn...
    -dqdxn(1)*If*(1+LAMBDA)*Sinu/rB/(1-lambda)/(stf-sba);
%     dFdxn(:,1) = (F(:,2)-F(:,1))/dxn;%dFdxn(ii,M) = (F(ii,M)-GSGf(ii))/dxn;
%     dLadxn(1) = (La(2)-La(1))/dxn;
    Sx(1) = -(be(2)-be(1))/(x(2)-x(1));
    dHdx(1) = (H(2)-H(1))/(x(2)-x(1));

    % cases 2 to N points
    dqdxn(2:noNodes) = au*(qs(2:noNodes)-qs(1:noNodes-1))/dxn+(1-au)*(qs(3:noNodes+1)-qs(2:noNodes))/dxn;
    dqdx(2:noNodes) = au*(qs(2:noNodes)-qs(1:noNodes-1))/dx+(1-au)*(qs(3:noNodes+1)-qs(2:noNodes))/dx;
    dndt(2:noNodes) = (xn(2:noNodes)*rstf+(1-xn(2:noNodes))*rsba)/(stf-sba).*(be(3:noNodes+1)-be(2:noNodes))/dxn...
    -dqdxn(2:noNodes)*If*(1+LAMBDA)*Sinu/rB/(1-lambda)/(stf-sba);
%     dFdxn(:,2:noNodes) = au*(F(:,2:noNodes)-F(:,1:noNodes-1))/dxn+(1-au)*(F(:,3:noNodes+1)-F(:,2:noNodes))/dxn;
%     dLadxn(2:noNodes) = au*(La(2:noNodes)-La(1:noNodes-1))/dxn+(1-au)*(La(3:noNodes+1)-La(2:noNodes))/dxn;
    Sx(2:noNodes) = -(be(2:noNodes)-be(1:noNodes-1))./(x(2:noNodes)-x(1:noNodes-1));
    dHdx(2:noNodes) = (H(2:noNodes)-H(1:noNodes-1))./(x(2:noNodes)-x(1:noNodes-1));
    %      keyboard   

    be = be+dndt*dt;
    sp = be+H;
    Sfrix = kc^(1/3)/ar^2*qw^2/g./H.^(10/3);
    dSfdx(noNodes+1) = (Sfrix(noNodes+1)-Sfrix(noNodes))/(x(noNodes+1)-x(noNodes));
    dSfdx(1) = (Sfrix(2)-Sfrix(1))/(x(2)-x(1));
    dSfdx(2:noNodes) = (Sfrix(2:noNodes)-Sfrix(1:noNodes-1))./(x(2:noNodes)-x(1:noNodes-1));
    dSxdx(noNodes+1) = (Sx(noNodes+1)-Sx(noNodes))/(x(noNodes+1)-x(noNodes));
    dSxdx(1) = (Sx(2)-Sx(1))/(x(2)-x(1));
    dSxdx(2:noNodes) = (Sx(2:noNodes)-Sx(1:noNodes-1))./(x(2:noNodes)-x(1:noNodes-1));

%     % calculate updated fractions of grain size in active layer
%     M=noNodes+1;
%     for i=1:(noNodes+1)
%         F(:,M) = F(:,M)+(xn(M)*rstf+(1-xn(M))*rsba)/(stf-sba)*dFdxn(M)*dt...
%         -1/La(M)*(F(:,M)-Fs(:,M))*((La(M)-Lai(M))/dt...
%         -(xn(M)*rstf+(1-xn(M))*rsba)/(stf-sba)*dLadxn(M))*dt...
%         -If*(1+LAMBDA)*Sinu/(La(M)*(1-lambda)*(stf-sba)*rB)*(dqdxn(M)*f(:,M)-Fs(:,M)*dqdxn(M))*dt;
%         F(F<0) = eps;
%         M=M-1;
%     end
%     Lai=La;% update old active layer Lai to new La for next cycle
%     Phia=zeros(noNodes+1); % reset the arithmetic mean to be zero for calculation of next node.
%     PhiaFs=zeros(noNodes+1);
%     Dv=zeros(noNodes+1); % reset the deviation to be zero for calculation of next node.
%     DgPhi = log2(1000*Dg);

    seaLevel = seaLevel+seaLevelRate; % update seaLevel
    H(end)=seaLevel-be(end);% update river mouth flow depth

    ssbt = ssb; % temporary ssb for mass balance calculation
    stft = stf; % temporary stf for mass balance calculation
    bebaset = bebase; % temporary bebase for mass balance calculation
    ssb = ssb+rssb*dt; % update ssb
    stf = stf+rstf*dt; % update stf
    sba = sba+rsba*dt; % update sba
    STF(tspan) = stf;

    x = xn.*(stf-sba)+sba; % convert normalized x distance back to dimensional
    dx = dxn*(stf-sba);
%     dDgdx = [(-Dg(2)+Dg(1))/(x(2)-x(1)),(Dg(1:noNodes-1)-Dg(3:noNodes+1))/(x(2)-x(1))/2,(Dg(noNodes)-Dg(noNodes+1))/(x(2)-x(1))];

    % calculate dSdx        
    dSdx(2:noNodes+1) = (Sx(2:noNodes+1)-Sx(1:noNodes))/dx;
    dSdx(1) = dSdx(2);

    % calculate bed elevation increment
    if rstf>0 %case of shoreline progradation
        BEInt = interp1(x,be,xp,'spline');  
        aggradationRate = BEInt-bep;% unit depends on tf value. m/yr if tf = 1, 
    else %case of shoreline retreat
        BEInt = interp1(xp,bep,x,'spline');
        aggradationRate = be - BEInt;
    end
    bep = be; % update BEP so for next iteration of calculation
    xp = x;
    
    if tspan/100 - round(tspan/100) == 0
        figure(1)
        hold on
        plot(x./1e3,be+H,'-b')
        plot(x./1e3,be,'-k')
        plot([stf/1e3,stf/1e3+100],[seaLevel,seaLevel],'-b')
    end
    
    migRate(tspan) = rstf;%shoreline migration rate
end

figure(1)
hold on
plot(x./1e3,be+H,'-b')
plot([stf/1e3,stf/1e3+100],[seaLevel,seaLevel],'-b')
xlim([-100,100]);
xlabel('distance (km)');
ylabel('elevation (m)');

toc

