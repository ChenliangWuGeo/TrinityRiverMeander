function output = calculateSinu2(x_t,y_t,accuS_t,b)
    x = x_t;
    y = y_t;
    accuS = accuS_t;
    b2 = 2*b;
    
    curvature = getCurv(x,y);

    
    [pks1, loca1]=findpeaks(curvature);
    [pks2, loca2]=findpeaks(-curvature);

    locaAll = [loca1,loca2];
    locaAll = sort(locaAll);

    [~,lgth] = size(locaAll);
    if mod(lgth,2) ~= 0 %remainder after divided by 2, in case of odd number
        locaAll = locaAll(1:end-1);
        lgth = lgth - 1;
    end

    xOver = nan(1,lgth);
    yOver = nan(1,lgth);
    accuSOver = nan(1,lgth);

    intrinsicLength = nan(1,lgth-1);
    directLength = nan(1,lgth-1);

    for i = 1  : lgth-1
        locaTemp = locaAll(i):locaAll(i+1);
        xTemp = x(locaTemp);
        yTemp = y(locaTemp);
        accuSTemp = accuS(locaTemp);
        curvTemp = curvature(locaTemp);
        xOver(i) = interp1(curvTemp,xTemp,[0]);
        yOver(i) = interp1(curvTemp,yTemp,[0]);
        accuSOver(i) = interp1(curvTemp,accuSTemp,[0]);
    end

    xOver = xOver(~isnan(xOver));
    yOver = yOver(~isnan(yOver));
    accuSOver = accuSOver(~isnan(accuSOver));
    RK = max(accuS)-accuS;
    RKOver = max(accuS)-accuSOver;
    
    intrinsicLength = accuSOver(2:end)-accuSOver(1:end-1);
    directLength = sqrt(...
        (xOver(2:end)-xOver(1:end-1)).^2 +...
        (yOver(2:end)-yOver(1:end-1)).^2);
   
    %to ensure local perturbations excluded.
    accuSOver(intrinsicLength<b2*3) = [];
    xOver(intrinsicLength<b2*3) = [];
    yOver(intrinsicLength<b2*3) = [];
    RKOver(intrinsicLength<b2*3) = [];
    intrinsicLength = accuSOver(2:end)-accuSOver(1:end-1);
    directLength = sqrt(...
        (xOver(2:end)-xOver(1:end-1)).^2 +...
        (yOver(2:end)-yOver(1:end-1)).^2);

    sinu2 = intrinsicLength./directLength;
    aveSinu = movmean(sinu2,10);
    xtemp = [0:0.05:3];
%     aveSinu = interp1(fliplr(RKOver(2:end)/river.Lb/1e3),fliplr(aveSinu),...
%         xtemp);
    
    %calculate sinuosity using moving window
    n = 20;
    Lpath = n * b2*1e3/2;
    Ldist = sqrt( (x(n+1:end) - x(1:end-n)).^2+...
        (y(n+1:end)- y(1:end-n)).^2);
    sinu1 = Lpath./Ldist;
    sinu1RK = RK(n+1:end);
    
    %store output
    output.sinu1 = sinu1;
    output.sinu2 = sinu2;
    output.sinu1RK = sinu1RK;
    output.xOver = xOver;
    output.yOver = yOver;
    output.accuSOver = accuSOver;
    output.xy = [x;y];
    output.curvature = curvature;
    output.accuS = accuS;
    output.RK = RK;
    output.RKOver = RKOver;
%     output.Lb = river.Lb;
%     output.xyi = [xi;yi];
%     output.aveSinu = aveSinu;
    
    % figure(1);hold on
    % plot(x/1e3,y/1e3,'-k','linewidth',1.5)
    % set(gca,'DataAspectRatio', [1 1 1]);
    % ylabel('km');
    % xlabel('km');
    % % plot(x_t(locaAll)/1e3,y_t(locaAll)/1e3,'ro')
    % plot(xOver/1e3,yOver/1e3,'ob')
end

function output = getCurv(x,y)
%     [~,temp] = size(x);   
    dx1 = x;%nan(1,temp);
    dx2 = x;%nan(1,temp);
    dy1 = x;%nan(1,temp);
    dy2 = x;%nan(1,temp);
    
    ds = x;
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
    
    output = curvature;
end
% %this is the function version of findCrossOver.m
% function output = calculateSinu2(curvature_t,x_t,y_t,accuS_t)
% 
%     [pks1, loca1]=findpeaks(curvature_t);
%     [pks2, loca2]=findpeaks(-curvature_t);
% 
%     locaAll = [loca1,loca2];
%     locaAll = sort(locaAll);
% 
%     [~,lgth] = size(locaAll);
%     if mod(lgth,2) ~= 0 %remainder after divided by 2, in case of odd number
%         locaAll = locaAll(1:end-1);
%         lgth = lgth - 1;
%     end
% 
%     xOver = nan(1,lgth);
%     yOver = nan(1,lgth);
%     accuSOver = nan(1,lgth);
% 
%     intrinsicLength = nan(1,lgth-1);
%     directLength = nan(1,lgth-1);
% 
%     for i = 1  : lgth-1
%         locaTemp = locaAll(i):locaAll(i+1);
%         xTemp = x_t(locaTemp);
%         yTemp = y_t(locaTemp);
%         accuSTemp = accuS_t(locaTemp);
%         curvTemp = curvature_t(locaTemp);
%         xOver(i) = interp1(curvTemp,xTemp,[0]);
%         yOver(i) = interp1(curvTemp,yTemp,[0]);
%         accuSOver(i) = interp1(curvTemp,accuSTemp,[0]);
%     end
% 
%     xOver = xOver(~isnan(xOver));
%     yOver = yOver(~isnan(yOver));
%     accuSOver = accuSOver(~isnan(accuSOver));
% 
%     intrinsicLength = accuSOver(2:end)-accuSOver(1:end-1);
%     directLength = sqrt(...
%         (xOver(2:end)-xOver(1:end-1)).^2 +...
%         (yOver(2:end)-yOver(1:end-1)).^2);
% 
%     sinu2 = intrinsicLength./directLength;
%     
%     output = sinu2;
%     % figure(1);hold on
%     % plot(x/1e3,y/1e3,'-k','linewidth',1.5)
%     % set(gca,'DataAspectRatio', [1 1 1]);
%     % ylabel('km');
%     % xlabel('km');
%     % % plot(x_t(locaAll)/1e3,y_t(locaAll)/1e3,'ro')
%     % plot(xOver/1e3,yOver/1e3,'ob')
% end

