clear
load('progradationCaseResults.mat')
plotTime = [1,500,1000];
plotOffSet = [6,0,-8];

xRange = max(X,[],'all');
x = X-5000;
y = Y;
x = x/xRange;%normalize x distance by valley length
x = 1 - x;%flip x so that 0 is mouth and 1 is upstream
y = y/xRange;%normalize
xOver = XOver/xRange;
xOver = 1 - xOver;
segment = nan(3,1000);

for i = 1:1000
    temp = SINU(i,:);
    segment(1,i) = mean(temp(xOver(i,:)<=0.33));
    segment(2,i) = mean(temp(xOver(i,:)>0.33 & xOver(i,:)<=0.67));
    segment(3,i) = mean(temp(xOver(i,:)>0.67 & xOver(i,:)<=1.00));
end


X(X==0)=nan;%for some reason, there are 0 values
figure(1)
ax1 = axes;
set(ax1,'Position',[0.1 0.62 .7 .35]);

mR = abs(migrRate_t)/20*30-15;
h = fill(x_sort/1e3-5,mR,'r');
set(h,'FaceColor',[.8 .8 .8],'EdgeColor','none');
xlim([0 80])

hold on

for j = 1:3
    t = plotTime(j);    
    x = X-5000;
    y = Y;
    offSet = plotOffSet(j);
    if  j ==3
        ph = plot(coXY(1,:)/1e3-5,coXY(2,:)/1e3 + offSet,'-','color',[.6 .6 .6]);
    end
    p(j) = plot(x(t,:)/1e3,y(t,:)/1e3 + offSet,'-k','linewidth',1);
    ylim([-15000/1e3 15000/1e3]);
    txt = sprintf('t = %d yr',plotTime(j));
    if j == 3
        text(64,-5,txt);
    elseif j == 2
        text(58,plotOffSet(j)+2,txt);
    else
        text(23,plotOffSet(j)+2.5,txt);
    end
   
    set(gca,'DataAspectRatio', [1 1 1]);  
    
    if j == 3
        plot([0 80],[15 15],'-k','LineWidth',1)
        plot([80 80],[-15 15],'-k','LineWidth',1);

        text(5,14,'A','VerticalAlignment','top','HorizontalAlignment','center','fontsize',10)

        %plot right axis
        plot([79.1 80],[0 0],'-k','LineWidth',.5)
        text((81),(-15),'0','VerticalAlignment','middle','HorizontalAlignment','left','fontsize',9,'color',[.7 .7 .7])
        text((81),(0),'10','VerticalAlignment','middle','HorizontalAlignment','left','fontsize',9,'color',[.7 .7 .7])
        text((81),(15),'20','VerticalAlignment','middle','HorizontalAlignment','left','fontsize',9,'color',[.7 .7 .7])
        h = text((88),(-8),{'lateral migration';'    rate (m/yr)'},'VerticalAlignment','middle','HorizontalAlignment','left','fontsize',10);
        set(h,'Rotation',90,'color',[.7 .7 .7]);

        xlabel('{\itx} (km)','fontsize',10);
        ylabel('{\ity} (km)','fontsize',10);
        set(0,'DefaultAxesTitleFontWeight','normal');
    end
    set(gca,'layer','top','XColor',[0 0 0],'YColor',[0 0 0],'LineWidth',1);

    
    clearvars aveSinu2 lowerLimit upperLimit pos
end
box off
legend ([p(3),ph],{'channel','abandoned channel'},'box','off');

