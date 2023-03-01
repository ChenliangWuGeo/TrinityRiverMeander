clear
load('EGradientResults.mat')
plotTime = [1,500,1000];
plotOffSet = [7000,1000,-7000]/1e3;
x = X;
y = Y;

figure(1)
ax2 = axes;
set(ax2,'Position',[0.1 0.22 .7 .35]);
set(gcf,'Unit','inch','Position',[1 1 5 5])

% yyaxis right
mR = abs(migrRate_t)/20*30-15;
h = fill(x_sort/1e3,mR,'r');
set(h,'FaceColor',[.8 .8 .8],'EdgeColor','none');
ylabel({'lateral migration';'rate (m/yr)'});

hold on

xlim([0 80]);
ylim([-15 15]);

for j = 1:3
    t = plotTime(j);
    offSet = plotOffSet(j);
    if  j ==3
        plot(coXY(1,:)/1e3-5,coXY(2,:)/1e3 + offSet,'-','color',[.6 .6 .6]);
    end

    plot(x(t,:)/1e3,y(t,:)/1e3 + offSet,'-k','linewidth',1);
    txt = sprintf('t = %d yr',plotTime(j));
    if j == 3
        text(64,plotOffSet(j)+3.5,txt);
    else
        text(64,plotOffSet(j)+2.5,txt);
    end
   
    set(gca,'DataAspectRatio', [1 1 1]);  
    
    if j == 3
        %plot upper axis
        plot([0 80],[15 15],'-k','LineWidth',1)
        plot([80 80],[-15 15],'-k','LineWidth',1);
        for n = 1:4
            m = 4 - n + 1;
            plot(([m m]*1/3-1/3)*80,[14.1 15],'-k','LineWidth',1)
            tempLabel = num2str((n)*2);
            text((m*1/3-1/3)*80,(15),tempLabel,'VerticalAlignment','bottom','HorizontalAlignment','center','fontsize',10)
        end
        
        %plot right axis
        plot([79.1 80],[0 0],'-k','LineWidth',.5)
        text((81),(-15),'0','VerticalAlignment','middle','HorizontalAlignment','left','fontsize',9,'color',[.7 .7 .7])
        text((81),(0),'10','VerticalAlignment','middle','HorizontalAlignment','left','fontsize',9,'color',[.7 .7 .7])
        text((81),(15),'20','VerticalAlignment','middle','HorizontalAlignment','left','fontsize',9,'color',[.7 .7 .7])
        h = text((88),(-8),{'lateral migration';'    rate (m/yr)'},'VerticalAlignment','middle','HorizontalAlignment','left','fontsize',10);
        set(h,'Rotation',90,'color',[.7 .7 .7]);

        
        text(40,17,'{\itE} (10^{-7})','VerticalAlignment','bottom','HorizontalAlignment','center','fontsize',9)
        
        text(5,14,'B','VerticalAlignment','top','HorizontalAlignment','center','fontsize',10)

        xlabel('{\itx} (km)','fontsize',10);
        ylabel('{\ity} (km)','fontsize',10);
        set(0,'DefaultAxesTitleFontWeight','normal');
    end

    set(gca,'layer','top','XColor',[0 0 0],'YColor',[0 0 0],'LineWidth',1);

    clearvars aveSinu2 lowerLimit upperLimit pos
end
box off


