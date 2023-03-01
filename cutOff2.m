% calculate cutoff for the last (noNode-k-5)points
function [x,y,T1,noNode] = cutOff2(noNode,k,x,y,b,t,T1)
    for i = noNode-k-5:(noNode-5)
        if i>noNode-3
            break
        end
        T1(t)=T1(t);

        p = [x(i+3:end);y(i+3:end)]';
        pq = [x(i),y(i)];

        [idx,dist] = dsearchn(p,pq);
        if dist < 2*b
    %             figure(1);hold on
    %             subplot(2,1,1),hold on
    %             plot(x(i:i+2+idx),y(i:i+2+idx),'-','color',[.6 .8 1],'linewidth',1);
    %             plot(x(i)/1e3,y(i)/1e3,'or','markerfacecolor','k','markersize',1.5)
    %             plot(x(i+4+idx)/1e3,y(i+4+idx)/1e3,'or','markerfacecolor','r','markersize',1.5)

            x = [x(1:i),x(i+4+idx:end)];
            y = [y(1:i),y(i+4+idx:end)];
            T1(t) = t;   
    %             keyboard
            clearvars p
        end

        [~,noNode] = size(x);

    end
end