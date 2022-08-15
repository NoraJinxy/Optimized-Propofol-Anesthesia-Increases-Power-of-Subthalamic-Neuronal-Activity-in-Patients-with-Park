% function draw_radar1(data1,data2,data3,data4,lim,labels)
function draw_radar1(data1,data2,lim,labels)

%     n=length(data1);
%     limit=[lim(1):0.2:lim(2)];
%     baseline=zeros(n,1)+1;
%     point_baseline=zeros(n,2);
%     Awake=data2./data1;
%     point_Awake=zeros(n,2);
%     TCI=data4./data3;
%     point_TCI=zeros(n,2);
%     A=120;
%     a=A/(numel(limit)-1);
    
    n=length(data1);
    limit=[lim(1):0.2:lim(2)];
    baseline=zeros(1,n)+1;
    point_baseline=zeros(n,2);
    Awake=data1+baseline;
    point_Awake=zeros(n,2);
    TCI=data2+baseline;
    point_TCI=zeros(n,2);
    A=120;
    a=A/(numel(limit)-1);

    set(gca,'units','normal','pos',[0 0 1 1]);
    axis off
    axis equal
    hold on
    theta_last=pi/2;
    for i=1:n
        theta=2*pi/n*i+pi/2;
        plot([0,A*cos(theta)],[0,A*sin(theta)],'k-','linewidth',2);
        for j=1:numel(limit)-1
           plot([j*a*cos(theta_last),j*a*cos(theta)],[j*a*sin(theta_last),j*a*sin(theta)],'--','linewidth',0.75,'color',[0.5,0.5,0.5]);
        end
        theta_last=theta;
        adj_baseline=(baseline(i)-lim(1))/(lim(2)-lim(1))*A;
        point_baseline(i,:)=[adj_baseline*cos(theta);adj_baseline*sin(theta)];
        adj_Awake=(Awake(i)-lim(1))/(lim(2)-lim(1))*A;
        point_Awake(i,:)=[adj_Awake*cos(theta);adj_Awake*sin(theta)];
        adj_TCI=(TCI(i)-lim(1))/(lim(2)-lim(1))*A;
        point_TCI(i,:)=[adj_TCI*cos(theta);adj_TCI*sin(theta)];
        text_around((A+10)*cos(theta),(A+10)*sin(theta),labels{i},theta);
    end
    plot([point_baseline(:,1);point_baseline(1,1)],[point_baseline(:,2);point_baseline(1,2)],'k-','linewidth',1.5);
    fill(point_baseline(:,1),point_baseline(:,2),[0 0 0.80392]);
    plot([point_Awake(:,1);point_Awake(1,1)],[point_Awake(:,2);point_Awake(1,2)],'k-','linewidth',1.5);
    fill(point_Awake(:,1),point_Awake(:,2),[0.52941 0.80784 0.92157]);
    alpha(0.5);
    plot([point_TCI(:,1);point_TCI(1,1)],[point_TCI(:,2);point_TCI(1,2)],'k-','linewidth',1.5);
    fill(point_TCI(:,1),point_TCI(:,2),[1 1 1])
    alpha(0.5);
    theta=2*pi/n+pi/2;
        for j=1:numel(limit)
            text_around((j-1)*a*cos(theta),(j-1)*a*sin(theta),num2str(limit(j)),theta+pi/2,7);
        end
end

function text_around(x,y,txt,theta,fontsize)
    if nargin==4
        fontsize=10;
    end
    section=mod(theta+pi/12,2*pi);
    if section>pi+pi/6
        %上对齐
        if section>1.5*pi+pi/6
            %左对齐
            text(x,y,txt,'VerticalAlignment','cap','HorizontalAlignment','left','Fontsize',fontsize);
        elseif section>1.5*pi
            %中对齐
            text(x,y,txt,'VerticalAlignment','cap','HorizontalAlignment','center','Fontsize',fontsize);
        else
            %右对齐
            text(x,y,txt,'VerticalAlignment','cap','HorizontalAlignment','right','Fontsize',fontsize);
        end
    elseif section>pi
        %中、右对齐
        text(x,y,txt,'VerticalAlignment','middle','HorizontalAlignment','right','Fontsize',fontsize);
    elseif section>pi/6
        %下对齐
        if section>0.5*pi+pi/6
            %右对齐
            text(x,y,txt,'VerticalAlignment','bottom','HorizontalAlignment','right','Fontsize',fontsize);
        elseif section>0.5*pi
            %中对齐
            text(x,y,txt,'VerticalAlignment','bottom','HorizontalAlignment','center','Fontsize',fontsize);
        else
            %左对齐
            text(x,y,txt,'VerticalAlignment','bottom','HorizontalAlignment','left','Fontsize',fontsize);
        end
    else
        %中、左对齐
        text(x,y,txt,'VerticalAlignment','middle','HorizontalAlignment','left','Fontsize',fontsize);
    end
end