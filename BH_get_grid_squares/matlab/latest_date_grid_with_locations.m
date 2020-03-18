clear
close
load inputs/colocbh_no.txt
load inputs/BH_locats.txt
load inputs/bhlocats.txt

set(0,'DefaultFigureColor',[1 1 1])

make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.05], [0.1 0.01], [0.1 0.01]);
if ~make_it_tight,  clear subplot;  end

nparts=173;
bh_long = (colocbh_no(:,1)-0.5)*5;
bh_lat = 90 -  (colocbh_no(:,2)-0.5)*5;
%bh_long2 = (colocbh_no2(:,1)-0.5)*5;
%bh_lat2 = 90 -  (colocbh_no2(:,2)-0.5)*5;

nparts=nparts;

for ic=1:nparts
 if(bh_long(ic)>180)
     bh_long(ic)=bh_long(ic) -360;
 end
end


% The number of MCMC iterations to sample from:
    itmax=50000;

hold on
for ic=1:nparts
    pos = [bh_long(ic)-2.5 bh_lat(ic)-2.5 5 5];
    rectangle('position',pos)
end

ax0=subplot(2,2,1)
hold all



xlim([-150 160])
ylim([-50 77.5])
box on
pbaspect([2 1 1])

hTitle=title('Borehole density')
hXLabel = xlabel('lon'                     );
hYLabel = ylabel('lat'                      );

set([hTitle, hXLabel, hYLabel], ...
    'FontName'   , 'Helvetica');









set(gca,'color',[0.8 0.8 0.8])



nparts=173;

partition_lonlat=colocbh_no;
partition_lonlat(:,1)=(colocbh_no(:,1))*5 - 2.5 ;
partition_lonlat(:,2)=(18 - colocbh_no(:,2)+1)*5- 2.5;
bhlocats(:,1) = bhlocats(:,1)*100*5 - 2.5 ;
bhlocats(:,2) = 90 - bhlocats(:,2)*37*5 +2.5;
 
glob_partition_lonlat=zeros(nparts,2);
glob_partition_lonlat(1:nparts,:)=partition_lonlat;


bh_long = (colocbh_no(:,1)-0.5)*5;
bh_lat = 90 -  (colocbh_no(:,2)-0.5)*5;

for ic=1:nparts
 if(bh_long(ic)>180)
     bh_long(ic)=bh_long(ic) -360;
 end
 if(partition_lonlat(ic,1)>180)
     partition_lonlat(ic,1)=partition_lonlat(ic,1) -360;
 end
 if(bhlocats(ic,1)>180)
     bhlocats(ic,1)=bhlocats(ic,1) -360;
 end
end

partition_count=zeros(nparts);
dellon=2.5;
dellat=2.5;
for ib=1:1012
for ic=1:nparts
   if(BH_locats(ib,1)>=partition_lonlat(ic,1)-dellon & BH_locats(ib,1)<partition_lonlat(ic,1)+dellon)
       if(BH_locats(ib,2)>=partition_lonlat(ic,2)-dellat & BH_locats(ib,2)<partition_lonlat(ic,2)+dellat)
           
           partition_count(ic) =partition_count(ic)+1;
       end 
   end
    
end
end

partition_count(:,1);
cmph=colormap(ax0, parula(7));
%
hold all
for ic=1:nparts
    pos = [partition_lonlat(ic,1)-2.5 partition_lonlat(ic,2)-2.5 5 5];
    
    if( partition_count(ic,1)<=2)
        color = cmph(1,:);
    end
    if(partition_count(ic,1)>2 & partition_count(ic,1)<=4)
        color = cmph(2,:);
    end
    if(partition_count(ic,1)>4 & partition_count(ic,1)<=6)
        color = cmph(3,:);
    end
    if(partition_count(ic,1)>6 & partition_count(ic,1)<=8)
        color = cmph(4,:);
    end
   if(partition_count(ic,1)>8 & partition_count(ic,1)<=10)
        color = cmph(5,:);
    end
    if(partition_count(ic,1)>10 & partition_count(ic,1)<=12)
        color = cmph(6,:);
    end
    if(partition_count(ic,1)>12)
        color = cmph(7,:);
    end
    %rectangle('position',pos,'FaceColor',color,'EdgeColor',[0.5 0.5 0.5],'LineWidth',0.1)
    rectangle('position',pos,'FaceColor',color,'EdgeColor','none')
    
end

hcb1=colorbar;
set(gca, 'CLim', [0, 14]);
set(hcb1,'YTick',[,2,4,6,8,10,12,],'Location','eastoutside')
load inputs/world.dat
plot(world(:,1),world(:,2),'color','black')

%plot(bh_long,bh_lat,'linestyle','none','marker','s','markersize',15,'color','blue')
%plot(BH_locats(:,1),BH_locats(:,2),'linestyle','none','marker','.','markersize',8,'color','black','linewidth',0.25)
%hold on
%plot(BH_locats2(:,1),BH_locats2(:,2),'linestyle','none','marker','o','markersize',3,'color','black','linewidth',1)


% ------------------------------------
ax1=subplot(2,2,2)

% Now I need the latest logging date in each gridcell:
% This is the logging date minus 10 years in each case
borehole_date = zeros(nparts,2);
borehole_z0 = zeros(nparts,2);
borehole_zmax = zeros(nparts,2);
load inputs/minsofar.dat
borehole_date(1:nparts,2) = minsofar;

load inputs/minsofar2.dat
borehole_z0(1:nparts,2)= minsofar2;


load inputs/zmax.dat
borehole_zmax(1:nparts,2)= zmax;


hold all
load inputs/world.dat
plot(world(:,1),world(:,2),'color','black')

xlim([-180 180])
ylim([-50 80])
box on
pbaspect([2 2 1])
cmph=colormap(ax1, parula(6));
%
hold all
for ic=1:nparts
    pos = [partition_lonlat(ic,1)-2.5 partition_lonlat(ic,2)-2.5 5 5];
    
    if( borehole_date(ic,2)<=1965)
        color = cmph(1,:);
    end
    if(borehole_date(ic,2)>1965 & borehole_date(ic,2)<=1975)
        color = cmph(2,:);
    end
    if(borehole_date(ic,2)>1975 & borehole_date(ic,2)<=1985)
        color = cmph(3,:);
    end
    if(borehole_date(ic,2)>1985 & borehole_date(ic,2)<=1995)
        color = cmph(4,:);
    end
   if(borehole_date(ic,2)>1995 & borehole_date(ic,2)<=2005)
        color = cmph(5,:);
    end
    if(borehole_date(ic,2)>2005 & borehole_date(ic,2)<=2015)
        color = cmph(6,:);
    end
    %rectangle('position',pos,'FaceColor',color,'EdgeColor',[0.5 0.5 0.5],'LineWidth',0.1)
    rectangle('position',pos,'FaceColor',color,'EdgeColor','none')
    
end

hcb2=colorbar;
set(gca, 'CLim', [1955, 2015]);
set(hcb2,'YTick',[,1965,1975,1985,1995,2005,2015],'Location','eastoutside')
%cbh=colorbar(gca);
%set(gca,'XTick',[1955:10:2015],'XTickLabel',{'','1965','1975','1985','1995','2005',''});

xlim([-150 160])
ylim([-50 77.5])
box on
pbaspect([2 1 1])
hTitle=title('Latest logging date (year CE)')
hXLabel = xlabel('lon'                     );
hYLabel = ylabel('lat'                      );

set([hTitle, hXLabel, hYLabel], ...
    'FontName'   , 'Helvetica');


set(gca,'color',[0.8 0.8 0.8])








ax2=subplot(2,2,3)
hold all
load inputs/world.dat
plot(world(:,1),world(:,2),'color','black')

xlim([-180 180])
ylim([-50 80])
box on
pbaspect([2 1 1])
%cmp=colormap(ax2,hot(11        ))
%cmp=colormap(ax2,ametrine(11,'invert',1))
cmp=colormap(ax2,parula(8))
hcb3=colorbar;
hold all
for ic=1:nparts
    pos = [partition_lonlat(ic,1)-2.5 partition_lonlat(ic,2)-2.5 5 5];
    
    if( borehole_z0(ic,2)<=10)
        color = cmp(1,:);
    end
    if(borehole_z0(ic,2)>10 & borehole_z0(ic,2)<=20)
        color = cmp(2,:);
    end
    if(borehole_z0(ic,2)>20 & borehole_z0(ic,2)<=30)
        color = cmp(3,:);
    end
    if(borehole_z0(ic,2)>30 & borehole_z0(ic,2)<=40)
        color = cmp(4,:);
    end
   if(borehole_z0(ic,2)>40 & borehole_z0(ic,2)<=50)
        color = cmp(5,:);
    end
    if(borehole_z0(ic,2)>50 & borehole_z0(ic,2)<=60)
        color = cmp(6,:);
    end
    if(borehole_z0(ic,2)>60 & borehole_z0(ic,2)<=70)
        color = cmp(7,:);
    end
    if(borehole_z0(ic,2)>70)
        color = cmp(8,:);
    end
   
    rectangle('position',pos,'FaceColor',color,'EdgeColor','none')
    
end


set(gca, 'CLim', [0, 80]);
set(hcb3,'YTick',[,10,20,30,40,50,60,70,],'Location','eastoutside')
xlim([-150 160])
ylim([-50 77.5])
box on
pbaspect([2 1 1])
hTitle=title('Profile top (m)')
hXLabel = xlabel('lon'                     );
hYLabel = ylabel('lat'                      );

set([hTitle, hXLabel, hYLabel], ...
    'FontName'   , 'Helvetica');
set(gca,'color',[0.8 0.8 0.8])









ax3=subplot(2,2,4)
hold all
load inputs/world.dat
plot(world(:,1),world(:,2),'color','black')

xlim([-180 180])
ylim([-50 80])
box on
pbaspect([2 1 1])
cmp3=colormap(ax3,parula(7        ))
hcb4=colorbar
hold all
for ic=1:nparts
    pos = [partition_lonlat(ic,1)-2.5 partition_lonlat(ic,2)-2.5 5 5];
    
    if( borehole_zmax(ic,2)<=250)
        color = cmp3(1,:);
    end
    if(borehole_zmax(ic,2)>250 & borehole_zmax(ic,2)<=350)
        color = cmp3(2,:);
    end
    if(borehole_zmax(ic,2)>350 & borehole_zmax(ic,2)<=450)
        color = cmp3(3,:);
    end
    if(borehole_zmax(ic,2)>450 & borehole_zmax(ic,2)<=550)
        color = cmp3(4,:);
    end
   if(borehole_zmax(ic,2)>550 & borehole_zmax(ic,2)<=650)
        color = cmp3(5,:);

    end
    if(borehole_zmax(ic,2)>650 & borehole_zmax(ic,2)<=750)
        color = cmp3(6,:);
    end
   
    
    if(borehole_zmax(ic,2)>750)
        color = cmp3(7,:);
    end
    rectangle('position',pos,'FaceColor',color,'EdgeColor','none')
    
end


set(gca, 'CLim', [150, 850]);
set(hcb4,'YTick',[,250,350,450,550,650,750,],'Location','eastoutside')
xlim([-150 160])
ylim([-50 77.5])
box on
pbaspect([2 1 1.5])
hTitle=title('Profile bottom (m)')
hXLabel = xlabel('lon'                     );
hYLabel = ylabel('lat'                      );

set([hTitle, hXLabel, hYLabel], ...
    'FontName'   , 'Helvetica');
set(gca,'color',[0.8 0.8 0.8])

ax0.TitleFontSizeMultiplier = 0.9;
ax1.TitleFontSizeMultiplier = 0.9;
ax2.TitleFontSizeMultiplier = 0.9;
ax3.TitleFontSizeMultiplier = 0.9;


x0=100;
y0=100;
width=900;
height=500

set(gcf,'position',[x0,y0,width,height])
fig = gcf;
fig.InvertHardcopy = 'off';
print -painters -depsc2 -r2500 plots/logging_date_depths_density_all.eps

%saveas(gca,'plots/myfigure.pdf')