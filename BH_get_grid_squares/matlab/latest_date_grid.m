clear
close
load inputs/colocbh_no.txt
load inputs/BH_locats.txt
load inputs/bhlocats.txt

nparts=173;

partition_lonlat=colocbh_no;
partition_lonlat(:,1)=(colocbh_no(:,1))*5 - 2.5 ;
partition_lonlat(:,2)=(18 - colocbh_no(:,2)+1)*5- 2.5;
bhlocats(:,1) = bhlocats(:,1)*100*5 - 2.5 ;
bhlocats(:,2) = 90 - bhlocats(:,2)*37*5 +2.5;
 



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


glob_partition_lonlat=zeros(nparts,2);
glob_partition_lonlat(1:nparts,:)=partition_lonlat;


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


ax1=subplot(3,1,1)
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
    rectangle('position',pos,'FaceColor',color)
    
end

hcb1=colorbar;
set(gca, 'CLim', [1955, 2015]);

%cbh=colorbar(gca);
%set(gca,'XTick',[1955:10:2015],'XTickLabel',{'','1965','1975','1985','1995','2005',''});

xlim([-150 160])
ylim([-50 77.5])
box on
pbaspect([2 1 1])
title('Latest logging date (year AD)')













ax2=subplot(3,1,2)
hold all
load inputs/world.dat
plot(world(:,1),world(:,2),'color','black')

xlim([-180 180])
ylim([-50 80])
box on
pbaspect([2 1 1])
%cmp=colormap(ax2,hot(11        ))
%cmp=colormap(ax2,ametrine(11,'invert',1))
cmp=colormap(ax2,parula(11))
hcb2=colorbar;
hold all
for ic=1:nparts
    pos = [partition_lonlat(ic,1)-2.5 partition_lonlat(ic,2)-2.5 5 5];
    
    if( borehole_z0(ic,2)<=10)
        color = cmp(1,:);
    end
    if(borehole_z0(ic,2)>10 & borehole_z0(ic,2)<=15)
        color = cmp(2,:);
    end
    if(borehole_z0(ic,2)>15 & borehole_z0(ic,2)<=20)
        color = cmp(3,:);
    end
    if(borehole_z0(ic,2)>20 & borehole_z0(ic,2)<=25)
        color = cmp(4,:);
    end
   if(borehole_z0(ic,2)>25 & borehole_z0(ic,2)<=30)
        color = cmp(5,:);
    end
    if(borehole_z0(ic,2)>30 & borehole_z0(ic,2)<=35)
        color = cmp(6,:);
    end
    if(borehole_z0(ic,2)>35 & borehole_z0(ic,2)<=40)
        color = cmp(7,:);
    end
    if(borehole_z0(ic,2)>40 & borehole_z0(ic,2)<=45)
        color = cmp(8,:);
    end
    if(borehole_z0(ic,2)>45 & borehole_z0(ic,2)<=50)
        color = cmp(9,:);
    end
    if(borehole_z0(ic,2)>50 & borehole_z0(ic,2)<=55)
        color = cmp(10,:);
    end
    if(borehole_z0(ic,2)>60)
        color = cmp(11,:);
    end
    rectangle('position',pos,'FaceColor',color)
    
end


set(gca, 'CLim', [0, 65]);

xlim([-150 160])
ylim([-50 77.5])
box on
pbaspect([2 1 1])
title('Profile top (m)')










ax3=subplot(3,1,3)
hold all
load inputs/world.dat
plot(world(:,1),world(:,2),'color','black')

xlim([-180 180])
ylim([-50 80])
box on
pbaspect([2 1 1])
cmp3=colormap(ax3,parula(13        ))
colorbar
hold all
for ic=1:nparts
    pos = [partition_lonlat(ic,1)-2.5 partition_lonlat(ic,2)-2.5 5 5];
    
    if( borehole_zmax(ic,2)<=250)
        color = cmp3(1,:);
    end
    if(borehole_zmax(ic,2)>250 & borehole_zmax(ic,2)<=300)
        color = cmp3(2,:);
    end
    if(borehole_zmax(ic,2)>300 & borehole_zmax(ic,2)<=350)
        color = cmp3(3,:);
    end
    if(borehole_zmax(ic,2)>350 & borehole_zmax(ic,2)<=400)
        color = cmp3(4,:);
    end
   if(borehole_zmax(ic,2)>400 & borehole_zmax(ic,2)<=450)
        color = cmp3(5,:);
    end
    if(borehole_zmax(ic,2)>450 & borehole_zmax(ic,2)<=500)
        color = cmp3(6,:);
    end
    if(borehole_zmax(ic,2)>500 & borehole_zmax(ic,2)<=550)
        color = cmp3(7,:);
    end
    if(borehole_zmax(ic,2)>550 & borehole_zmax(ic,2)<=600)
        color = cmp3(8,:);
    end
    if(borehole_zmax(ic,2)>600 & borehole_zmax(ic,2)<=650)
        color = cmp3(9,:);
    end
    if(borehole_zmax(ic,2)>650 & borehole_zmax(ic,2)<=700)
        color = cmp3(10,:);
    end
    if(borehole_zmax(ic,2)>700 & borehole_zmax(ic,2)<=750)
        color = cmp3(11,:);
    end
    if(borehole_zmax(ic,2)>750 & borehole_zmax(ic,2)<=800)
        color = cmp3(12,:);
    end
    if(borehole_zmax(ic,2)>800)
        color = cmp3(13,:);
    end
    rectangle('position',pos,'FaceColor',color)
    
end


set(gca, 'CLim', [200, 850]);

xlim([-150 160])
ylim([-50 77.5])
box on
pbaspect([2 1 1.5])
title('Profile bottom (m)')

print -painters -depsc2 -r2500 plots/logging_date_depths_all.eps
