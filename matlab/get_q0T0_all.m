
% ----------------------------------------------------------
% Now get the heatflow and Teq
% ----------------------------------------------------------

clear
set(0,'DefaultFigureColor',[1 1 1])

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


root='/Users/Peter/Documents/WORK/NH_SH_bh/glob_bh_new/Samples/';
extension='H.txt';
root2='X';

q0T0 = zeros(nparts,2);   % 104 partitions x 120 timesteps (5yrs)
%stdev_G = zeros(104,120);   % 104 partitions x 120 timesteps (5yrs)

for ic=1:nparts
      
      s= num2str(ic-1);
    filenameH=strcat(root,s,extension);
    varname = strcat(root2,s,'H')
    %q0T0 = textread(filenameH, '%f','delimiter', '\n');
    load(filenameH);
    varhere=eval(varname);
    sq0 = size(varhere(1,:))/2;
    sq0 = sq0(2);
    for is=1:sq0
      q0T0(ic,1) =q0T0(ic,1)+ (1.00/sq0)*mean(varhere(:,1+(is-1)*2));
      q0T0(ic,2) =q0T0(ic,2)+ (1.00/sq0)*mean(varhere(:,2+(is-1)*2));
    end

    
end

title('Site locations')
make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.05], [0.1 0.01], [0.1 0.01]);
if ~make_it_tight,  clear subplot;  end


ax0=subplot(1,2,1)
hold all
load inputs/world.dat
plot(world(:,1),world(:,2),'color','black')

xlim([-180 180])
ylim([7.5 77.5])
box on
pbaspect([2 1 1])
cmp=colormap(jet(7))

hold all
for ic=1:nparts
    pos = [partition_lonlat(ic,1)-2.5 partition_lonlat(ic,2)-2.5 5 5];
    acolor = (q0T0(ic,2) +10)/35
    bcolor = acolor;
    ccolor = bcolor;
    if( q0T0(ic,2)<=-5)
        color = cmp(1,:);
    end
    if(q0T0(ic,2)>-5 & q0T0(ic,2)<=0)
        color = cmp(2,:);
    end
    if(q0T0(ic,2)>0 & q0T0(ic,2)<=5)
        color = cmp(3,:);
    end
    if(q0T0(ic,2)>5 & q0T0(ic,2)<=10)
        color = cmp(4,:);
    end
    if(q0T0(ic,2)>10 & q0T0(ic,2)<=15)
        color = cmp(5,:);
    end
    if(q0T0(ic,2)>15 & q0T0(ic,2)<=20)
        color = cmp(6,:);
    end
    if(q0T0(ic,2)>20 )
        color = cmp(7,:);
    end
    
    rectangle('position',pos,'FaceColor',color,'EdgeColor','none')
    
end
colormap(gray(7))
cbh0=colorbar
set(gca, 'CLim', [-10, 25]);
set(cbh0,'YTick',[,-5,0,5,10,15,20,],'Location','eastoutside')
set(gca,'color',[0.8 0.8 0.8])
xlim([-150 160])
ylim([-50 77.5])
box on
pbaspect([2 1 1])
title('Pre-reconstruction mean temperature (\circC)')




% NOw heatflow q0
ax1=subplot(1,2,2)

plot(world(:,1),world(:,2),'color','black')

xlim([-180 180])
ylim([-50 77.5])
box on



hold all
colormap(jet(7))
cbh2=colorbar

for ic=1:nparts
    pos = [partition_lonlat(ic,1)-2.5 partition_lonlat(ic,2)-2.5 5 5];
     acolor = (q0T0(ic,1) )/0.14;
     
    if(q0T0(ic,1)>0 & q0T0(ic,1)<=0.02)
        color = cmp(1,:);
    end
    if(q0T0(ic,1)>0.02 & q0T0(ic,1)<=0.04)
        color = cmp(2,:);
    end
    if(q0T0(ic,1)>0.04 & q0T0(ic,1)<=0.06)
        color = cmp(3,:);
    end
    if(q0T0(ic,1)>0.06 & q0T0(ic,1)<=0.08)
        color = cmp(4,:);
    end
    if(q0T0(ic,1)>0.08 & q0T0(ic,1)<=0.1)
        color = cmp(5,:);
    end
    if(q0T0(ic,1)>0.1 & q0T0(ic,1)<=0.12)
        color = cmp(6,:);
    end
    if(q0T0(ic,1)>0.12 )
        color = cmp(7,:);
    end
    
    rectangle('position',pos,'FaceColor',color,'EdgeColor','none')
    
end

set(gca, 'CLim', [0, 0.14]);
set(cbh2,'YTick',[0,0.02,0.04,0.06,0.08,0.1,0.12,],'Location','eastoutside')

set(gca,'color',[0.8 0.8 0.8])
xlim([-150 160])
ylim([-50 77.5])
box on
pbaspect([2 1 1])

title('Heat Flow (Wm^-^2)')


ax0.TitleFontSizeMultiplier = 0.9;
ax1.TitleFontSizeMultiplier = 0.9;



x0=100;
y0=100;
width=900;
height=500

set(gcf,'position',[x0,y0,width,height])
fig = gcf;
fig.InvertHardcopy = 'off';

print -painters -depsc2 -r2500 plots/q0T0_all.eps
