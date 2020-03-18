
% ----------------------------------------------------------
% Now get the likelihood hyper-parameter
% ----------------------------------------------------------

clear
 
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
extension='NA.txt';
root2='X';

q0T0 = zeros(nparts,2);   % 104 partitions x 120 timesteps (5yrs)

for ic=1:nparts

    s= num2str(ic-1);
    filenameH=strcat(root,s,extension);
    varname = strcat(root2,s,'NA')
    %q0T0 = textread(filenameH, '%f','delimiter', '\n');
    load(filenameH);
    varhere=eval(varname);
    sq0 = size(varhere(1,:));
    sq0 = sq0(2);
    for is=1:sq0
      q0T0(ic,1) =q0T0(ic,1)+ (1.00/sq0)*mean(varhere(:,1));
      %q0T0(ic,2) =q0T0(ic,2)+ (1.00/sq0)*mean(varhere(:,2+(is-1)*2));
    end

 
end

ax0=figure
hold all
load inputs/world.dat
plot(world(:,1),world(:,2),'color','black')

xlim([-180 180])
ylim([-50 77.5])


hold all
cmp=colormap(jet(7))
cbh=colorbar

for ic=1:nparts
    pos = [partition_lonlat(ic,1)-2.5 partition_lonlat(ic,2)-2.5 5 5];
     acolor = (q0T0(ic,1) )/0.35;
     
    if(q0T0(ic,1)>0 & q0T0(ic,1)<=0.05)
        color = cmp(1,:);
    end
    if(q0T0(ic,1)>0.05 & q0T0(ic,1)<=0.1)
        color = cmp(2,:);
    end
    if(q0T0(ic,1)>0.1 & q0T0(ic,1)<=0.15)
        color = cmp(3,:);
    end
    if(q0T0(ic,1)>0.15 & q0T0(ic,1)<=0.2)
        color = cmp(4,:);
    end
    if(q0T0(ic,1)>0.2 & q0T0(ic,1)<=0.25)
        color = cmp(5,:);
    end
    if(q0T0(ic,1)>0.25 & q0T0(ic,1)<=0.3)
        color = cmp(6,:);
    end
    if(q0T0(ic,1)>0.3 )
        color = cmp(7,:);
    end
    
    rectangle('position',pos,'FaceColor',color,'EdgeColor','none')
    
end

save outputs/NA.txt q0T0 -ascii

set(gca, 'CLim', [0,0.35]);

set(gca,'color',[0.8 0.8 0.8])
set(cbh,'YTick',[0.0:0.05:0.3])
xlim([-150 160])
ylim([-50 77.5])
box on
pbaspect([2 1 1])

title('Likelihood hyper-parameter (\circC)')
ax0.TitleFontSizeMultiplier = 0.9;



x0=100;
y0=100;
width=900;
height=500

set(gcf,'position',[x0,y0,width,height])
fig = gcf;
fig.InvertHardcopy = 'off';

print -painters -depsc2 -r2500 plots/NA_all.eps
