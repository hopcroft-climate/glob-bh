close all
clear

%cd ~/DataShare/WORK/glob_bh

% Code to read the samples from the RJ-MCMC algorithm.
% These are text/ascii files with s x time points, with temperatures
% W. These are interpolated here (as is done in the RJ-MCMC algorithm).
% Here the latest measurement date in each gridcell is read in and used
% to truncate the inferred temperature history.
% The averaging to regions is done after subtracting the CE 1955-180 
% mean value in each gridcell.
% Peter Hopcroft
% January 2019

load inputs/colocbh_no.txt
load inputs/BH_locats.txt
load inputs/bhlocats.txt

% the number of gridcells that we are reading in:
nparts=173;

% The number of MCMC iterations to sample from:
itmax=100000;
    
    
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
glob_partition_lonlat(1:nparts,:)=partition_lonlat(1:nparts,:);

hold all
plot(bh_long,bh_lat,'linestyle','none','marker','s')
plot(BH_locats(:,1),BH_locats(:,2),'linestyle','none','marker','p')
load inputs/world.dat
plot(world(:,1),world(:,2),'color','black')

root='~/Documents/WORK/NH_SH_bh/glob_bh_new_hcoef/Samples/';
%root='C:\Users\hopcrofp\Documents\NH_SH_bh\glob_bh_new_lcoef\Samples\';

extensionpp='ALPHA.txt';
root2='X';


pp_array=zeros(nparts,6);
pp_array_std=zeros(nparts,6);

for ic=1:nparts
   s= num2str(ic-1);
   filenamepp=strcat(root,s,extensionpp);
   varname = strcat(root2,s,'ALPHA')
   
   %pp = textread(filenamepp, '%f','delimiter', '\n');
   load(filenamepp);
   varhere=eval(varname);
   for ic2=1:6
    pp_array(ic,ic2) = mean(varhere(50000:100000,ic2));
    pp_array_std(ic,ic2) = std(varhere(50000:100000,ic2));
   end
end


save outputs_lcoef/ALPHA_array.txt pp_array -ascii
save outputs_lcoef/ALPHA_array_std.txt pp_array_std -ascii


% ====================================================
% ====================================================

pp_diff=zeros(nparts,1);
for ic=1:nparts
     pp_diff(ic,1) = (pp_array(ic,1))  - mean(pp_array(ic,2:3));
end



% Map the prior hyperparameter pp
figure

subplot(2,1,1)
hold all
lev=zeros(7);
lev(1) = -3;
lev(2)= -2;
lev(3) = -1;
lev(4) = 1;
lev(5) = 2
lev(6) = 3;
lev(7)=4;



plot(world(:,1),world(:,2),'color','black')


cmp=colormap(parula(8))
colorbar
for ic=1:nparts
    pos = [partition_lonlat(ic,1)-2.5 partition_lonlat(ic,2)-2.5 5 5];
    if( pp_diff(ic)<=lev(1))
        color = cmp(1,:);
    end
    if(pp_diff(ic)>lev(1) & pp_diff(ic)<=lev(2))
        color = cmp(2,:);
    end
    if(pp_diff(ic)>lev(2) & pp_diff(ic)<=lev(3))
        color = cmp(3,:);
    end
    if(pp_diff(ic)>lev(3) & pp_diff(ic)<=lev(4))
        color =cmp(4,:);
    end
    if(pp_diff(ic)>lev(4) & pp_diff(ic)<=lev(5))
        color = cmp(5,:);
    end
        if(pp_diff(ic)>lev(5) & pp_diff(ic)<=lev(6))
        color = cmp(6,:);
        end
     if(pp_diff(ic)>lev(6) & pp_diff(ic)<=lev(7))
        color = cmp(7,:);
      end
        if(pp_diff(ic)>lev(7) )
        color = cmp(8,:);
        end
   % if(borehole_date(ic,2) < 1980)
   %    color = [0.7 0.7 0.7 ];
   % end
    rectangle('position',pos,'FaceColor',color,'EdgeColor','None')
    
end

set(gca, 'CLim', [-4, 5]);
set(gca,'color',[0.8 0.8 0.8])
xlim([-150 160])
ylim([-50 77.5])


title('coupling coefficient (Wm^-2/\circC)')
box on
colorbar
pbaspect([2 1 1])
x0=100;
y0=100;
width=900;
height=500
set(gcf,'position',[x0,y0,width,height])
fig = gcf;
fig.InvertHardcopy = 'off';
print -painters -depsc2 -r2500 plots/alpha_glob.eps




%subplot(2,2,1)
%imagesc(pp_array)
%colorbar


subplot(2,1,2)
count_nh=0.0;
count_sh=0.0
mean_nh=zeros(6,1);
mean_sh=zeros(6,1);

std_nh=zeros(6,1);
std_sh=zeros(6,1);


for ic=1:nparts
    if(glob_partition_lonlat(ic,2)>30 & glob_partition_lonlat(ic,2)<70)
        mean_nh=mean_nh+pp_array(ic,1:6)';
        std_nh=std_nh+pp_array_std(ic,1:6)';
        count_nh=count_nh+1.0;
    end
    if(glob_partition_lonlat(ic,2)<0)
        mean_sh=mean_sh+pp_array(ic,1:6)';
        std_sh=std_sh+pp_array_std(ic,1:6)';
        count_sh=count_sh+1.0;
    end
end
for ic2=1:6
   mean_nh(ic2) = mean_nh(ic2)/count_nh;
   mean_sh(ic2) = mean_sh(ic2)/count_sh;
   std_nh(ic2) = std_nh(ic2)/count_nh;
   std_sh(ic2) = std_sh(ic2)/count_sh;
end
%figure
centuries=[2015,1900,1800,1700,1600,1500];
hold all
stairs(centuries,mean_nh,'linewidth',2.0,'color','blue')
stairs(centuries,mean_nh+std_nh,'linewidth',1.0,'color','blue')
stairs(centuries,mean_nh-std_nh,'linewidth',1.0,'color','blue')
xlim([1500 2020])
%ylim([1 4])
hold all
stairs(centuries,mean_sh,'linewidth',2.0,'color','red')
stairs(centuries,mean_sh+std_sh,'linewidth',1.0,'color','red')
stairs(centuries,mean_sh-std_sh,'linewidth',1.0,'color','red')

legend('NH extra-tropics','SH','fontsize',16','location','northwest')
ylabel('coupling coefficient (Wm^-2/\circC)')
xlabel('year CE')
set(gcf,'position',[x0,y0,width,height])
fig = gcf;
fig.InvertHardcopy = 'off';
print -painters -depsc2 -r2500 plots/ALPHA_lcoef.eps

