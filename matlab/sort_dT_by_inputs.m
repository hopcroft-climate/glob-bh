close all
clear

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
itmax=10;

load inputs/minsofar.dat
borehole_date(1:nparts,2) = minsofar;

load inputs/minsofar2.dat
borehole_z0(1:nparts,2)= minsofar2;

load inputs/zmax.dat
borehole_zmax(1:nparts,2)= zmax;


    
    
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




% save these histories to a text file

load outputs/mean_G.txt 
load outputs/stdev_G.txt
load outputs/NA.txt
load outputs/pp_array.txt
load outputs/q0T0.txt



yearsAD=[1405:5:2015];
glob=zeros(123,1);
glob_p2s=zeros(123,1);
glob_m2s=zeros(123,1);
glob_p3s=zeros(123,1);
glob_m3s=zeros(123,1);
SH=zeros(123,1);
SH_p2s=zeros(123,1);
SH_m2s=zeros(123,1);
SH_p3s=zeros(123,1);
SH_m3s=zeros(123,1);
% Now I need the latest logging date in each gridcell:

borehole_date = zeros(nparts,2);
borehole_z0 = zeros(nparts,2);
load inputs/minsofar.dat
borehole_date(:,2) = minsofar(1:nparts);

load inputs/minsofar2.dat;
borehole_z0(:,2) = minsofar2(1:nparts);

% Subtract the depth from the date:
% Characteristic time period tau will propate
% sqrt(tau*kappa), where kappa diffusivity
% or l^2 /kappa for depth l

end_index =zeros(nparts,1);
for ic=1:nparts
  % depth_time = 0.3*((borehole_z0(ic,2))^2/1e-6)/(3600*24*365.25)
  % borehole_date(ic,2) = borehole_date(ic,2)  -depth_time;
   end_index(ic) =   int16(( 0.5+  borehole_date(ic,2) -2015  + 615) /5) 
   % Here deselect partitions with overly narrow prior limits
   %if (pp_array(ic) < 0.25)
   %    end_index(ic)=0;
   %end
end

count=zeros(123,1);
countSH=zeros(123,1);
A=6378185.0;
pi=3.141592654;
grid_lat_rad=2.5;
grid_lon_rad=2.5;
latrad=2*pi*partition_lonlat(:,2)/360.0;
lonrad=2*pi*partition_lonlat(:,1)/360.0;

f6090clim=(1980-1400)/5;
i6090clim=(1955-1400)/5;

for ic=1:nparts
    
    mean_G(ic,end_index(ic)+1:123) = NaN;
    clim6090=mean( mean_G(ic,i6090clim:f6090clim));
    if(isnan(clim6090))
        %
    else
    
    if(partition_lonlat(ic,2)>30 & partition_lonlat(ic,2)<60)  % Only extra-tropical sites
 area = (A * A*((cos(latrad(ic) - grid_lat_rad) - cos(latrad(ic)+grid_lat_rad))*2*pi*5)/360);
 
  
  
 for it=1:end_index(ic)
     
  glob(it) = glob(it) + area *  (mean_G(ic,it) - clim6090);
  glob_p2s(it) = glob_p2s(it) + area * (mean_G(ic,it)- clim6090 + 2*stdev_G(ic,it)) ;
  glob_m2s(it) = glob_m2s(it) + area * (mean_G(ic,it)- clim6090 - 2*stdev_G(ic,it)) ;
  glob_p3s(it) = glob_p3s(it) + area * (mean_G(ic,it)- clim6090 + 3*stdev_G(ic,it)) ;
  glob_m3s(it) = glob_m3s(it) + area * (mean_G(ic,it)- clim6090 - 3*stdev_G(ic,it)) ;
  count(it) = count(it)+area;
 end
    end
    end % if isnan(clim6090)
end

for ic=1:nparts
 mean_G(ic,end_index(ic)+1:123) = NaN;
 clim6090=mean( mean_G(ic,i6090clim:f6090clim));
if(isnan(clim6090))
        %
    else if(partition_lonlat(ic,2)<0)  % All SH sites
 area = (-A * A*((cos(latrad(ic) - grid_lat_rad) - cos(latrad(ic)+grid_lat_rad))*2*pi*5)/360);
 
 for it=1:end_index(ic)
  SH(it) = SH(it) + area *  (mean_G(ic,it) - clim6090);
  SH_p2s(it) = SH_p2s(it) + area * (mean_G(ic,it)- clim6090 + 2*stdev_G(ic,it)) ;
  SH_m2s(it) = SH_m2s(it) + area * (mean_G(ic,it)- clim6090 - 2*stdev_G(ic,it)) ;
  SH_p3s(it) = SH_p3s(it) + area * (mean_G(ic,it)- clim6090 + 3*stdev_G(ic,it)) ;
  SH_m3s(it) = SH_m3s(it) + area * (mean_G(ic,it)- clim6090 - 3*stdev_G(ic,it)) ;
  countSH(it) = countSH(it)+area;
 end
        end
end
end

for it=1:123
 glob(it) = glob(it) /count(it);
 glob_p2s(it) = glob_p2s(it) /count(it);
 glob_m2s(it) = glob_m2s(it) /count(it);
 glob_p3s(it) = glob_p3s(it) /count(it);
 glob_m3s(it) = glob_m3s(it) /count(it);
 SH(it) = SH(it) /countSH(it);
 SH_p2s(it) = SH_p2s(it) /countSH(it);
 SH_m2s(it) = SH_m2s(it) /countSH(it);
 SH_p3s(it) = SH_p3s(it) /countSH(it);
 SH_m3s(it) = SH_m3s(it) /countSH(it);
end
subplot(2,1,1)
hold all
plot(yearsAD,glob)
plot(yearsAD,glob_p2s)
plot(yearsAD,glob_m2s) 
plot(yearsAD,glob_p3s)
plot(yearsAD,glob_m3s)
xlim([1500 2020]); 
subplot(2,1,2)
hold all
plot(yearsAD,SH)
plot(yearsAD,SH_p2s)
plot(yearsAD,SH_m2s) 
plot(yearsAD,SH_p3s)
plot(yearsAD,SH_m3s)
xlim([1500 2020]); 



% --------------------------------------
% Now map plot the magnitude of warming
% --------------------------------------
iclim6090=i6090clim;
fclim6090=f6090clim;

% 1755-1800:
iclimLIA=71;
fclimLIA=iclimLIA+10;

warming_6090_m_1800_1850=zeros(nparts,1);
for ic=1:nparts
     if isnan(mean(mean_G(ic,iclim6090:fclim6090)))
        warming_6090_m_1800_1850(ic)=NaN;  
     else
        warming_6090_m_1800_1850(ic)= mean(mean_G(ic,iclim6090:fclim6090)) - mean(mean_G(ic,iclimLIA:fclimLIA));
     end
end

sd_warming_6090_m_1800_1850=zeros(nparts,1);
for ic=1:nparts
    sd_warming_6090_m_1800_1850(ic)= sqrt((mean(stdev_G(ic,iclim6090:fclim6090))/2)^2 + (mean(stdev_G(ic,iclimLIA:fclimLIA))/2)^2);
end

% Make it the 3 s.d.
sd_warming_6090_m_1800_1850=sd_warming_6090_m_1800_1850*3;

load inputs/world.dat

figure 
subplot(2,1,1)
hold all
plot(world(:,1),world(:,2),'color','black')

colormap(jet(7))
colorbar

cmp=colormap(jet(7))

for ic=1:nparts
    pos = [partition_lonlat(ic,1)-2.5 partition_lonlat(ic,2)-2.5 5 5];
    % acolor = (warming_6090_m_1800_1850(ic))/3.5;
     
    if(warming_6090_m_1800_1850(ic)>-0.6 & warming_6090_m_1800_1850(ic)<=-0.2)
        color = cmp(1,:);
    end
    if(warming_6090_m_1800_1850(ic)>-0.2 & warming_6090_m_1800_1850(ic)<=0.2)
        color = cmp(2,:);
    end
    if(warming_6090_m_1800_1850(ic)>0.2 & warming_6090_m_1800_1850(ic)<=0.6)
        color = cmp(3,:);
    end
    if(warming_6090_m_1800_1850(ic)>0.6 & warming_6090_m_1800_1850(ic)<=1.0)
        color = cmp(4,:);
    end
    if(warming_6090_m_1800_1850(ic)>1.0 & warming_6090_m_1800_1850(ic)<=1.4)
        color = cmp(5,:);
    end
    if(warming_6090_m_1800_1850(ic)>1.4 & warming_6090_m_1800_1850(ic)<=1.8)
        color = cmp(6,:);
    end
    if(warming_6090_m_1800_1850(ic)>1.8 )
        color = cmp(7,:);
    end
    if isnan(warming_6090_m_1800_1850(ic))
        color = [0.7 0.7 0.7 ];
    end
    rectangle('position',pos,'FaceColor',color)
    
end

set(gca, 'CLim',[-0.6, 2.2] );
cbh=colorbar(gca);
set(cbh,'XTick',[-0.6:0.4:2.2],'XTickLabel',{'-0.6','-0.2','0.2','0.6','1.0','1.4','1.8',''});

xlim([-150 160])
ylim([-50 77.5])

title('mean temperature change (\circC)')
box on
pbaspect([1.5 1 1])


subplot(2,1,2)
hold all
colormap(jet(7))
colorbar
plot(world(:,1),world(:,2),'color','black')

cmp=colormap(jet(7))

for ic=1:nparts
    pos = [partition_lonlat(ic,1)-2.5 partition_lonlat(ic,2)-2.5 5 5];
    
    if(sd_warming_6090_m_1800_1850(ic)>0 & sd_warming_6090_m_1800_1850(ic)<=0.05)
        color = cmp(1,:);
    end
    if(sd_warming_6090_m_1800_1850(ic)>0.05 & sd_warming_6090_m_1800_1850(ic)<=0.1)
        color = cmp(2,:);
    end
    if(sd_warming_6090_m_1800_1850(ic)>0.1 & sd_warming_6090_m_1800_1850(ic)<=0.15)
        color = cmp(3,:);
    end
    if(sd_warming_6090_m_1800_1850(ic)>0.15 & sd_warming_6090_m_1800_1850(ic)<=0.2)
        color = cmp(4,:);
    end
        if(sd_warming_6090_m_1800_1850(ic)>0.2 & sd_warming_6090_m_1800_1850(ic)<=0.25)
        color = cmp(5,:);
        end
            if(sd_warming_6090_m_1800_1850(ic)>0.25 & sd_warming_6090_m_1800_1850(ic)<=0.3)
        color = cmp(6,:);
            end
    
        if(sd_warming_6090_m_1800_1850(ic)>0.3 )
        color = cmp(7,:);
        end
   if isnan(warming_6090_m_1800_1850(ic))
        color = [0.7 0.7 0.7 ];
    end
    
        
    rectangle('position',pos,'FaceColor',color)
    
end

set(gca, 'CLim', [0, 0.35]);

xlim([-150 160])
ylim([-50 77.5])

title('3 s.d. temperature change (\circC)')
box on
pbaspect([1.5 1 1])
print -painters -depsc2 -r2500 plots/dT_glob.eps


grid_dT=zeros(73,37);
grid_dT_sd=zeros(73,37);
for ic=1:nparts
    pos = [partition_lonlat(ic,1)-2.5 partition_lonlat(ic,2)-2.5 5 5];
    i = colocbh_no(ic,1);
    j = colocbh_no(ic,2);
    
    grid_dT(i,j) = warming_6090_m_1800_1850(ic);
    grid_dT_sd(i,j) = sd_warming_6090_m_1800_1850(ic);
end




figure
ylim_val=0.6;
nxplots=3;
nyplots=3;
iplot=1;
hold on

nbins=40;

subplot(nxplots,nyplots,iplot)
hold on
g = [-3:0.25:3];

%[y1,x1]=hist(warming_6090_m_1800_1850,nbins)
%[h g ]=hist(warming_6090_m_1800_1850,nbins)
%plot(g,h,'color','red','linewidth',1.5);
%[h g ] = hist(warming_6090_m_1800_1850,nbins);
hold all;
%plot(g,h,'color',[0.6 0.6 0.6],'linewidth',1.5);
%bar(g,h,'facecolor',[0.6 0.6 0.6]);
pd0=fitdist(warming_6090_m_1800_1850,'Normal');
y=pdf(pd0,g);
plot(g,y,'Linewidth',2,'color',[0.6 0.6 0.6])
xlim([-2,3])
box on
a=nanmean(warming_6090_m_1800_1850)
hold on
line([a, a],[0, 30],'color',[0.5 0.5 0.5])
ylim([0 ylim_val])
pbaspect([1.5 1 1])
title('All gridcells')

% ------------------------------------------------
% Now subselect by NA
iplot=iplot+1;
mean_NA=mean(NA(:,1))
range=max(NA(:,1)) -min(NA(:,1))
lowq=0.33*range+min(NA(:,1));

warming_6090_m_1800_1850_lowNA = zeros(nparts,1);
warming_6090_m_1800_1850_lowNA(:) = NaN;
for ic=1:nparts

    if(NA(ic,1)<lowq)
       warming_6090_m_1800_1850_lowNA(ic) = warming_6090_m_1800_1850(ic);
    end
end

subplot(nxplots,nyplots,iplot)
%[y2,x2] = hist(warming_6090_m_1800_1850_lowNA,nbins)
%[h g ] = hist(warming_6090_m_1800_1850_lowNA,nbins)
%plot(g,h,'color','red','linewidth',1.5);
pd=fitdist(warming_6090_m_1800_1850_lowNA,'Normal');
y=pdf(pd,g);
plot(g,y,'Linewidth',2,'color','red')

%bar(g,h,'facecolor','red');
%[h g ] = hist(warming_6090_m_1800_1850,nbins);
hold all;
%plot(g,h,'color',[0.6 0.6 0.6],'linewidth',1.5);xlim([-2,3])
%bar(g,h,'facecolor',[0.6 0.6 0.6]);
pd=fitdist(warming_6090_m_1800_1850,'Normal');
y=pdf(pd,g);
plot(g,y,'Linewidth',2,'color',[0.6 0.6 0.6])


title('Low noise profiles')
a2=nanmean(warming_6090_m_1800_1850_lowNA)
hold on
line([a2, a2],[0, 30],'color','red')
line([a, a],[0, 30],'color',[0.5 0.5 0.5])
ylim([0 ylim_val])
pbaspect([1.5 1 1])
pd=fitdist(warming_6090_m_1800_1850,'Normal');
y=pdf(pd,g);
plot(g,y,'Linewidth',2,'color',[0.6 0.6 0.6])
box on
xlim([-2,3])
% ------------------------------------------------
% No subselect by PP
iplot=iplot+1;
mean_PP=mean(pp_array)
range=max(pp_array(:,1)) -min(pp_array(:,1))
lowq=0.33*range+min(pp_array(:,1));

warming_6090_m_1800_1850_lowPP = zeros(nparts,1);
warming_6090_m_1800_1850_lowPP(:) = NaN;
for ic=1:nparts
    if(pp_array(ic,1)<lowq)
       warming_6090_m_1800_1850_lowPP(ic) = warming_6090_m_1800_1850(ic);
    end
end
subplot(nxplots,nyplots,iplot)
%[y3,x3]=hist(warming_6090_m_1800_1850_lowPP,nbins);
%[h g ] = hist(warming_6090_m_1800_1850_lowPP,nbins);
%plot(g,h,'color','red','linewidth',1.5);
%[h g ] = hist(warming_6090_m_1800_1850,nbins);
hold all;
%plot(g,h,'color',[0.6 0.6 0.6],'linewidth',1.5);xlim([-2,3]);
pd=fitdist(warming_6090_m_1800_1850_lowPP,'Normal');
y=pdf(pd,g);
plot(g,y,'Linewidth',2,'color','red')


title('Low prior on climate variability')
a3=nanmean(warming_6090_m_1800_1850_lowPP)
hold on
line([a3, a3],[0, 30],'color','red')
line([a, a],[0, 30],'color',[0.5 0.5 0.5])
ylim([0 ylim_val])
pbaspect([1.5 1 1])
pd=fitdist(warming_6090_m_1800_1850,'Normal');
y=pdf(pd,g);
plot(g,y,'Linewidth',2,'color',[0.6 0.6 0.6])
box on
xlim([-2,3])
% ------------------------------------------------
% No subselect by q0
iplot=iplot+1;
mean_q0=mean(q0T0(:,1));
range=max(q0T0(:,1)) -min(q0T0(:,1))
lowq=0.33*range+min(q0T0(:,1));

warming_6090_m_1800_1850_lowq0 = zeros(nparts,1);
warming_6090_m_1800_1850_lowq0(:) = NaN;
for ic=1:nparts
    if(q0T0(ic,1)<lowq(1,1))
       warming_6090_m_1800_1850_lowq0(ic) = warming_6090_m_1800_1850(ic);
    end
end
subplot(nxplots,nyplots,iplot)
%[y4,x4]=hist(warming_6090_m_1800_1850_lowq0,nbins);
%[h g]=hist(warming_6090_m_1800_1850_lowq0,nbins);
%plot(g,h,'color','red','linewidth',1.5);
%[h g ] = hist(warming_6090_m_1800_1850,nbins);
hold all;
%plot(g,h,'color',[0.6 0.6 0.6],'linewidth',1.5);xlim([-2,3]);
pd=fitdist(warming_6090_m_1800_1850_lowq0,'Normal');
y=pdf(pd,g);
plot(g,y,'Linewidth',2,'color','red')

title('Low basal heat flux')
a4=nanmean(warming_6090_m_1800_1850_lowq0)
hold on
line([a4, a4],[0, 30],'color','red')
line([a, a],[0, 30],'color',[0.5 0.5 0.5])
ylim([0 ylim_val])
pbaspect([1.5 1 1])
ylabel('probability density')
pd=fitdist(warming_6090_m_1800_1850,'Normal');
y=pdf(pd,g);
plot(g,y,'Linewidth',2,'color',[0.6 0.6 0.6])
box on
xlim([-2,3])
% ------------------------------------------------
% No subselect by low T0
iplot=iplot+1;
mean_T0=15.0;

warming_6090_m_1800_1850_lowT0 = zeros(nparts,1);
warming_6090_m_1800_1850_lowT0(:) = NaN;
for ic=1:nparts
    if(q0T0(ic,2)<mean_T0(1,1))
       warming_6090_m_1800_1850_lowT0(ic) = warming_6090_m_1800_1850(ic);
    end
end
subplot(nxplots,nyplots,iplot)
%[y5,x5]=hist(warming_6090_m_1800_1850_lowT0,nbins);
%[h g ]=hist(warming_6090_m_1800_1850_lowT0,nbins);
%plot(g,h,'color','red','linewidth',1.5);
xlim([-2,3]);
%[h g ] = hist(warming_6090_m_1800_1850,nbins);
hold all;
%plot(g,h,'color',[0.6 0.6 0.6],'linewidth',1.5);title('Low site temperature')
pd=fitdist(warming_6090_m_1800_1850_lowT0,'Normal');
y=pdf(pd,g);
plot(g,y,'Linewidth',2,'color','red')

a5=nanmean(warming_6090_m_1800_1850_lowT0)
hold on
line([a5, a5],[0, 30],'color','red')
line([a, a],[0, 30],'color',[0.5 0.5 0.5])
ylim([0 ylim_val])
pbaspect([1.5 1 1])
pd=fitdist(warming_6090_m_1800_1850,'Normal');
y=pdf(pd,g);
plot(g,y,'Linewidth',2,'color',[0.6 0.6 0.6])
title('Low mean temperature')
box on
xlim([-2,3])
% ------------------------------------------------
% No subselect by high T0
iplot=iplot+1;
mean_T0=15.0;

warming_6090_m_1800_1850_lowT0 = zeros(nparts,1);
warming_6090_m_1800_1850_lowT0(:) = NaN;
for ic=1:nparts
    if(q0T0(ic,2)>mean_T0(1,1))
       warming_6090_m_1800_1850_lowT0(ic) = warming_6090_m_1800_1850(ic);
    end
end
subplot(nxplots,nyplots,iplot)
%[y6,x6]=hist(warming_6090_m_1800_1850_lowT0,nbins);
%[h g]=hist(warming_6090_m_1800_1850_lowT0,nbins);
%plot(g,h,'color','red','linewidth',1.5);xlim([-2,3]);
%[h g ] = hist(warming_6090_m_1800_1850,nbins);
hold all;
%plot(g,h,'color',[0.6 0.6 0.6],'linewidth',1.5);title('High site temperature')
pd=fitdist(warming_6090_m_1800_1850_lowT0,'Normal');
y=pdf(pd,g);
plot(g,y,'Linewidth',2,'color','red')


title('High mean temperature')

a6=nanmean(warming_6090_m_1800_1850_lowT0)
hold on
line([a6, a6],[0, 30],'color','red')
line([a, a],[0, 30],'color',[0.5 0.5 0.5])
ylim([0 ylim_val])
pbaspect([1.5 1 1])

pd=fitdist(warming_6090_m_1800_1850,'Normal');
y=pdf(pd,g);
plot(g,y,'Linewidth',2,'color',[0.6 0.6 0.6])


box on
xlim([-2,3])
% box on------------------------------------------------
% Now subselect by high borehole_z0
iplot=iplot+1;
mean_z0=25;

warming_6090_m_1800_1850_lowz0 = zeros(nparts,1);
warming_6090_m_1800_1850_lowz0(:) = NaN;
for ic=1:nparts
    if(borehole_z0(ic,2)<mean_z0(1,1))
       warming_6090_m_1800_1850_lowz0(ic) = warming_6090_m_1800_1850(ic);
    end
end
subplot(nxplots,nyplots,iplot)
%[y7,x7]=hist(warming_6090_m_1800_1850_lowz0,nbins);
%[h g ] = hist(warming_6090_m_1800_1850_lowz0,nbins);
%plot(g,h,'color','red','linewidth',1.5);
xlim([-2,3]);
pd=fitdist(warming_6090_m_1800_1850_lowz0,'Normal');
y=pdf(pd,g);
plot(g,y,'Linewidth',2,'color','red')

%[h g ] = hist(warming_6090_m_1800_1850,nbins);
hold all;
%plot(g,h,'color',[0.6 0.6 0.6],'linewidth',1.5);
pd=fitdist(warming_6090_m_1800_1850,'Normal');
y=pdf(pd,g);
plot(g,y,'Linewidth',2,'color',[0.6 0.6 0.6])

title('Shallow top <25m')



a7=nanmean(warming_6090_m_1800_1850_lowz0);
hold on
line([a7, a7],[0, 30],'color','red')
line([a, a],[0, 30],'color',[0.5 0.5 0.5])
ylim([0 ylim_val])
pbaspect([1.5 1 1])
box on
xlim([-2,3])
% ------------------------------------------------
% No subselect by high borehole_zmax
iplot=iplot+1;
mean_zmax=450.0

warming_6090_m_1800_1850_lowzmax = zeros(nparts,1);
warming_6090_m_1800_1850_lowzmax(:) = NaN;
for ic=1:nparts
    if(borehole_zmax(ic,2)>mean_zmax(1,1))
       warming_6090_m_1800_1850_lowzmax(ic) = warming_6090_m_1800_1850(ic);
    end
end
subplot(nxplots,nyplots,iplot)
%[y8, x8] =hist(warming_6090_m_1800_1850_lowzmax,nbins);
%[h g]=hist(warming_6090_m_1800_1850_lowzmax,nbins);
%plot(g,h,'color','red','linewidth',1.5);
%[h g ] = hist(warming_6090_m_1800_1850,nbins);
hold all;
pd=fitdist(warming_6090_m_1800_1850_lowzmax,'Normal');
y=pdf(pd,g);
plot(g,y,'Linewidth',2,'color','red')
%plot(g,h,'color',[0.6 0.6 0.6],'linewidth',1.5);
%pd=fitdist(warming_6090_m_1800_1850,'Normal');
y=pdf(pd0,g);
plot(g,y,'Linewidth',2,'color',[0.6 0.6 0.6])

xlim([-2,3]);
title('Deeper profiles > 450m')

a8=nanmean(warming_6090_m_1800_1850_lowzmax)
hold on
line([a8, a8],[0, 30],'color','red')
line([a, a],[0, 30],'color',[0.5 0.5 0.5])
ylim([0 ylim_val])
pbaspect([1.5 1 1])
xlabel('warming (K)')
box on
xlim([-2,3])
% ------------------------------------------------
% Now subselect by date
iplot=iplot+1;
mean_date=1985;

warming_6090_m_1800_1850_lowdate = zeros(nparts,1);
warming_6090_m_1800_1850_lowdate(:) = NaN;
for ic=1:nparts
    if(borehole_date(ic,2)>mean_date(1,1))
       warming_6090_m_1800_1850_lowdate(ic) = warming_6090_m_1800_1850(ic);
    end
end
subplot(nxplots,nyplots,iplot)
%[h g ] = hist(warming_6090_m_1800_1850_lowdate,nbins);
%plot(g,h,'color','red','linewidth',1.5);
pd=fitdist(warming_6090_m_1800_1850_lowdate,'Normal');
y=pdf(pd,g);
plot(g,y,'Linewidth',2,'color','red')
hold all
%[h g ] = hist(warming_6090_m_1800_1850,nbins);
hold all;
%plot(g,h,'color',[0.6 0.6 0.6],'linewidth',1.5);
pd=fitdist(warming_6090_m_1800_1850,'Normal');
y=pdf(pd,g);
plot(g,y,'Linewidth',2,'color',[0.6 0.6 0.6])


xlim([-2,3]);
title('Profiles measure after CE 1985')


a9=nanmean(warming_6090_m_1800_1850_lowdate)
hold on
line([a9, a9],[0, 30],'color','red')
line([a, a],[0, 30],'color',[0.5 0.5 0.5])
ylim([0 ylim_val])
pbaspect([1.5 1 1])

pbaspect([1.5 1 1])
print -painters -depsc2 -r2500 plots/sort_dT_hist_pars.eps

figure
subplot(3,3,1)
hist(NA(:,1),nbins)
title({'Estimated observational',' noise (K)'})
pbaspect([1.5 1 1])
subplot(3,3,2)
hist(pp_array,nbins)
title({'Estimated climate','variability prior (K)'})
pbaspect([1.5 1 1])
subplot(3,3,3)
hist(q0T0(:,1),nbins)
title({'Basal heat','flux (Wm-2)'})
pbaspect([1.5 1 1])
subplot(3,3,4)
hist(q0T0(:,2),nbins)
title({'Long-term mean surface','temperature (\circC)'})
pbaspect([1.5 1 1])
subplot(3,3,5)
hist(borehole_z0(:,2),nbins)
title({'Minimum depth of','profile top (m)'})
pbaspect([1.5 1 1])
subplot(3,3,6)
hist(borehole_zmax(:,2),nbins)
title({'Max depth of','profile base (m)'})
pbaspect([1.5 1 1])
subplot(3,3,7)
hist(borehole_date(:,2),nbins)
title({'Latest date of profile','measurment (yr CE)'})


pbaspect([1.5 1 1])
print -painters -depsc2 -r2500 plots/hist_pars.eps