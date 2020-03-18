close
clear


close
clear
load world.dat
% Code to read the samples from the RJ-MCMC algorithm.
% These are text/ascii files with s x time points, with temperatures
% W. These are interpolated here (as is done in the RJ-MCMC algorithm).
% Here the latest measurement date in each gridcell is read in and used
% to truncate the inferred temperature history.
% The averaging to regions is done after subtracting the CE 1955-180 
% mean value in each gridcell.
% Peter Hopcroft
% January 2019

% Add in correction for SAT-GST offset from Melo_Auilar et al 2018, their 
% figure 9.


sat_gst_corr=ncread('~/Desktop/Fig_9_data/all_mean_delT.nc','delT')
size(sat_gst_corr)

imagesc(sat_gst_corr)

load colocbh_no.txt
load BH_locats.txt
load bhlocats.txt

% the number of gridcells that we are reading in:
nparts=173;

% The number of MCMC iterations to sample from:
itmax=500;
    
    
partition_lonlat=colocbh_no;
partition_lonlat(:,1)=(colocbh_no(:,1))*5 - 2.5 ;
partition_lonlat(:,2)=(18 - colocbh_no(:,2)+1)*5- 2.5;
bhlocats(:,1) = bhlocats(:,1)*100*5 - 2.5 ;
bhlocats(:,2) = 90 - bhlocats(:,2)*37*5 +2.5;
% Load and plot the GST_sat correction:



corr_warming_6090_m_1800_1850=zeros(nparts);
corr=zeros(nparts);
for ic=1:nparts
 
    lathere=partition_lonlat(ic,2);
    lonhere=partition_lonlat(ic,1);
    ilat=round((lathere+90)/5 ,0)
    
    if(lonhere<0) lonhere=lonhere+360;
    end
        
    ilon=round(lonhere/5 ,0)
    
    
   
    corr(ic) = sat_gst_corr(ilon,ilat);
end


figure
plot(partition_lonlat(1:nparts,2),corr(1:nparts),'linestyle','none','marker','o')

% -----map corrected warming ----------------
figure 

hold all
plot(world(:,1),world(:,2),'color','black')

colormap(jet(7))
colorbar

cmp=colormap(jet(7))

for ic=1:nparts
    pos = [partition_lonlat(ic,1)-2.5 partition_lonlat(ic,2)-2.5 5 5];
    % acolor = (warming_6090_m_1800_1850(ic))/3.5;
     
    if(corr(ic)<=-0.2)
        color = cmp(1,:);
    end
    if(corr(ic)>-0.2 & corr(ic)<=0.0)
        color = cmp(2,:);
    end
    if(corr(ic)>0.0 & corr(ic)<=0.2)
        color = cmp(3,:);
    end
    if(corr(ic)>0.2 & corr(ic)<=0.4)
        color = cmp(4,:);
    end
    if(corr(ic)>0.4 & corr(ic)<=0.6)
        color = cmp(5,:);
    end
    if(corr(ic)>0.6 & corr(ic)<=0.8)
        color = cmp(6,:);
    end
    if(corr(ic)>0.8 )
        color = cmp(7,:);
    end
    if isnan(corr(ic))
        color = [0.7 0.7 0.7 ];
    end
    rectangle('position',pos,'FaceColor',color)
    
end

set(gca, 'CLim',[-0.4, 1.0] );
cbh=colorbar(gca);
set(cbh,'XTick',[-0.6:0.4:2.2],'XTickLabel',{'-0.6','-0.2','0.2','0.6','1.0','1.4','1.8',''});

xlim([-150 160])
ylim([-50 77.5])

title('mean temperature change (\circC)')
box on
pbaspect([2 1 1])




print -painters -depsc2 -r2500 glob_corr.eps






% ====================================================
% Corrected timeseries: (corrected for SAT_GST offset from Melar-Aguilo et
% al

yearsAD=[1405:5:2015];
glob=zeros(123,1);
SH=zeros(123,1);
temp=zeros(123,1);
% Now I need the latest logging date in each gridcell:

borehole_date = zeros(nparts,2);
borehole_z0 = zeros(nparts,2);
load minsofar.dat
borehole_date(:,2) = minsofar(1:nparts);

load minsofar2.dat;
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
   end_index
end


%corr(:)=0.06;

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


y1800 = (1800 - 1400)/5;
figure
hold on
for ic=1:nparts
    
    mean_G(ic,end_index(ic)+1:123) = NaN;
    clim6090=mean( mean_G(ic,i6090clim:f6090clim));
    if(isnan(clim6090))
        %
    else
    
    if(partition_lonlat(ic,2)>30 & partition_lonlat(ic,2)<60)  % Only extra-tropical sites
 area = (A * A*((cos(latrad(ic) - grid_lat_rad) - cos(latrad(ic)+grid_lat_rad))*2*pi*5)/360);
 
 
  
 for it=1:end_index(ic)
     if(it>y1800)
      glob(it) = glob(it) +  area *  ( (corr(ic)*(it-y1800)/(200/5)));
      temp(it) =  (corr(ic)*(it-y1800)/(200/5))
   else
         glob(it) = glob(it) + area *  (0.0);
         temp(it) =0.0;
     end
     
  
  count(it) = count(it)+area;
 end
 
 plot(temp)
 corr(ic)
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
   if(it>y1800)
       SH(it) = SH(it) + area *  ((corr(ic)*(it-y1800)/(200/5)) );
   else
    SH(it) = SH(it) + area *  (0.0);
       
   end
  countSH(it) = countSH(it)+area;
 end
        end
end
end

for it=1:123
 glob(it) = glob(it) /count(it);
 SH(it) = SH(it) /countSH(it);
end
 figure
subplot(2,1,1)
hold all
plot(yearsAD,glob)
xlim([1500 2020]); 
subplot(2,1,2)
hold all
plot(yearsAD,SH)
xlim([1500 2020]); 



save Tglob_new_corr_test.txt glob -ascii
save TSH_new_corr_test.txt SH -ascii
