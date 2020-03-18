
close
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

% Add in correction for SAT-GST offset from Melo_Auilar et al 2018, their 
% figure 9.


sat_gst_corr=ncread('Fig_9_data/all_mean_delT.nc','delT')  % This is the 2-century change so divide by 2:

sat_gst_corr=sat_gst_corr/2.;

size(sat_gst_corr)

sat_gst_corr=sat_gst_corr  %  because this is SAT-GST

% test
sat_gst_corr(:,:)=0.1;

imagesc(sat_gst_corr)

load colocbh_no.txt
load BH_locats.txt
load bhlocats.txt

% the number of gridcells that we are reading in:
nparts=173;

% The number of MCMC iterations to sample from:
itmax=50;
    
    
partition_lonlat=colocbh_no;
partition_lonlat(:,1)=(colocbh_no(:,1))*5 - 2.5 ;
partition_lonlat(:,2)=(18 - colocbh_no(:,2)+1)*5- 2.5;
bhlocats(:,1) = bhlocats(:,1)*100*5 - 2.5 ;
bhlocats(:,2) = 90 - bhlocats(:,2)*37*5 +2.5;




for ic=1:nparts
 
    lathere=partition_lonlat(ic,2);
    lonhere=partition_lonlat(ic,1);
    ilat=round((lathere+90)/5 ,0)
    
    if(lonhere<0) lonhere=lonhere+360;
    end
    ilon=round(lonhere/5 ,0)
    
    
    if(ilon<0)
        stop
    end
    if(ilat<0)
        stop
    end
    
    %warming_60_90_m_1800_1850(ic) = warming_60_90_m_1800_1850(ic) -  sat_gst_corr(ilat,ilon);
end
    





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
load world.dat
plot(world(:,1),world(:,2),'color','black')



root='Samples_ver0/';
extension='G.txt';

extensiont='T.txt';
extensions='S.txt';
extensionpp='PP.txt';
root2='X';
extensions2='G';
mean_G = zeros(nparts,123);   % nparts partitions x 123 timesteps (5yrs)
stdev_G = zeros(nparts,123);   % nparts partitions x 123 timesteps (5yrs)
times5yr = [5:5:615];
pp_array=zeros(nparts,1);
i=0;
for ic=1:nparts

    s= num2str(ic-1);
    filenameG=strcat(root,s,extension)
    filenamet=strcat(root,s,extensiont);
    filenames=strcat(root,s,extensions);
    filenamepp=strcat(root,s,extensionpp);
    
    
    GST = textread(filenameG, '%f','delimiter', '\n');
    t = textread(filenamet, '%f','delimiter', '\n');
    s = textread(filenames, '%f','delimiter', '\n');
    pp = textread(filenamepp, '%f','delimiter', '\n');
    % now interpolate to 5 year timesteps
    
    i=0;
    j=0;

    for it=1:itmax
        
        i = j+1;
        j = i + s(it) -1;
        tmcmc = flipud(t(i:j));
        wmcmc = flipud(GST(i:j));
       mean_G(ic,:) =mean_G(ic,:)+(1/itmax)* interp1(tmcmc,wmcmc,times5yr); 
    end
     i=0;
    j=0;
     for it=1:itmax
        i = j+1;
        j = i + s(it) -1;
        tmcmc = flipud(t(i:j));
        wmcmc = flipud(GST(i:j));
       dev2_G = (interp1(tmcmc,wmcmc,times5yr) - mean_G(ic,:)).^2;
       stdev_G(ic,:) = stdev_G(ic,:) + (1/itmax)* dev2_G;
       pp_array(ic) = mean(pp);
     end
end

stdev_G = sqrt(stdev_G);
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



save Tglob_new.txt glob -ascii
save Tglob_new_p2s.txt glob_p2s -ascii
save Tglob_new_m2s.txt glob_m2s -ascii

save Tglob_new_p3s.txt glob_p3s -ascii
save Tglob_new_m3s.txt glob_m3s -ascii

save TSH_new.txt SH -ascii
save TSH_new_p2s.txt SH_p2s -ascii
save TSH_new_m2s.txt SH_m2s -ascii

save TSH_new_p3s.txt SH_p3s -ascii
save TSH_new_m3s.txt SH_m3s -ascii



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


% Write out the  in the mean change:

%nccreate('mean_del_GST.nc','del_GST')
%nccreate('mean_del_GST.nc','sd_GST')
%ncwrite('mean_del_GST.nc','del_GST',warming_6090_m_1800_1850)
%ncwrite('mean_del_GST.nc','sd_GST',sd_warming_6090_m_1800_1850)


corr_warming_6090_m_1800_1850=zeros(nparts);
corr=zeros(nparts);
for ic=1:nparts
 
    lathere=partition_lonlat(ic,2);
    lonhere=partition_lonlat(ic,1);
    ilat=round((lathere+90)/5 ,0)
    
    if(lonhere<0) lonhere=lonhere+360;
    end
        
    ilon=round(lonhere/5 ,0)
    
    
    corr_warming_6090_m_1800_1850(ic) = warming_6090_m_1800_1850(ic) +  sat_gst_corr(ilon,ilat);
    
    if(lathere>40)
    corr(ic) = 0.25
    else
        corr(ic)=0.
    end
end


figure
plot(partition_lonlat(1:nparts,2),corr(1:nparts),'linestyle','none','marker','o')


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
pbaspect([2 1 1])


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
pbaspect([2 1 1])
print -painters -depsc2 -r2500 dT_glob.eps



% -----map corrected warming ----------------
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
     
    if(corr_warming_6090_m_1800_1850(ic)>-0.6 & corr_warming_6090_m_1800_1850(ic)<=-0.2)
        color = cmp(1,:);
    end
    if(corr_warming_6090_m_1800_1850(ic)>-0.2 & corr_warming_6090_m_1800_1850(ic)<=0.2)
        color = cmp(2,:);
    end
    if(corr_warming_6090_m_1800_1850(ic)>0.2 & corr_warming_6090_m_1800_1850(ic)<=0.6)
        color = cmp(3,:);
    end
    if(corr_warming_6090_m_1800_1850(ic)>0.6 & corr_warming_6090_m_1800_1850(ic)<=1.0)
        color = cmp(4,:);
    end
    if(corr_warming_6090_m_1800_1850(ic)>1.0 & corr_warming_6090_m_1800_1850(ic)<=1.4)
        color = cmp(5,:);
    end
    if(corr_warming_6090_m_1800_1850(ic)>1.4 & corr_warming_6090_m_1800_1850(ic)<=1.8)
        color = cmp(6,:);
    end
    if(corr_warming_6090_m_1800_1850(ic)>1.8 )
        color = cmp(7,:);
    end
    if isnan(corr_warming_6090_m_1800_1850(ic))
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
pbaspect([2 1 1])


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
pbaspect([2 1 1])
print -painters -depsc2 -r2500 dT_glob_corr.eps

% -----------------------------------------------------------------






% -----------------------------------------------------------------






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
subplot(2,1,1)
imagesc(grid_dT)

subplot(2,1,2)
imagesc(grid_dT_sd)


%nccreate('dT.nc','grid_dT')
%nccreate('dT_sd.nc','grid_dT_sd')
%ncwrite('dT.nc','grid_dT',grid_dT);
%ncwrite('dT_sd.nc','grid_dT_sd',grid_dT_sd);





% ====================================================
% Corrected timeseries: (corrected for SAT_GST offset from Melar-Aguilo et
% al
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
count(:)=0;

%corr(:)=0.0;

f6090clim=(1980-1400)/5;
i6090clim=(1955-1400)/5;

y1600 = 1+(1600 - 1400)/5;
y1700 = 1+(1700 - 1400)/5;
y1750 = 1+(1750 - 1400)/5;
y1700 = 1+(1700 - 1400)/5;
y1900 = 1+(1900 - 1400)/5;

ycorr= y1750;


for ic=1:nparts
    
    mean_G(ic,end_index(ic)+1:123) = NaN;
    clim6090=mean( mean_G(ic,i6090clim:f6090clim));
    if(isnan(clim6090))
        %
    else
    
    if(partition_lonlat(ic,2)>30 & partition_lonlat(ic,2)<60)  % Only extra-tropical sites
 area = (A * A*((cos(latrad(ic) - grid_lat_rad) - cos(latrad(ic)+grid_lat_rad))*2*pi*5)/360);
 
  
  
 for it=1:end_index(ic)
     if(it>ycorr)
  glob(it) = glob(it) + area *  (mean_G(ic,it) - clim6090 + (corr(ic)*(it-ycorr)/(100/5)));
  glob_p2s(it) = glob_p2s(it) + area * (mean_G(ic,it)- clim6090 + 2*stdev_G(ic,it)+ (corr(ic)*(it-ycorr)/(100/5))) ;
  glob_m2s(it) = glob_m2s(it) + area * (mean_G(ic,it)- clim6090 - 2*stdev_G(ic,it)+ (corr(ic)*(it-ycorr)/(100/5))) ;
  glob_p3s(it) = glob_p3s(it) + area * (mean_G(ic,it)- clim6090 + 3*stdev_G(ic,it)+ (corr(ic)*(it-ycorr)/(100/5))) ;
  glob_m3s(it) = glob_m3s(it) + area * (mean_G(ic,it)- clim6090 - 3*stdev_G(ic,it)+ (corr(ic)*(it-ycorr)/(100/5))) ;
  
 (it-ycorr)/(100/5)
  
     else
         glob(it) = glob(it) + area *  (mean_G(ic,it) - clim6090);
     glob_p2s(it) = glob_p2s(it) + area * (mean_G(ic,it)- clim6090 + 2*stdev_G(ic,it)) ;
  glob_m2s(it) = glob_m2s(it) + area * (mean_G(ic,it)- clim6090 - 2*stdev_G(ic,it)) ;
  glob_p3s(it) = glob_p3s(it) + area * (mean_G(ic,it)- clim6090 + 3*stdev_G(ic,it)) ;
  glob_m3s(it) = glob_m3s(it) + area * (mean_G(ic,it)- clim6090 - 3*stdev_G(ic,it)) ;
  
     end
     
  
  count(it) = count(it)+area;
 end
    end
    end % if isnan(clim6090)
end


countSH(:)=0;
for ic=1:nparts
 mean_G(ic,end_index(ic)+1:123) = NaN;
 clim6090=mean( mean_G(ic,i6090clim:f6090clim));
if(isnan(clim6090))
        %
    else if(partition_lonlat(ic,2)<0)  % All SH sites
 area = (-A * A*((cos(latrad(ic) - grid_lat_rad) - cos(latrad(ic)+grid_lat_rad))*2*pi*5)/360);
 
 for it=1:end_index(ic)
   if(it>ycorr)
      SH(it) = SH(it)     + area * (mean_G(ic,it) + clim6090 +  (corr(ic)*(it-ycorr)/(100/5)) );
  SH_p2s(it) = SH_p2s(it) + area * (mean_G(ic,it)- clim6090 + 2*stdev_G(ic,it)+ (corr(ic)*(it-ycorr)/(100/5))) ;
  SH_m2s(it) = SH_m2s(it) + area * (mean_G(ic,it)- clim6090 - 2*stdev_G(ic,it)+ (corr(ic)*(it-ycorr)/(100)/5)) ;
  SH_p3s(it) = SH_p3s(it) + area * (mean_G(ic,it)- clim6090 + 3*stdev_G(ic,it)+ (corr(ic)*(it-ycorr)/(100)/5)) ;
  SH_m3s(it) = SH_m3s(it) + area * (mean_G(ic,it)- clim6090 - 3*stdev_G(ic,it)+ (corr(ic)*(it-ycorr)/(100)/5)) ;
   else
    SH(it) = SH(it) + area *  (mean_G(ic,it) - clim6090);
  SH_p2s(it) = SH_p2s(it) + area * (mean_G(ic,it)- clim6090 + 2*stdev_G(ic,it)) ;
  SH_m2s(it) = SH_m2s(it) + area * (mean_G(ic,it)- clim6090 - 2*stdev_G(ic,it)) ;
  SH_p3s(it) = SH_p3s(it) + area * (mean_G(ic,it)- clim6090 + 3*stdev_G(ic,it)) ;
  SH_m3s(it) = SH_m3s(it) + area * (mean_G(ic,it)- clim6090 - 3*stdev_G(ic,it)) ;
       
   end
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



save outputs/Tglob_new_corr.txt glob -ascii
save outputs/Tglob_new_p2s_corr.txt glob_p2s -ascii
save outputs/Tglob_new_m2s_corr.txt glob_m2s -ascii

save outputs/Tglob_new_p3s_corr.txt glob_p3s -ascii
save outputs/Tglob_new_m3s_corr.txt glob_m3s -ascii

save outputs/TSH_new_corr.txt SH -ascii
save outputs/TSH_new_p2s_corr.txt SH_p2s -ascii
save outputs/TSH_new_m2s_corr.txt SH_m2s -ascii

save outputs/TSH_new_p3s_corr.txt SH_p3s -ascii
save outputs/TSH_new_m3s_corr.txt SH_m3s -ascii
% ====================================================

% Map the prior hyperparameter pp
figure
hold all
plot(world(:,1),world(:,2),'color','black')

cmp=colormap(jet(7))

for ic=1:nparts
    pos = [partition_lonlat(ic,1)-2.5 partition_lonlat(ic,2)-2.5 5 5];
    if(pp_array(ic)>0 & pp_array(ic)<=0.75)
        color = cmp(1,:);
    end
    if(pp_array(ic)>0.75 & pp_array(ic)<=1.0)
        color = cmp(2,:);
    end
    if(pp_array(ic)>1.0 & pp_array(ic)<=1.25)
        color = cmp(3,:);
    end
    if(pp_array(ic)>1.25 & pp_array(ic)<=1.5)
        color = cmp(4,:);
    end
        if(pp_array(ic)>1.5 & pp_array(ic)<=1.75)
        color = cmp(5,:);
        end
            if(pp_array(ic)>1.75 & pp_array(ic)<=2.0)
        color = cmp(6,:);
            end
        if(pp_array(ic)>2.0 )
        color = cmp(7,:);
        end
   % if(borehole_date(ic,2) < 1980)
   %    color = [0.7 0.7 0.7 ];
   % end
    rectangle('position',pos,'FaceColor',color)
    
end

set(gca, 'CLim', [0, 3.5]);

xlim([-150 160])
ylim([-50 77.5])


title('Prior hyper-parameter (\circC)')
box on
colorbar
pbaspect([2 1 1])
print -painters -depsc2 -r2500 hyperprior_glob.eps
