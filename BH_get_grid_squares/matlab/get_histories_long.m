close
clear

load colocbh_no.txt
load BH_locats.txt
load bhlocats.txt

nparts=173;
% The number of MCMC iterations to sample from:
    itmax=500;
    
    
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
load world.dat
plot(world(:,1),world(:,2),'color','black')


ntimes=403;
root='Samples/';
extension='G.txt';

extensiont='T.txt';
extensions='S.txt';
extensionpp='PP.txt';
root2='X';
extensions2='G';
mean_G = zeros(nparts,ntimes);   % nparts partitions x ntimes timesteps (5yrs)
stdev_G = zeros(nparts,ntimes);   % nparts partitions x ntimes timesteps (5yrs)
times5yr = [5:5:2015];
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
yearsAD=[5:5:2015];
glob=zeros(ntimes,1);
glob_p2s=zeros(ntimes,1);
glob_m2s=zeros(ntimes,1);
glob_p3s=zeros(ntimes,1);
glob_m3s=zeros(ntimes,1);
SH=zeros(ntimes,1);
SH_p2s=zeros(ntimes,1);
SH_m2s=zeros(ntimes,1);
SH_p3s=zeros(ntimes,1);
SH_m3s=zeros(ntimes,1);
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
   depth_time = 0.3*((borehole_z0(ic,2))^2/1e-6)/(3600*24*365.25)
   borehole_date(ic,2) = borehole_date(ic,2)  -depth_time;
   end_index(ic) =   int16(( 0.5+  borehole_date(ic,2) -2015  + 2015) /5) 
   % Here deselect partitions with overly narrow prior limits
  % if (pp_array(ic) < 0.25)
  %     end_index(ic)=0;
  % end
   
end

count=zeros(ntimes,1);
countSH=zeros(ntimes,1);
A=6378185.0;
pi=3.141592654;
grid_lat_rad=2.5;
grid_lon_rad=2.5;
latrad=2*pi*partition_lonlat(:,2)/360.0;
lonrad=2*pi*partition_lonlat(:,1)/360.0;

f6090clim=(1980-1400)/5;
i6090clim=(1955-1400)/5;


for ic=1:nparts
    
    mean_G(ic,end_index(ic)+1:403) = NaN;
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
 
 if(partition_lonlat(ic,2)<0)  % Only extra-tropical sites
 area = (-A * A*((cos(latrad(ic) - grid_lat_rad) - cos(latrad(ic)+grid_lat_rad))*2*pi*5)/360);
 
 for it=1:end_index(ic)
  SH(it) = SH(it) + area *  mean_G(ic,it);
  SH_p2s(it) = SH_p2s(it) + area * (mean_G(ic,it) + 2*stdev_G(ic,it)) ;
  SH_m2s(it) = SH_m2s(it) + area * (mean_G(ic,it) - 2*stdev_G(ic,it)) ;
  SH_p3s(it) = SH_p3s(it) + area * (mean_G(ic,it) + 3*stdev_G(ic,it)) ;
  SH_m3s(it) = SH_m3s(it) + area * (mean_G(ic,it) - 3*stdev_G(ic,it)) ;
  countSH(it) = countSH(it)+area;
 end
 end
end

for ic=1:nparts
 mean_G(ic,end_index(ic)+1:403) = NaN;
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
for it=1:403
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



save Tglob_new_long.txt glob -ascii
save Tglob_new_p2s_long.txt glob_p2s -ascii
save Tglob_new_m2s_long.txt glob_m2s -ascii

save Tglob_new_p3s_long.txt glob_p3s -ascii
save Tglob_new_m3s_long.txt glob_m3s -ascii

save TSH_new_long.txt SH -ascii
save TSH_new_p2s_long.txt SH_p2s -ascii
save TSH_new_m2s_long.txt SH_m2s -ascii

save TSH_new_p3s_long.txt SH_p3s -ascii
save TSH_new_m3s_long.txt SH_m3s -ascii



% --------------------------------------
% Now map plot the magnitude of warming
% --------------------------------------


warming_6090_m_1800_1850=zeros(nparts,1);
for ic=1:nparts
    warming_6090_m_1800_1850(ic)= mean(mean_G(ic,112:115)) - mean(mean_G(ic,80:86));
end

sd_warming_6090_m_1800_1850=zeros(nparts,1);
for ic=1:nparts
    sd_warming_6090_m_1800_1850(ic)= sqrt((mean(stdev_G(ic,112:115))/2)^2 + (mean(stdev_G(ic,80:86))/2)^2);
end

% Make it the 2 s.d.
sd_warming_6090_m_1800_1850=sd_warming_6090_m_1800_1850*2;

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
     
    if(warming_6090_m_1800_1850(ic)>-1 & warming_6090_m_1800_1850(ic)<=-0.5)
        color = cmp(1,:);
    end
    if(warming_6090_m_1800_1850(ic)>-0.5 & warming_6090_m_1800_1850(ic)<=0)
        color = cmp(2,:);
    end
    if(warming_6090_m_1800_1850(ic)>0 & warming_6090_m_1800_1850(ic)<=0.5)
        color = cmp(3,:);
    end
    if(warming_6090_m_1800_1850(ic)>0.5 & warming_6090_m_1800_1850(ic)<=1.0)
        color = cmp(4,:);
    end
    if(warming_6090_m_1800_1850(ic)>1.0 & warming_6090_m_1800_1850(ic)<=1.5)
        color = cmp(5,:);
    end
    if(warming_6090_m_1800_1850(ic)>1.5 & warming_6090_m_1800_1850(ic)<=2.0)
        color = cmp(6,:);
    end
    if(warming_6090_m_1800_1850(ic)>2.0 )
        color = cmp(7,:);
    end
    if(borehole_date(ic,2) < 1980)
       color = [0.7 0.7 0.7 ];
    end
    rectangle('position',pos,'FaceColor',color)
    
end

set(gca, 'CLim', [-1, 2.5]);

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
    if(borehole_date(ic,2) < 1980)
       color = [0.7 0.7 0.7 ];
    end
    
        
    rectangle('position',pos,'FaceColor',color)
    
end

set(gca, 'CLim', [0, 0.35]);

xlim([-150 160])
ylim([-50 77.5])

title('2 s.d. temperature change (\circC)')
box on
pbaspect([2 1 1])
print -painters -depsc2 -r2500 dT_glob_long.eps


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
print -painters -depsc2 -r2500 hyperprior_glob_long.eps