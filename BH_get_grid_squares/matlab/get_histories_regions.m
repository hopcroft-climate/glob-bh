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

load colocbh_no.txt
load BH_locats.txt
load bhlocats.txt

% the number of gridcells that we are reading in:
nparts=173;

% The number of MCMC iterations to sample from:
itmax=10000;
    
    
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

% 6 regions: North America, Europe, Asia, S America, S Africa, Australia:
glob=zeros(123,6);
glob_p2s=zeros(123,6);
glob_m2s=zeros(123,6);
glob_p3s=zeros(123,6);
glob_m3s=zeros(123,6);
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

count=zeros(123,6);
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
    area = (A * A*((cos(latrad(ic) - grid_lat_rad) - cos(latrad(ic)+grid_lat_rad))*2*pi*5)/360);
    
    if(partition_lonlat(ic,2)>30 & partition_lonlat(ic,2)<60)
         if(partition_lonlat(ic,1)>210-360 & partition_lonlat(ic,1)<320-360)% Only North America
 
    for it=1:end_index(ic)
     % North America
  glob(it,1) = glob(it,1) + area *  (mean_G(ic,it) - clim6090);
  glob_p2s(it,1) = glob_p2s(it,1) + area * (mean_G(ic,it)- clim6090 + 2*stdev_G(ic,it)) ;
  glob_m2s(it,1) = glob_m2s(it,1) + area * (mean_G(ic,it)- clim6090 - 2*stdev_G(ic,it)) ;
  glob_p3s(it,1) = glob_p3s(it,1) + area * (mean_G(ic,it)- clim6090 + 3*stdev_G(ic,it)) ;
  glob_m3s(it,1) = glob_m3s(it,1) + area * (mean_G(ic,it)- clim6090 - 3*stdev_G(ic,it)) ;
  count(it,1) = count(it,1)+area;
 end
         end
    end
    
    % Europe
    if(partition_lonlat(ic,2)>30 & partition_lonlat(ic,2)<60)
         if(partition_lonlat(ic,1)>0 & partition_lonlat(ic,1)<50 )% Europe
 for it=1:end_index(ic)
     % Europe
  glob(it,2) = glob(it,2) + area *  (mean_G(ic,it) - clim6090);
  glob_p2s(it,2) = glob_p2s(it,2) + area * (mean_G(ic,it)- clim6090 + 2*stdev_G(ic,it)) ;
  glob_m2s(it,2) = glob_m2s(it,2) + area * (mean_G(ic,it)- clim6090 - 2*stdev_G(ic,it)) ;
  glob_p3s(it,2) = glob_p3s(it,2) + area * (mean_G(ic,it)- clim6090 + 3*stdev_G(ic,it)) ;
  glob_m3s(it,2) = glob_m3s(it,2) + area * (mean_G(ic,it)- clim6090 - 3*stdev_G(ic,it)) ;
  count(it,2) = count(it,2)+area;
 end
         end  % lat
    end % lon
           if(partition_lonlat(ic,2)>30 & partition_lonlat(ic,2)<60)
               if(partition_lonlat(ic,1)>-20 & partition_lonlat(ic,1)<=0 )% Europe
 for it=1:end_index(ic)
     % Europe
  glob(it,2) = glob(it,2) + area *  (mean_G(ic,it) - clim6090);
  glob_p2s(it,2) = glob_p2s(it,2) + area * (mean_G(ic,it)- clim6090 + 2*stdev_G(ic,it)) ;
  glob_m2s(it,2) = glob_m2s(it,2) + area * (mean_G(ic,it)- clim6090 - 2*stdev_G(ic,it)) ;
  glob_p3s(it,2) = glob_p3s(it,2) + area * (mean_G(ic,it)- clim6090 + 3*stdev_G(ic,it)) ;
  glob_m3s(it,2) = glob_m3s(it,2) + area * (mean_G(ic,it)- clim6090 - 3*stdev_G(ic,it)) ;
  count(it,2) = count(it,2)+area;
 end
         end  % lat
    end % lon
    
      % Asia
               if(partition_lonlat(ic,2)>30 & partition_lonlat(ic,2)<60)  % Asia
                   if(partition_lonlat(ic,1)>50 )% Asia
 for it=1:end_index(ic)
     % Europe
  glob(it,3) = glob(it,3) + area *  (mean_G(ic,it) - clim6090);
  glob_p2s(it,3) = glob_p2s(it,3) + area * (mean_G(ic,it)- clim6090 + 2*stdev_G(ic,it)) ;
  glob_m2s(it,3) = glob_m2s(it,3) + area * (mean_G(ic,it)- clim6090 - 2*stdev_G(ic,it)) ;
  glob_p3s(it,3) = glob_p3s(it,3) + area * (mean_G(ic,it)- clim6090 + 3*stdev_G(ic,it)) ;
  glob_m3s(it,3) = glob_m3s(it,3) + area * (mean_G(ic,it)- clim6090 - 3*stdev_G(ic,it)) ;
  count(it,3) = count(it,3)+area;
 end
         end  % lat
    end % lon

    
    
    
        % South America
               if(partition_lonlat(ic,2)<0)  
                   if(partition_lonlat(ic,1)>-100 & partition_lonlat(ic,1)<-30 )
 for it=1:end_index(ic)
  glob(it,4) = glob(it,4) + area *  (mean_G(ic,it) - clim6090);
  glob_p2s(it,4) = glob_p2s(it,4) + area * (mean_G(ic,it)- clim6090 + 2*stdev_G(ic,it)) ;
  glob_m2s(it,4) = glob_m2s(it,4) + area * (mean_G(ic,it)- clim6090 - 2*stdev_G(ic,it)) ;
  glob_p3s(it,4) = glob_p3s(it,4) + area * (mean_G(ic,it)- clim6090 + 3*stdev_G(ic,it)) ;
  glob_m3s(it,4) = glob_m3s(it,4) + area * (mean_G(ic,it)- clim6090 - 3*stdev_G(ic,it)) ;
  count(it,4) = count(it,4)+area;
 end
         end  % lat
    end % lon

    
    
     % South Africa
               if(partition_lonlat(ic,2)<0)  
                   if(partition_lonlat(ic,1)>0 & partition_lonlat(ic,1)<50 )
 for it=1:end_index(ic)
  glob(it,5) = glob(it,5) + area *  (mean_G(ic,it) - clim6090);
  glob_p2s(it,5) = glob_p2s(it,5) + area * (mean_G(ic,it)- clim6090 + 2*stdev_G(ic,it)) ;
  glob_m2s(it,5) = glob_m2s(it,5) + area * (mean_G(ic,it)- clim6090 - 2*stdev_G(ic,it)) ;
  glob_p3s(it,5) = glob_p3s(it,5) + area * (mean_G(ic,it)- clim6090 + 3*stdev_G(ic,it)) ;
  glob_m3s(it,5) = glob_m3s(it,5) + area * (mean_G(ic,it)- clim6090 - 3*stdev_G(ic,it)) ;
  count(it,5) = count(it,5)+area;
 end
         end  % lat
    end % lon
    
    
    
    
    % Australia
               if(partition_lonlat(ic,2)<0)  
                   if(partition_lonlat(ic,1)>100 & partition_lonlat(ic,1)<180 )
 for it=1:end_index(ic)
  glob(it,6) = glob(it,6) + area *  (mean_G(ic,it) - clim6090);
  glob_p2s(it,6) = glob_p2s(it,6) + area * (mean_G(ic,it)- clim6090 + 2*stdev_G(ic,it)) ;
  glob_m2s(it,6) = glob_m2s(it,6) + area * (mean_G(ic,it)- clim6090 - 2*stdev_G(ic,it)) ;
  glob_p3s(it,6) = glob_p3s(it,6) + area * (mean_G(ic,it)- clim6090 + 3*stdev_G(ic,it)) ;
  glob_m3s(it,6) = glob_m3s(it,6) + area * (mean_G(ic,it)- clim6090 - 3*stdev_G(ic,it)) ;
  count(it,6) = count(it,6)+area;
 end
         end  % lat
    end % lon
    
    
    end % if isnan(clim6090)
end


for ic=1:6
for it=1:123
 glob(it,ic) = glob(it,ic) /count(it,ic);
 glob_p2s(it,ic) = glob_p2s(it,ic) /count(it,ic);
 glob_m2s(it,ic) = glob_m2s(it,ic) /count(it,ic);
 glob_p3s(it,ic) = glob_p3s(it,ic) /count(it,ic);
 glob_m3s(it,ic) = glob_m3s(it,ic) /count(it,ic);
end

subplot(2,3,ic)
hold all
plot(yearsAD,glob(:,ic))
plot(yearsAD,glob_p2s(:,ic))
plot(yearsAD,glob_m2s(:,ic)) 
plot(yearsAD,glob_p3s(:,ic))
plot(yearsAD,glob_m3s(:,ic))
xlim([1500 2020]); 
end % ic


save Tregions_new.txt glob -ascii
save Tregion_new_p2s.txt glob_p2s -ascii
save Tregions_new_m2s.txt glob_m2s -ascii

save Tregions_new_p3s.txt glob_p3s -ascii
save Tregions_new_m3s.txt glob_m3s -ascii


