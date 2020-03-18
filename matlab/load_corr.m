
close all
clear



% Add in correction for SAT-GST offset from Melo_Auilar et al 2018, their 
% figure 9.


sat_gst_corr=ncread('Fig_9_data/all_mean_delT.nc','delT')  % This is the 2century change so divide by 2:

sat_gst_corr=sat_gst_corr/2.;

size(sat_gst_corr)

sat_gst_corr=sat_gst_corr  %  because this is SAT-GST
B=sat_gst_corr;
B(abs(sat_gst_corr) < 0.1)= NaN;

%sat_gst_corr=B;

imagesc(B)

figure

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




% --------------------------------------
% Now map plot the magnitude of warming
% --------------------------------------


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
figure
plot(partition_lonlat(1:nparts,1),corr(1:nparts),'linestyle','none','marker','o')
