
clear
close all
%close
%close
%cd Out_ver0
%load TGLOBAL.txt
%load T_NAM.txt
%load T_EU.txt
%load T_ASIA.txt


load inputs/nhex_crutemp_20yr.dat
load inputs/shex_crutemp_20yr.dat
years=[1850:1850+200-1];


times=[1405:5:2015];




f6090clim=(1980-1400)/5;
i6090clim=(1955-1400)/5;
iclim6090=i6090clim;
fclim6090=f6090clim;
% 1755-1800:
iLIAclim=71;
fLIAclim=iLIAclim+10;

% Read in the series calculated in MATLAB

yearsAD=[1405:5:2015];
load outputs_lcoef/Tglob_new.txt
load outputs_lcoef/Tglob_new_p2s.txt
load outputs_lcoef/Tglob_new_m2s.txt
load outputs_lcoef/Tglob_new_p3s.txt
load outputs_lcoef/Tglob_new_m3s.txt

load outputs_lcoef/TSH_new.txt
load outputs_lcoef/TSH_new_p2s.txt
load outputs_lcoef/TSH_new_m2s.txt
load outputs_lcoef/TSH_new_p3s.txt
load outputs_lcoef/TSH_new_m3s.txt

%Tglob_new=Tglob_new_all;
%Tglob_new_m2s=Tglob_new_m2s_all;
%Tglob_new_p2s=Tglob_new_p2s_all;
%Tglob_new_m3s=Tglob_new_m3s_all;
%Tglob_new_p3s=Tglob_new_p3s_all;

mean_TGLOBAL=Tglob_new;
TGLOB_p2s=Tglob_new_p2s;
TGLOB_m2s=Tglob_new_m2s;

TGLOB_p3s=Tglob_new_p3s;
TGLOB_m3s=Tglob_new_m3s;
 
mean_TSH=TSH_new;
TSH_p2s=TSH_new_p2s;
TSH_m2s=TSH_new_m2s;

TSH_p3s=TSH_new_p3s;
TSH_m3s=TSH_new_m3s; 
 
 

 


%1750-1800 to 6090 difference:
X30yr_amplitude = mean(mean_TGLOBAL(i6090clim:f6090clim)) - mean(mean_TGLOBAL(90:93));
X6090_2sd=mean(TGLOB_p2s(i6090clim:f6090clim,:)) - mean(mean_TGLOBAL(i6090clim:f6090clim));
X1850_2sd=mean(TGLOB_p2s(90:93,:)) - mean(mean_TGLOBAL(90:93));
X30yr_2sd = sqrt((X6090_2sd/2)*(X6090_2sd/2) + (X1850_2sd/2) *(X1850_2sd/2));
X30yr_2sd =X30yr_2sd *2;

X30yr_max_amplitude = mean(mean_TGLOBAL(i6090clim:f6090clim)) - mean(mean_TGLOBAL(iLIAclim:fLIAclim));
X1825_2sd=mean(TGLOB_p2s(iLIAclim:fLIAclim,:)) - mean(mean_TGLOBAL(iLIAclim:fLIAclim));
X30yr_2sd = sqrt((X6090_2sd/2)*(X6090_2sd/2) + (X1850_2sd/2) *(X1850_2sd/2));
X30yr_2sd =X30yr_2sd *2;


XSH30yr_amplitude = mean(mean_TSH(i6090clim:f6090clim)) - mean(mean_TSH(90:93));
XSH6090_2sd=mean(TSH_p2s(i6090clim:f6090clim,:)) - mean(mean_TSH(i6090clim:f6090clim));
XSH1850_2sd=mean(TSH_p2s(90:93,:)) - mean(mean_TGLOBAL(90:93));
XSH30yr_2sd = sqrt((XSH6090_2sd/2)*(XSH6090_2sd/2) + (XSH1850_2sd/2) *(XSH1850_2sd/2));
XSH30yr_2sd =XSH30yr_2sd *2;

XSH30yr_max_amplitude = mean(mean_TSH(i6090clim:f6090clim)) - mean(mean_TSH(iLIAclim:fLIAclim));
XSH1825_2sd=mean(TSH_p2s(iLIAclim:fLIAclim,:)) - mean(mean_TSH(iLIAclim:fLIAclim));
XSH30yr_2sd = sqrt((XSH6090_2sd/2)*(XSH6090_2sd/2) + (XSH1850_2sd/2) *(XSH1850_2sd/2));
XSH30yr_2sd =XSH30yr_2sd *2;

% Same in CRU:
  mean(nhex_crutemp_20yr(116:116+24-1)) - mean(nhex_crutemp_20yr(12:12+30-1))
 mean(shex_crutemp_20yr(116:116+24-1)) - mean(shex_crutemp_20yr(12:12+30-1))



iclim6090_new=iclim6090+2;
fclim6090_new=fclim6090+2;

subplot(2,1,1)
hold all
%hmean=plot(yearsAD(1:fclim6090=),mean_TGLOBAL(1:fclim6090,1)-mean(mean_TGLOBAL(iclim6090:fclim6090)),'color','blue','linewidth',2.0)
%plot(yearsAD(1:fclim6090),TGLOB_p2s(1:fclim6090)-mean(mean_TGLOBAL(iclim6090:fclim6090)),'color','blue')
%plot(yearsAD(1:fclim6090),TGLOB_m2s(1:fclim6090)-mean(mean_TGLOBAL(iclim6090:fclim6090)),'color','blue')
%plot(yearsAD(1:118),TGLOB_p3s(1:118)-mean(mean_TGLOBAL(iclim6090:fclim6090)),'color','blue')
%plot(yearsAD(1:118),TGLOB_m3s(1:118)-mean(mean_TGLOBAL(iclim6090:fclim6090)),'color','blue')

% Plot all years
hmean=plot(yearsAD,mean_TGLOBAL-mean(mean_TGLOBAL(iclim6090_new:fclim6090_new)),'color',[0.15 0.5 0.6],'linewidth',2.0)
plot(yearsAD,TGLOB_p2s-mean(mean_TGLOBAL(iclim6090_new:fclim6090_new)),'color',[0.15 0.5 0.6])
plot(yearsAD,TGLOB_m2s-mean(mean_TGLOBAL(iclim6090_new:fclim6090_new)),'color',[0.15 0.5 0.6])

% Plot all years
%hmean2=plot(yearsAD,mean_TGLOBAL_gIgcorr-mean(mean_TGLOBAL_gIgcorr(iclim6090_new:fclim6090_new)),'color',[0.8 0.5 0.6],'linewidth',2.0)
%plot(yearsAD,TGLOB_p2s_gIgcorr-mean(mean_TGLOBAL_gIgcorr(iclim6090_new:fclim6090_new)),'color',[0.8 0.5 0.6])
%plot(yearsAD,TGLOB_m2s_gIgcorr-mean(mean_TGLOBAL_gIgcorr(iclim6090_new:fclim6090_new)),'color',[0.8 0.5 0.6])



hcru=plot(years(12:150),nhex_crutemp_20yr(12:150),'color',[0.2 0.2 0.2],'linewidth',2.0)
%xlim([1500 2020])
%ylim([-1 1])
box on

% Approximate 5 deg Pollack Smerdon 2004:
%pollack04 = [ 1500 1600 1700 1800 1900 2000; 
%    -0.8  -0.75  -0.7  -0.6  -0.35  0.25] ;
    
% Replace with Hang et al 2000 for consistency with SH:
pollack04 = [ 1500 1600 1700 1800 1900 2000; 
    -1.05 -1.0  -0.9  -0.75  -0.55  0.0] ;
    
ipollack04 = interp1([1500:100:2000],[-1.05 -1.0  -0.9  -0.75  -0.55  0.0],[1501:2000]); 
offsetpollack04=mean(ipollack04(461:490))
pollack04(2,:)=pollack04(2,:)-offsetpollack04;
        

    
pollack04 = pollack04';
pollack04_corr = pollack04;
%pollack04_corr(6,2) = 0.0;
%pollack04_corr(6,1) = 1950;

%hpollack=plot(pollack04(:,1),pollack04(:,2),'color',[0 0.6 0.4],'linewidth',2.0)
hps04=plot(pollack04_corr(:,1),pollack04_corr(:,2),'color',[0.3 0.3 0.3],'linewidth',2.0)


%cd ../nhbh/AbrametalPAGES2k2016_Supp2_inputdata/
cd AbrametalPAGES2k2016_Supp2_inputdata/


load archive_p2k_regional_reconstructions.mat

pages2k  = zeros(504,1);
%pages2k(1:501) = (34.4E6)*archive_p2k_regional_reconstructions(1).data.rec(:);
pages2k(1:504) = pages2k(1:504)+(13.0E6)*archive_p2k_regional_reconstructions(2).data.rec(:);
pages2k(1:490) = pages2k(1:490)+(31.1E6)*archive_p2k_regional_reconstructions(3).data.rec(:);
pages2k(1:480) = pages2k(1:480)+(12.5E6)*archive_p2k_regional_reconstructions(4).data.rec(:);

pages2k(1:480) = pages2k(1:480)/( 13.0E6 + 31.1E6 + 12.5E6);
pages2k(481:490) = pages2k(481:490)/( 13.0E6 + 31.1E6 );
pages2k(491:504) = pages2k(491:504)/( 13.0E6 );

tpages2k=[1500:2003];
tpages2k=fliplr(tpages2k);
% now smooth with 25 year moving average:
%pages2k_25yr = smooth(pages2k,25)
pages2k_25yr = pages2k;


%plot(tpages2k(1:501),pages2k(1:501))
hpages = plot(tpages2k(13:489),pages2k_25yr(13:489),'linewidth',0.5,'color',[0.5 0.5 0.5])
%legend(hcru,'CRUTEM4')
legend([hcru,hpages,hps04,hmean],'CRUTEM4','PAGES 2k (2013)','Huang et al., 2000','Posterior mean \pm2 s.d.','location','northwest')



%plot(yearsAD,Tglob_new-0.4,'linewidth',2.0,'color','red')
%plot(yearsAD,Tglob_new_p2s-0.4,'linewidth',2.0,'color','red','linestyle','--')
%plot(yearsAD,Tglob_new_m2s-0.4,'linewidth',2.0,'color','red','linestyle','--')

xlabel('year CE')
ylabel('temperature (\circC) w.r.t. 1960-1990')
ylim([-1.0 0.5])
xlim([1500 2010])
pbaspect([2 1 1])
cd ..







pwd
cd ../glob_bh



subplot(2,1,2)
hold all
%hmean=plot(yearsAD(1:fclim6090),mean_TSH(1:fclim6090,1)-mean(mean_TSH(iclim6090:fclim6090)),'color','blue','linewidth',2.0)
%plot(yearsAD(1:fclim6090),TSH_p2s(1:fclim6090)-mean(mean_TSH(iclim6090:fclim6090)),'color','blue')
%plot(yearsAD(1:fclim6090),TSH_m2s(1:fclim6090)-mean(mean_TSH(iclim6090:fclim6090)),'color','blue')
%plot(yearsAD(1:118),TGLOB_p3s(1:118)-mean(mean_TGLOBAL(iclim6090:fclim6090)),'color','blue')
%plot(yearsAD(1:118),TGLOB_m3s(1:118)-mean(mean_TGLOBAL(iclim6090:fclim6090)),'color','blue')
% Plot all years
hmean=plot(yearsAD,mean_TSH-mean(mean_TSH(iclim6090_new:fclim6090_new)),'color',[0.15 0.5 0.6],'linewidth',2.0)
plot(yearsAD,TSH_p2s-mean(mean_TSH(iclim6090_new:fclim6090_new)),'color',[0.15 0.5 0.6])
plot(yearsAD,TSH_m2s-mean(mean_TSH(iclim6090_new:fclim6090_new)),'color',[0.15 0.5 0.6])

%hmean2=plot(yearsAD,mean_TSH_gIgcorr-mean(mean_TSH_gIgcorr(iclim6090_new:fclim6090_new)),'color',[0.9 0.5 0.6],'linewidth',2.0)
%plot(yearsAD,TSH_p2s_gIgcorr-mean(mean_TSH_gIgcorr(iclim6090_new:fclim6090_new)),'color',[0.9 0.5 0.6])
%plot(yearsAD,TSH_m2s_gIgcorr-mean(mean_TSH_gIgcorr(iclim6090_new:fclim6090_new)),'color',[0.9 0.5 0.6])

hcru=plot(years(12:150),shex_crutemp_20yr(12:150),'color',[0.2 0.2 0.2],'linewidth',2.0)
%xlim([1500 2020])
%ylim([-0.9 0.5])
box on

% Approximate 5 deg Pollack Smerdon 2004:
pollack04 = [ 1500 1600 1700 1800 1900 2000; 
    -0.8 -0.8  -0.75 -0.6  -0.4  0.0] ;
ipollack04 = interp1([1500:100:2000],[-0.8 -0.8  -0.75 -0.6  -0.4  0.0],[1501:2000]); 
offsetpollack04=mean(ipollack04(461:490))
pollack04(2,:)=pollack04(2,:)-offsetpollack04;
    
pollack04 = pollack04';
pollack04_corr = pollack04;
%pollack04_corr(6,2) = 0.0;
%pollack04_corr(6,1) = 1950;

%hpollack=plot(pollack04(:,1),pollack04(:,2),'color',[0 0.6 0.4],'linewidth',2.0)
plot(pollack04(:,1),pollack04_corr(:,2),'color',[0.3 0.3 0.3],'linewidth',2.0)


cd AbrametalPAGES2k2016_Supp2_inputdata/
load archive_p2k_regional_reconstructions.mat
pages2k  = zeros(504,1);
%pages2k(1:501) = (34.4E6)*archive_p2k_regional_reconstructions(1).data.rec(:);
pages2k(1:502) = pages2k(1:502)+37.9E6*archive_p2k_regional_reconstructions(8).data.rec(:);
%pages2k(1:496) = pages2k(1:496)+(20.0E6)*archive_p2k_regional_reconstructions(9).data.rec(:);
%pages2k(1:480) = pages2k(1:480)+(12.5E6)*archive_p2k_regional_reconstructions(4).data.rec(:);
pages2k(1:496) = pages2k(1:496)/( 37.9E6 );
pages2k(497:502) = pages2k(497:502)/( 37.9E6 );
tpages2k=[1500:2003];
tpages2k=fliplr(tpages2k);
% now smooth with 25 year moving average:
%pages2k_25yr = smooth(pages2k,20)
pages2k_25yr = pages2k;


tpages2k=[1500:2003];
tpages2k=fliplr(tpages2k);
% now smooth with 25 year moving average:
%pages2k_25yr = smooth(pages2k,25)
pages2k_25yr = pages2k;


%plot(tpages2k(1:501),pages2k(1:501))
hpages = plot(tpages2k(13:489),pages2k_25yr(13:489),'linewidth',0.5,'color',[0.5 0.5 0.5])
%legend(hcru,'CRUTEM4')
legend([hcru,hpages,hps04,hmean],'CRUTEM4','PAGES 2k (2013)','Huang et al., 2000','Posterior mean \pm2 s.d.','location','northwest')

cd ../../glob_bh
 
%plot(yearsAD,Tglob_new-0.4,'linewidth',2.0,'color','red')
%plot(yearsAD,Tglob_new_p2s-0.4,'linewidth',2.0,'color','red','linestyle','--')
%plot(yearsAD,Tglob_new_m2s-0.4,'linewidth',2.0,'color','red','linestyle','--')

xlabel('year CE')
ylabel('temperature (\circC) w.r.t. 1960-1990')
ylim([-1.0 0.5])
xlim([1500 2010])
pbaspect([2 1 1])

print -painters -depsc2 -r2500 tas_glob_lcoef.eps

