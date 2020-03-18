close
clear
%cd ../../Desktop/nhbh/
load inputs/colocbh_no.txt
load inputs/BH_locats.txt
load inputs/bhlocats.txt
%colocbh_no2=load('profiles_set2/colocbh_no.txt')
nparts=173;
bh_long = (colocbh_no(:,1)-0.5)*5;
bh_lat = 90 -  (colocbh_no(:,2)-0.5)*5;
%bh_long2 = (colocbh_no2(:,1)-0.5)*5;
%bh_lat2 = 90 -  (colocbh_no2(:,2)-0.5)*5;

nparts=nparts;

for ic=1:nparts
 if(bh_long(ic)>180)
     bh_long(ic)=bh_long(ic) -360;
 end
end


% The number of MCMC iterations to sample from:
    itmax=50000;

hold on
for ic=1:nparts
    pos = [bh_long(ic)-2.5 bh_lat(ic)-2.5 5 5];
    rectangle('position',pos)
end

load inputs/world.dat
plot(world(:,1),world(:,2),'color','black')

%plot(bh_long,bh_lat,'linestyle','none','marker','s','markersize',15,'color','blue')
plot(BH_locats(:,1),BH_locats(:,2),'linestyle','none','marker','.','markersize',10,'color','black','linewidth',1)
%hold on
%plot(BH_locats2(:,1),BH_locats2(:,2),'linestyle','none','marker','o','markersize',3,'color','black','linewidth',1)


xlim([-150 160])
ylim([-50 77.5])
box on
pbaspect([2 1 1])

print -painters -depsc2 -r2500 bh_locats_all.eps

%  Now load in the 96 new profiles:

%load profiles_set2/BH_locats.txt
%load profiles_set2/colocbh_no.txt
