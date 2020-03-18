
ntimes=2;
ts = [0:100:500];
ts=ts*365.*24.*60*60.;

Ts=zeros(ntimes,1);

Ts(2)=1.0;


ndepths=20;

ndepths=50;
z=[10:10:500];
flipud(z);

A = zeros(ndepths,ntimes);

kappa=1E-6;

for iz=1:ndepths
     for it=2:ntimes
         kappat = (kappa*ts(it))^0.5;
         kappatm1 = (kappa*ts(it-1))^0.5;
         A(iz,it) =  erfc(z(iz)/(2.0*kappat)) - erfc(z(iz)/(2.0*kappatm1));
         
     end
end

Tz=zeros(ndepths,1);
Tz(:)=0.0;
for iz=1:ndepths
     for it=2:ntimes
         Tz(iz) = Tz(iz) + Ts(it)*A(iz,it);
     end
end


plot(z,Tz')


% Now alter the thermal conductivity

A(:,:)=0.0;
Tz(:)=0.0;
kappa2=kappa*1.5;
        
for iz=1:ndepths
     for it=2:ntimes
         kappat = (kappa2*ts(it))^0.5;
         kappatm1 = (kappa2*ts(it-1))^0.5;
         A(iz,it) =  erfc(z(iz)/(2.0*kappat)) - erfc(z(iz)/(2.0*kappatm1));
         
     end
end

Tz=zeros(ndepths,1);
Tz(:)=0.0;
for iz=1:ndepths
     for it=2:ntimes
         Tz(iz) = Tz(iz) + Ts(it)*A(iz,it);
     end
end

hold all
plot(z,Tz') 
