
%%%%%%%
% Inicialización puntos
%%%%%%% 

%%% semaforos=[40.337563 -3.756103; 40.334738 -3.757735; 40.340432 -3.756063; 40.336111 -3.758840];
%%% lat1=semaforos(1,1);
%%% long1=semaforos(1,2);
%%% lat2=semaforos(2,1);
%%% long2=semaforos(2,2);
%%% distsem=acos(sin(lat1*pi/180)*sin(lat2*pi/180)+cos(lat1*pi/180)*cos(lat2*pi/180)*cos(long2*pi/180-long2*pi/180))*6371000;

semaforos=[40.337495 -3.756176; 40.334676 -3.757981; 40.336179 -3.758749; 40.333329 -3.761056; 40.332658 -3.761433; 40.336871 -3.757811];
rotondas=[40.342896 -3.756548; 40.347349 -3.755183; 40.334535 -3.760159; 40.335031 -3.760151; 40.322657 -3.716250; 40.327287 -3.717826];
nullclass=[40.324589 -3.715099; 40.338980 -3.755624; 40.337382 -3.757131; 40.335848 -3.759195; 40.334031 -3.760639; 40.334733 -3.759064];
cruces=[40.324266 -3.714687; 40.322100 -3.716716; 40.332153 -3.762075; 40.326559 -3.718480; 40.332369 -3.762926; 40.321107 -3.717869];

clear vsem;
clear vrot;
clear vcru;
clear vnull;

%%%%%%%
% Semaforos
%%%%%%% 

vsem = zeros (1,1,61);
ini=0;

for i=100:length(Velocidad)-100
for s=1:length(semaforos)
lat1=semaforos(s,1);
long1=semaforos(s,2);
lat2=Latitud(i);
long2=Longitud(i);
distsem=acos(sin(lat1*pi/180)*sin(lat2*pi/180)+cos(lat1*pi/180)*cos(lat2*pi/180)*cos(long2*pi/180-long1*pi/180))*6371000;
if (distsem<5)
v=Velocidad(i-30:i+30);
if ini==0
vsem(s,length(vsem(1,:,1)),:)=v;
ini=1;
else
vsem(s,length(vsem(1,:,1))+1,:)=v;
end
end
end
end 

%%%matsem=zeros(100,length(semaforos));
clear matsem;
for i=1:length(semaforos)
tmp=0;
for j=1:length(vsem(1,:,1))
if (vsem(i,j,5)~=0)
if tmp==0
tmp=vsem(i,j,:);
tmp=tmp(:)';
else
tmp2=vsem(i,j,:);
tmp2=tmp2(:)';
tmp=[tmp;tmp2];
end
end
end
tmp=tmp';
tmp2=0;
a=size(tmp);
th=mean(mean(tmp))-mean(std(tmp));
for k=1:a(2)
if(mean(tmp(:,k))>th)
if tmp2==0
tmp2=tmp(31,k);
else
tmp2=[tmp2; tmp(31,k)];
end
end
end
tmp=-1*ones(100,1);
for l=1:length(tmp2);
if l<101
tmp(l)=tmp2(l);
end
end
matsem(:,i)=tmp;
end


%%%%%%%
% rotondas
%%%%%%% 

vrot = zeros (1,1,61);
ini=0;

for i=100:length(Velocidad)-100
for s=1:length(rotondas)
lat1=rotondas(s,1);
long1=rotondas(s,2);
lat2=Latitud(i);
long2=Longitud(i);
distrot=acos(sin(lat1*pi/180)*sin(lat2*pi/180)+cos(lat1*pi/180)*cos(lat2*pi/180)*cos(long2*pi/180-long1*pi/180))*6371000;
if (distrot<5)
v=Velocidad(i-30:i+30);
if ini==0
vrot(s,length(vrot(1,:,1)),:)=v;
ini=1;
else
vrot(s,length(vrot(1,:,1))+1,:)=v;
end
end
end
end 

clear matrot;
for i=1:length(rotondas)
tmp=0;
for j=1:length(vrot(1,:,1))
if (vrot(i,j,5)~=0)
if tmp==0
tmp=vrot(i,j,:);
tmp=tmp(:)';
else
tmp2=vrot(i,j,:);
tmp2=tmp2(:)';
tmp=[tmp;tmp2];
end
end
end
tmp=tmp';
tmp2=0;
a=size(tmp);
th=mean(mean(tmp))-mean(std(tmp));
for k=1:a(2)
if(mean(tmp(:,k))>th)
if tmp2==0
tmp2=tmp(31,k);
else
tmp2=[tmp2; tmp(31,k)];
end
end
end
tmp=-1*ones(100,1);
for l=1:length(tmp2)
if l<101
tmp(l)=tmp2(l);
end
end
matrot(:,i)=tmp;
end


%%%%%%%
% cruces
%%%%%%% 

vcru = zeros (1,1,61);
ini=0;

for i=100:length(Velocidad)-100
for s=1:length(cruces)
lat1=cruces(s,1);
long1=cruces(s,2);
lat2=Latitud(i);
long2=Longitud(i);
distcru=acos(sin(lat1*pi/180)*sin(lat2*pi/180)+cos(lat1*pi/180)*cos(lat2*pi/180)*cos(long2*pi/180-long1*pi/180))*6371000;
if (distcru<5)
v=Velocidad(i-30:i+30);
if ini==0
vcru(s,length(vcru(1,:,1)),:)=v;
ini=1;
else
vcru(s,length(vcru(1,:,1))+1,:)=v;
end
end
end
end 

clear matcru;
for i=1:length(cruces)
tmp=0;
for j=1:length(vcru(1,:,1))
if (vcru(i,j,5)~=0)
if tmp==0
tmp=vcru(i,j,:);
tmp=tmp(:)';
else
tmp2=vcru(i,j,:);
tmp2=tmp2(:)';
tmp=[tmp;tmp2];
end
end
end
tmp=tmp';
tmp2=0;
a=size(tmp);
th=mean(mean(tmp))-mean(std(tmp));
for k=1:a(2)
if(mean(tmp(:,k))>th)
if tmp2==0
tmp2=tmp(31,k);
else
tmp2=[tmp2; tmp(31,k)];
end
end
end
tmp=-1*ones(100,1);
for l=1:length(tmp2);
if(l<101)
tmp(l)=tmp2(l);
end
end
matcru(:,i)=tmp;
end


%%%%%%%
% nullclass
%%%%%%% 

vnull = zeros (1,1,61);
ini=0;

for i=100:length(Velocidad)-100
for s=1:length(nullclass)
lat1=nullclass(s,1);
long1=nullclass(s,2);
lat2=Latitud(i);
long2=Longitud(i);
distnull=acos(sin(lat1*pi/180)*sin(lat2*pi/180)+cos(lat1*pi/180)*cos(lat2*pi/180)*cos(long2*pi/180-long1*pi/180))*6371000;
if (distnull<10)
v=Velocidad(i-30:i+30);
if ini==0
vnull(s,length(vnull(1,:,1)),:)=v;
ini=1;
else
vnull(s,length(vnull(1,:,1))+1,:)=v;
end
end
end
end 

clear matnull;
for i=1:length(nullclass)
tmp=0;
for j=1:length(vnull(1,:,1))
if (vnull(i,j,5)~=0)
if tmp==0
tmp=vnull(i,j,:);
tmp=tmp(:)';
else
tmp2=vnull(i,j,:);
tmp2=tmp2(:)';
tmp=[tmp;tmp2];
end
end
end
tmp=tmp';
tmp2=0;
a=size(tmp);
th=mean(mean(tmp))-mean(std(tmp));
for k=1:a(2)
if(mean(tmp(:,k))>th)
if tmp2==0
tmp2=tmp(31,k);
else
tmp2=[tmp2; tmp(31,k)];
end
end
end
tmp=-1*ones(100,1);
for l=1:length(tmp2);
if(l<101)
tmp(l)=tmp2(l);
end
end
matnull(:,i)=tmp;
end


%%%%%%%%%%%%%
%%% pdfs
%%%%%%%%%%%%%

pdfsem=zeros(15,length(semaforos));
for j=1:length(semaforos)
for i=1:100
if (matsem(i,j)>=0 && matsem(i,j)<1)
pdfsem(1,j)=pdfsem(1,j)+1;
end
if (matsem(i,j)>=1 && matsem(i,j)<2)
pdfsem(2,j)=pdfsem(2,j)+1;
end
if (matsem(i,j)>=2 && matsem(i,j)<3)
pdfsem(3,j)=pdfsem(3,j)+1;
end
if (matsem(i,j)>=3 && matsem(i,j)<4)
pdfsem(4,j)=pdfsem(4,j)+1;
end
if (matsem(i,j)>=4 && matsem(i,j)<5)
pdfsem(5,j)=pdfsem(5,j)+1;
end
if (matsem(i,j)>=5 && matsem(i,j)<6)
pdfsem(6,j)=pdfsem(6,j)+1;
end
if (matsem(i,j)>=6 && matsem(i,j)<7)
pdfsem(7,j)=pdfsem(7,j)+1;
end
if (matsem(i,j)>=7 && matsem(i,j)<8)
pdfsem(8,j)=pdfsem(8,j)+1;
end
if (matsem(i,j)>=8 && matsem(i,j)<9)
pdfsem(9,j)=pdfsem(9,j)+1;
end
if (matsem(i,j)>=9 && matsem(i,j)<10)
pdfsem(10,j)=pdfsem(10,j)+1;
end
if (matsem(i,j)>=10 && matsem(i,j)<11)
pdfsem(11,j)=pdfsem(11,j)+1;
end
if (matsem(i,j)>=11 && matsem(i,j)<12)
pdfsem(12,j)=pdfsem(12,j)+1;
end
if (matsem(i,j)>=12 && matsem(i,j)<13)
pdfsem(13,j)=pdfsem(13,j)+1;
end
if (matsem(i,j)>=13 && matsem(i,j)<14)
pdfsem(14,j)=pdfsem(14,j)+1;
end
if (matsem(i,j)>=14 && matsem(i,j)<15)
pdfsem(15,j)=pdfsem(15,j)+1;
end
end
suma=sum(pdfsem(:,j));
for k=1:15
pdfsem(k,j)=pdfsem(k,j)/suma;
end
end


pdfrot=zeros(15,length(rotondas));
for j=1:length(rotondas)
for i=1:100
if (matrot(i,j)>=0 && matrot(i,j)<1)
pdfrot(1,j)=pdfrot(1,j)+1;
end
if (matrot(i,j)>=1 && matrot(i,j)<2)
pdfrot(2,j)=pdfrot(2,j)+1;
end
if (matrot(i,j)>=2 && matrot(i,j)<3)
pdfrot(3,j)=pdfrot(3,j)+1;
end
if (matrot(i,j)>=3 && matrot(i,j)<4)
pdfrot(4,j)=pdfrot(4,j)+1;
end
if (matrot(i,j)>=4 && matrot(i,j)<5)
pdfrot(5,j)=pdfrot(5,j)+1;
end
if (matrot(i,j)>=5 && matrot(i,j)<6)
pdfrot(6,j)=pdfrot(6,j)+1;
end
if (matrot(i,j)>=6 && matrot(i,j)<7)
pdfrot(7,j)=pdfrot(7,j)+1;
end
if (matrot(i,j)>=7 && matrot(i,j)<8)
pdfrot(8,j)=pdfrot(8,j)+1;
end
if (matrot(i,j)>=8 && matrot(i,j)<9)
pdfrot(9,j)=pdfrot(9,j)+1;
end
if (matrot(i,j)>=9 && matrot(i,j)<10)
pdfrot(10,j)=pdfrot(10,j)+1;
end
if (matrot(i,j)>=10 && matrot(i,j)<11)
pdfrot(11,j)=pdfrot(11,j)+1;
end
if (matrot(i,j)>=11 && matrot(i,j)<12)
pdfrot(12,j)=pdfrot(12,j)+1;
end
if (matrot(i,j)>=12 && matrot(i,j)<13)
pdfrot(13,j)=pdfrot(13,j)+1;
end
if (matrot(i,j)>=13 && matrot(i,j)<14)
pdfrot(14,j)=pdfrot(14,j)+1;
end
if (matrot(i,j)>=14 && matrot(i,j)<15)
pdfrot(15,j)=pdfrot(15,j)+1;
end
end
suma=sum(pdfrot(:,j));
for k=1:15
pdfrot(k,j)=pdfrot(k,j)/suma;
end
end


pdfcru=zeros(15,length(cruces));
for j=1:length(cruces)
for i=1:100
if (matcru(i,j)>=0 && matcru(i,j)<1)
pdfcru(1,j)=pdfcru(1,j)+1;
end
if (matcru(i,j)>=1 && matcru(i,j)<2)
pdfcru(2,j)=pdfcru(2,j)+1;
end
if (matcru(i,j)>=2 && matcru(i,j)<3)
pdfcru(3,j)=pdfcru(3,j)+1;
end
if (matcru(i,j)>=3 && matcru(i,j)<4)
pdfcru(4,j)=pdfcru(4,j)+1;
end
if (matcru(i,j)>=4 && matcru(i,j)<5)
pdfcru(5,j)=pdfcru(5,j)+1;
end
if (matcru(i,j)>=5 && matcru(i,j)<6)
pdfcru(6,j)=pdfcru(6,j)+1;
end
if (matcru(i,j)>=6 && matcru(i,j)<7)
pdfcru(7,j)=pdfcru(7,j)+1;
end
if (matcru(i,j)>=7 && matcru(i,j)<8)
pdfcru(8,j)=pdfcru(8,j)+1;
end
if (matcru(i,j)>=8 && matcru(i,j)<9)
pdfcru(9,j)=pdfcru(9,j)+1;
end
if (matcru(i,j)>=9 && matcru(i,j)<10)
pdfcru(10,j)=pdfcru(10,j)+1;
end
if (matcru(i,j)>=10 && matcru(i,j)<11)
pdfcru(11,j)=pdfcru(11,j)+1;
end
if (matcru(i,j)>=11 && matcru(i,j)<12)
pdfcru(12,j)=pdfcru(12,j)+1;
end
if (matcru(i,j)>=12 && matcru(i,j)<13)
pdfcru(13,j)=pdfcru(13,j)+1;
end
if (matcru(i,j)>=13 && matcru(i,j)<14)
pdfcru(14,j)=pdfcru(14,j)+1;
end
if (matcru(i,j)>=14 && matcru(i,j)<15)
pdfcru(15,j)=pdfcru(15,j)+1;
end
end
suma=sum(pdfcru(:,j));
for k=1:15
pdfcru(k,j)=pdfcru(k,j)/suma;
end
end


pdfnull=zeros(15,length(nullclass));
for j=1:length(nullclass)
for i=1:100
if (matnull(i,j)>=0 && matnull(i,j)<1)
pdfnull(1,j)=pdfnull(1,j)+1;
end
if (matnull(i,j)>=1 && matnull(i,j)<2)
pdfnull(2,j)=pdfnull(2,j)+1;
end
if (matnull(i,j)>=2 && matnull(i,j)<3)
pdfnull(3,j)=pdfnull(3,j)+1;
end
if (matnull(i,j)>=3 && matnull(i,j)<4)
pdfnull(4,j)=pdfnull(4,j)+1;
end
if (matnull(i,j)>=4 && matnull(i,j)<5)
pdfnull(5,j)=pdfnull(5,j)+1;
end
if (matnull(i,j)>=5 && matnull(i,j)<6)
pdfnull(6,j)=pdfnull(6,j)+1;
end
if (matnull(i,j)>=6 && matnull(i,j)<7)
pdfnull(7,j)=pdfnull(7,j)+1;
end
if (matnull(i,j)>=7 && matnull(i,j)<8)
pdfnull(8,j)=pdfnull(8,j)+1;
end
if (matnull(i,j)>=8 && matnull(i,j)<9)
pdfnull(9,j)=pdfnull(9,j)+1;
end
if (matnull(i,j)>=9 && matnull(i,j)<10)
pdfnull(10,j)=pdfnull(10,j)+1;
end
if (matnull(i,j)>=10 && matnull(i,j)<11)
pdfnull(11,j)=pdfnull(11,j)+1;
end
if (matnull(i,j)>=11 && matnull(i,j)<12)
pdfnull(12,j)=pdfnull(12,j)+1;
end
if (matnull(i,j)>=12 && matnull(i,j)<13)
pdfnull(13,j)=pdfnull(13,j)+1;
end
if (matnull(i,j)>=13 && matnull(i,j)<14)
pdfnull(14,j)=pdfnull(14,j)+1;
end
if (matnull(i,j)>=14 && matnull(i,j)<15)
pdfnull(15,j)=pdfnull(15,j)+1;
end
end
suma=sum(pdfnull(:,j));
for k=1:15
pdfnull(k,j)=pdfnull(k,j)/suma;
end
end

figure(1)
plot(pdfsem);
figure(2)
plot(pdfrot);
figure(3)
plot(pdfcru);
figure(4)
plot(pdfnull);

figure(5)
plot(mean(pdfsem'));
figure(6)
plot(mean(pdfrot'));
figure(7)
plot(mean(pdfcru'));
figure(8)
plot(mean(pdfnull'));


%%%%% filtrado pdfsem

%%% pdf=pdfsem;
%%% clear pdfsem;
%%% pdfsem(:,1)=pdf(:,1);
%%% pdfsem(:,2)=pdf(:,4);
%%% pdfsem(:,3)=pdf(:,5);
%%% pdfsem(:,4)=pdf(:,6);
%%% pdfsem(:,5)=pdf(:,7);



%%%%%%%%%% 
%%%%% validaciones con total variations pdfs 
%%%%%%%%%% 

%%%%% semaforos

clear tbss;
for i=1:length(pdfsem(1,:))
for j=1:length(pdfsem(1,:))
tbss(i,j)=sum(abs(pdfsem(:,i)-pdfsem(:,j)));
end
end

clear tbsc;
for i=1:length(pdfsem(1,:))
for j=1:length(pdfcru(1,:))
tbsc(i,j)=sum(abs(pdfsem(:,i)-pdfcru(:,j)));
end
end

clear tbsr;
for i=1:length(pdfsem(1,:))
for j=1:length(pdfrot(1,:))
tbsr(i,j)=sum(abs(pdfsem(:,i)-pdfrot(:,j)));
end
end

clear tbsn;
for i=1:length(pdfsem(1,:))
for j=1:length(pdfnull(1,:))
tbsn(i,j)=sum(abs(pdfsem(:,i)-pdfnull(:,j)));
end
end

%%%%% cruces

clear tbcc;
for i=1:length(pdfcru(1,:))
for j=1:length(pdfcru(1,:))
tbcc(i,j)=sum(abs(pdfcru(:,i)-pdfcru(:,j)));
end
end

clear tbcs;
for i=1:length(pdfcru(1,:))
for j=1:length(pdfsem(1,:))
tbcs(i,j)=sum(abs(pdfcru(:,i)-pdfsem(:,j)));
end
end

clear tbcr;
for i=1:length(pdfcru(1,:))
for j=1:length(pdfrot(1,:))
tbcr(i,j)=sum(abs(pdfcru(:,i)-pdfrot(:,j)));
end
end

clear tbcn;
for i=1:length(pdfcru(1,:))
for j=1:length(pdfnull(1,:))
tbcn(i,j)=sum(abs(pdfcru(:,i)-pdfnull(:,j)));
end
end

%%%%% rotondas

clear tbrr;
for i=1:length(pdfrot(1,:))
for j=1:length(pdfrot(1,:))
tbrr(i,j)=sum(abs(pdfrot(:,i)-pdfrot(:,j)));
end
end

clear tbrs;
for i=1:length(pdfrot(1,:))
for j=1:length(pdfsem(1,:))
tbrs(i,j)=sum(abs(pdfrot(:,i)-pdfsem(:,j)));
end
end

clear tbrc;
for i=1:length(pdfrot(1,:))
for j=1:length(pdfcru(1,:))
tbrc(i,j)=sum(abs(pdfrot(:,i)-pdfcru(:,j)));
end
end

clear tbrn;
for i=1:length(pdfrot(1,:))
for j=1:length(pdfnull(1,:))
tbrn(i,j)=sum(abs(pdfrot(:,i)-pdfnull(:,j)));
end
end

%%%%% null

clear tbnn;
for i=1:length(pdfnull(1,:))
for j=1:length(pdfnull(1,:))
tbnn(i,j)=sum(abs(pdfnull(:,i)-pdfnull(:,j)));
end
end

clear tbns;
for i=1:length(pdfnull(1,:))
for j=1:length(pdfsem(1,:))
tbns(i,j)=sum(abs(pdfnull(:,i)-pdfsem(:,j)));
end
end

clear tbnc;
for i=1:length(pdfnull(1,:))
for j=1:length(pdfcru(1,:))
tbnc(i,j)=sum(abs(pdfnull(:,i)-pdfcru(:,j)));
end
end

clear tbnr;
for i=1:length(pdfnull(1,:))
for j=1:length(pdfrot(1,:))
tbnr(i,j)=sum(abs(pdfnull(:,i)-pdfrot(:,j)));
end
end

%%% best mean

stbss=sort(tbss');
stbsr=sort(tbsr');
stbsc=sort(tbsc');
stbsn=sort(tbsn');
stbrr=sort(tbrr');
stbrs=sort(tbrs');
stbrc=sort(tbrc');
stbrn=sort(tbrn');
stbcs=sort(tbcs');
stbcr=sort(tbcr');
stbcc=sort(tbcc');
stbcn=sort(tbcn');
stbns=sort(tbns');
stbnr=sort(tbnr');
stbnc=sort(tbnc');
stbnn=sort(tbnn');



clear decsem;
for j=1:length(pdfsem(1,:))
m2=mean(stbss(1:length(stbss(:,1)),j));
decsem(j)=1;
if mean(stbsr(:,j))<m2
decsem(j)=2;
end
if mean(stbsc(:,j))<m2
decsem(j)=3;
end
if mean(stbsn(:,j))<m2
decsem(j)=4;
end
end

clear decrot;
for j=1:length(pdfrot(1,:))
m2=mean(stbrr(1:length(stbrr(:,1)),j));
decrot(j)=2;
if mean(stbrs(:,j))<m2
decrot(j)=1;
end
if mean(stbrc(:,j))<m2
decrot(j)=3;
end
if mean(stbrn(:,j))<m2
decrot(j)=4;
end
end

clear deccru;
for j=1:length(pdfcru(1,:))
m2=mean(stbcc(1:length(stbcc(:,1)),j));
deccru(j)=3;
if mean(stbcs(:,j))<m2
deccru(j)=1;
end
if mean(stbcr(:,j))<m2
deccru(j)=2;
end
if mean(stbcn(:,j))<m2
deccru(j)=4;
end
end

clear decnull;
for j=1:length(pdfnull(1,:))
m2=mean(stbnn(1:length(stbnn(:,1)),j));
decnull(j)=4;
if mean(stbns(:,j))<m2
decnull(j)=1;
end
if mean(stbnr(:,j))<m2
decnull(j)=2;
end
if mean(stbnc(:,j))<m2
decnull(j)=3;
end
end

decsem
decrot
deccru
decnull


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% momentos

X=[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14];

msem=zeros(1,length(pdfsem(1,:)));
m2sem=zeros(1,length(pdfsem(1,:)));
m3sem=zeros(1,length(pdfsem(1,:)));
m4sem=zeros(1,length(pdfsem(1,:)));

for i=1:length(pdfsem(1,:))
msem(i)=X*pdfsem(:,i);
m2sem(i)=((X-msem(i)).^2)*pdfsem(:,i);
m3sem(i)=((X-msem(i)).^3)*pdfsem(:,i);
m4sem(i)=((X-msem(i)).^4)*pdfsem(:,i);
stdsem(i)=sqrt(m2sem(i));
skewnesssem(i)=m3sem(i)/(stdsem(i)^3);
kurtosissem(i)=m4sem(i)/(stdsem(i)^4);
end

mrot=zeros(1,length(pdfrot(1,:)));
m2rot=zeros(1,length(pdfrot(1,:)));
m3rot=zeros(1,length(pdfrot(1,:)));
m4rot=zeros(1,length(pdfrot(1,:)));

for i=1:length(pdfrot(1,:))
mrot(i)=X*pdfrot(:,i);
m2rot(i)=((X-mrot(i)).^2)*pdfrot(:,i);
m3rot(i)=((X-mrot(i)).^3)*pdfrot(:,i);
m4rot(i)=((X-mrot(i)).^4)*pdfrot(:,i);
stdrot(i)=sqrt(m2rot(i));
skewnessrot(i)=m3rot(i)/(stdrot(i)^3);
kurtosisrot(i)=m4rot(i)/(stdrot(i)^4);
end

mcru=zeros(1,length(pdfcru(1,:)));
m2cru=zeros(1,length(pdfcru(1,:)));
m3cru=zeros(1,length(pdfcru(1,:)));
m4cru=zeros(1,length(pdfcru(1,:)));

for i=1:length(pdfcru(1,:))
mcru(i)=X*pdfcru(:,i);
m2cru(i)=((X-mcru(i)).^2)*pdfcru(:,i);
m3cru(i)=((X-mcru(i)).^3)*pdfcru(:,i);
m4cru(i)=((X-mcru(i)).^4)*pdfcru(:,i);
stdcru(i)=sqrt(m2cru(i));
skewnesscru(i)=m3cru(i)/(stdcru(i)^3);
kurtosiscru(i)=m4cru(i)/(stdcru(i)^4);
end

mnull=zeros(1,length(pdfnull(1,:)));
m2null=zeros(1,length(pdfnull(1,:)));
m3null=zeros(1,length(pdfnull(1,:)));
m4null=zeros(1,length(pdfnull(1,:)));

for i=1:length(pdfnull(1,:))
mnull(i)=X*pdfnull(:,i);
m2null(i)=((X-mnull(i)).^2)*pdfnull(:,i);
m3null(i)=((X-mnull(i)).^3)*pdfnull(:,i);
m4null(i)=((X-mnull(i)).^4)*pdfnull(:,i);
stdnull(i)=sqrt(m2null(i));
skewnessnull(i)=m3null(i)/(stdnull(i)^3);
kurtosisnull(i)=m4null(i)/(stdnull(i)^4);
end

datasem=[msem' stdsem' skewnesssem' kurtosissem' ones(length(pdfnull(1,:)),1)];
datarot=[mrot' stdrot' skewnessrot' kurtosisrot' 2*ones(length(pdfnull(1,:)),1)];
datacru=[mcru' stdcru' skewnesscru' kurtosiscru' 3*ones(length(pdfnull(1,:)),1)];
datanull=[mnull' stdnull' skewnessnull' kurtosisnull' 4*ones(length(pdfnull(1,:)),1)];

dataset= [datasem; datarot; datacru; datanull]




















%%%%%%%%%% 
%%%%% validaciones con total variations pdfs en bloques de 5 m/s
%%%%%%%%%% 

lpdfsem(1,:)=sum(pdfsem(1:5,:));
lpdfsem(2,:)=sum(pdfsem(6:10,:));
lpdfsem(3,:)=sum(pdfsem(11:15,:));
lpdfcru(1,:)=sum(pdfcru(1:5,:));
lpdfcru(2,:)=sum(pdfcru(6:10,:));
lpdfcru(3,:)=sum(pdfcru(11:15,:));
lpdfrot(1,:)=sum(pdfrot(1:5,:));
lpdfrot(2,:)=sum(pdfrot(6:10,:));
lpdfrot(3,:)=sum(pdfrot(11:15,:));
lpdfnull(1,:)=sum(pdfnull(1:5,:));
lpdfnull(2,:)=sum(pdfnull(6:10,:));
lpdfnull(3,:)=sum(pdfnull(11:15,:));

figure(1)
plot(lpdfsem)
figure(2)
plot(lpdfrot)
figure(3)
plot(lpdfcru)
figure(4)
plot(lpdfnull)

%%%%% semaforos

clear tbss;
for i=1:length(lpdfsem(1,:))
for j=1:length(lpdfsem(1,:))
tbss(i,j)=sum(abs(lpdfsem(:,i)-lpdfsem(:,j)));
end
end

clear tbsc;
for i=1:length(lpdfsem(1,:))
for j=1:length(lpdfcru(1,:))
tbsc(i,j)=sum(abs(lpdfsem(:,i)-lpdfcru(:,j)));
end
end

clear tbsr;
for i=1:length(lpdfsem(1,:))
for j=1:length(lpdfrot(1,:))
tbsr(i,j)=sum(abs(lpdfsem(:,i)-lpdfrot(:,j)));
end
end

clear tbsn;
for i=1:length(lpdfsem(1,:))
for j=1:length(lpdfnull(1,:))
tbsn(i,j)=sum(abs(lpdfsem(:,i)-lpdfnull(:,j)));
end
end

%%%%% cruces

clear tbcc;
for i=1:length(lpdfcru(1,:))
for j=1:length(lpdfcru(1,:))
tbcc(i,j)=sum(abs(lpdfcru(:,i)-lpdfcru(:,j)));
end
end

clear tbcs;
for i=1:length(lpdfcru(1,:))
for j=1:length(lpdfsem(1,:))
tbcs(i,j)=sum(abs(lpdfcru(:,i)-lpdfsem(:,j)));
end
end

clear tbcr;
for i=1:length(lpdfcru(1,:))
for j=1:length(lpdfrot(1,:))
tbcr(i,j)=sum(abs(lpdfcru(:,i)-lpdfrot(:,j)));
end
end

clear tbcn;
for i=1:length(lpdfcru(1,:))
for j=1:length(lpdfnull(1,:))
tbcn(i,j)=sum(abs(lpdfcru(:,i)-lpdfnull(:,j)));
end
end

%%%%% rotondas

clear tbrr;
for i=1:length(lpdfrot(1,:))
for j=1:length(lpdfrot(1,:))
tbrr(i,j)=sum(abs(lpdfrot(:,i)-lpdfrot(:,j)));
end
end

clear tbrs;
for i=1:length(lpdfrot(1,:))
for j=1:length(lpdfsem(1,:))
tbrs(i,j)=sum(abs(lpdfrot(:,i)-lpdfsem(:,j)));
end
end

clear tbrc;
for i=1:length(lpdfrot(1,:))
for j=1:length(lpdfcru(1,:))
tbrc(i,j)=sum(abs(lpdfrot(:,i)-lpdfcru(:,j)));
end
end

clear tbrn;
for i=1:length(lpdfrot(1,:))
for j=1:length(lpdfnull(1,:))
tbrn(i,j)=sum(abs(lpdfrot(:,i)-lpdfnull(:,j)));
end
end

%%%%% null

clear tbnn;
for i=1:length(lpdfnull(1,:))
for j=1:length(lpdfnull(1,:))
tbnn(i,j)=sum(abs(lpdfnull(:,i)-lpdfnull(:,j)));
end
end

clear tbns;
for i=1:length(lpdfnull(1,:))
for j=1:length(lpdfsem(1,:))
tbns(i,j)=sum(abs(lpdfnull(:,i)-lpdfsem(:,j)));
end
end

clear tbnc;
for i=1:length(lpdfnull(1,:))
for j=1:length(lpdfcru(1,:))
tbnc(i,j)=sum(abs(lpdfnull(:,i)-lpdfcru(:,j)));
end
end

clear tbnr;
for i=1:length(lpdfnull(1,:))
for j=1:length(lpdfrot(1,:))
tbnr(i,j)=sum(abs(lpdfnull(:,i)-lpdfrot(:,j)));
end
end

%%% k-NN k=3

stbss=sort(tbss');
stbsr=sort(tbsr');
stbsc=sort(tbsc');
stbsn=sort(tbsn');

clear decsem;
for j=1:length(lpdfsem(1,:))
m1=stbss(2,j);
m2=stbss(3,j);
decsem(j)=1;
if stbsr(2,j)<m2
decsem(j)=2;
end
if stbsc(2,j)<m2
decsem(j)=3;
end
if stbsn(2,j)<m2
decsem(j)=4;
end
end













%%%%%%%%%% 
%%%%% validaciones con distancias entre cdfs --> no funciona bien
%%%%%%%%%% 

for j=1:length(pdfsem(1,:))
for i=1:length(pdfsem(:,1))
cdfsem(i,j)=sum(pdfsem(1:i,j));
end
end

for j=1:length(pdfrot(1,:))
for i=1:length(pdfrot(:,1))
cdfrot(i,j)=sum(pdfrot(1:i,j));
end
end

for j=1:length(pdfcru(1,:))
for i=1:length(pdfcru(:,1))
cdfcru(i,j)=sum(pdfcru(1:i,j));
end
end

for j=1:length(pdfnull(1,:))
for i=1:length(pdfnull(:,1))
cdfnull(i,j)=sum(pdfnull(1:i,j));
end
end

%%%%% semaforos

clear tbss;
for i=1:length(cdfsem(1,:))
for j=1:length(cdfsem(1,:))
tbss(i,j)=sum(abs(cdfsem(:,i)-cdfsem(:,j)));
end
end

clear tbsc;
for i=1:length(cdfsem(1,:))
for j=1:length(cdfcru(1,:))
tbsc(i,j)=sum(abs(cdfsem(:,i)-cdfcru(:,j)));
end
end

clear tbsr;
for i=1:length(cdfsem(1,:))
for j=1:length(cdfrot(1,:))
tbsr(i,j)=sum(abs(cdfsem(:,i)-cdfrot(:,j)));
end
end

clear tbsn;
for i=1:length(cdfsem(1,:))
for j=1:length(cdfnull(1,:))
tbsn(i,j)=sum(abs(cdfsem(:,i)-cdfnull(:,j)));
end
end

%%%%% cruces

clear tbcc;
for i=1:length(cdfcru(1,:))
for j=1:length(cdfcru(1,:))
tbcc(i,j)=sum(abs(cdfcru(:,i)-cdfcru(:,j)));
end
end

clear tbcs;
for i=1:length(cdfcru(1,:))
for j=1:length(cdfsem(1,:))
tbcs(i,j)=sum(abs(cdfcru(:,i)-cdfsem(:,j)));
end
end

clear tbcr;
for i=1:length(cdfcru(1,:))
for j=1:length(cdfrot(1,:))
tbcr(i,j)=sum(abs(cdfcru(:,i)-cdfrot(:,j)));
end
end

clear tbcn;
for i=1:length(cdfcru(1,:))
for j=1:length(cdfnull(1,:))
tbcn(i,j)=sum(abs(cdfcru(:,i)-cdfnull(:,j)));
end
end

%%%%% rotondas

clear tbrr;
for i=1:length(cdfrot(1,:))
for j=1:length(cdfrot(1,:))
tbrr(i,j)=sum(abs(cdfrot(:,i)-cdfrot(:,j)));
end
end

clear tbrs;
for i=1:length(cdfrot(1,:))
for j=1:length(cdfsem(1,:))
tbrs(i,j)=sum(abs(cdfrot(:,i)-cdfsem(:,j)));
end
end

clear tbrc;
for i=1:length(cdfrot(1,:))
for j=1:length(cdfcru(1,:))
tbrc(i,j)=sum(abs(cdfrot(:,i)-cdfcru(:,j)));
end
end

clear tbrn;
for i=1:length(cdfrot(1,:))
for j=1:length(cdfnull(1,:))
tbrn(i,j)=sum(abs(cdfrot(:,i)-cdfnull(:,j)));
end
end

%%%%% null

clear tbnn;
for i=1:length(cdfnull(1,:))
for j=1:length(cdfnull(1,:))
tbnn(i,j)=sum(abs(cdfnull(:,i)-cdfnull(:,j)));
end
end

clear tbns;
for i=1:length(cdfnull(1,:))
for j=1:length(cdfsem(1,:))
tbns(i,j)=sum(abs(cdfnull(:,i)-cdfsem(:,j)));
end
end

clear tbnc;
for i=1:length(cdfnull(1,:))
for j=1:length(cdfcru(1,:))
tbnc(i,j)=sum(abs(cdfnull(:,i)-cdfcru(:,j)));
end
end

clear tbnr;
for i=1:length(cdfnull(1,:))
for j=1:length(cdfrot(1,:))
tbnr(i,j)=sum(abs(cdfnull(:,i)-cdfrot(:,j)));
end
end

%%% k-NN k=3

stbss=sort(tbss');
stbsr=sort(tbsr');
stbsc=sort(tbsc');
stbsn=sort(tbsn');

clear decsem;
for j=1:length(cdfsem(1,:))
m1=stbss(2,j);
m2=stbss(3,j);
decsem(j)=1;
if stbsr(2,j)<m2
decsem(j)=2;
end
if stbsc(2,j)<m2
decsem(j)=3;
end
if stbsn(2,j)<m2
decsem(j)=4;
end
end














%%%%%%%%%%%%%
%%% old
%%%%%%%%%%%%%

for s=1:length(rotondas)
lat1=rotondas(s,1);
long1=rotondas(s,2);
lat2=Latitud(i);
long2=Longitud(i);
distsem=acos(sin(lat1*pi/180)*sin(lat2*pi/180)+cos(lat1*pi/180)*cos(lat2*pi/180)*cos(long2*pi/180-long1*pi/180))*6371000;
if (distsem<50)
v=Velocidad(i-20:i+20);
a=Aceleracin(i-20:i+20);
matrizrotv(:,rc)=v;
matrizrota(:,rc)=a;
crot(1:3,rc)=[lat2; long2; i];
rc=rc+1;
end
end

for s=1:length(cruces)
lat1=cruces(s,1);
long1=cruces(s,2);
lat2=Latitud(i);
long2=Longitud(i);
distsem=acos(sin(lat1*pi/180)*sin(lat2*pi/180)+cos(lat1*pi/180)*cos(lat2*pi/180)*cos(long2*pi/180-long1*pi/180))*6371000;
if (distsem<50)
v=Velocidad(i-20:i+20);
a=Aceleracin(i-20:i+20);
matrizcruv(:,zc)=v;
matrizcrua(:,zc)=a;
ccru(1:3,zc)=[lat2; long2; i];
zc=zc+1;
end
end

end







atipicos=zeros(length(Velocidad),1);
clear v;
clear a;
clear matrizv;
clear matriza;
clear matrizrotv;
clear matrizrota;
clear matrizcruv;
clear matrizcrua;
thr=8;

clear coordat;
clear csem;
clear crot;
clear ccru;

sc=1;
rc=1;
cc=1;
zc=1;

for i=100:length(Velocidad)-100
v=Velocidad(i-10:i+10);
a=Aceleracin(i-10:i+10);
X=[v a];
Y=[Velocidad(i) Aceleracin(i)];
dm=0;
if (det(X'*X)>0)
dm=mahal(Y,X);
end
%%% dm=sqrt((Velocidad(i)-mean(v))^2/var(v)+(Aceleracin(i)-mean(a))^2/var(a));
%%% if (dm>thr && Velocidad(i)<10)
if (dm>thr)
atipicos(i)=1;
coordat(cc,:)=[Latitud(i) Longitud(i)];
cc=cc+1;
for s=1:length(semaforos)
lat1=semaforos(s,1);
long1=semaforos(s,2);
lat2=Latitud(i);
long2=Longitud(i);
distsem=acos(sin(lat1*pi/180)*sin(lat2*pi/180)+cos(lat1*pi/180)*cos(lat2*pi/180)*cos(long2*pi/180-long1*pi/180))*6371000;
if (distsem<50)
v=Velocidad(i-20:i+20);
a=Aceleracin(i-20:i+20);
matrizv(:,sc)=v;
matriza(:,sc)=a;
csem(1:3,sc)=[lat2; long2; i];
sc=sc+1;
end
end

for s=1:length(rotondas)
lat1=rotondas(s,1);
long1=rotondas(s,2);
lat2=Latitud(i);
long2=Longitud(i);
distsem=acos(sin(lat1*pi/180)*sin(lat2*pi/180)+cos(lat1*pi/180)*cos(lat2*pi/180)*cos(long2*pi/180-long1*pi/180))*6371000;
if (distsem<50)
v=Velocidad(i-20:i+20);
a=Aceleracin(i-20:i+20);
matrizrotv(:,rc)=v;
matrizrota(:,rc)=a;
crot(1:3,rc)=[lat2; long2; i];
rc=rc+1;
end
end

for s=1:length(cruces)
lat1=cruces(s,1);
long1=cruces(s,2);
lat2=Latitud(i);
long2=Longitud(i);
distsem=acos(sin(lat1*pi/180)*sin(lat2*pi/180)+cos(lat1*pi/180)*cos(lat2*pi/180)*cos(long2*pi/180-long1*pi/180))*6371000;
if (distsem<50)
v=Velocidad(i-20:i+20);
a=Aceleracin(i-20:i+20);
matrizcruv(:,zc)=v;
matrizcrua(:,zc)=a;
ccru(1:3,zc)=[lat2; long2; i];
zc=zc+1;
end
end

end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% entrenando autocodificadores 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

autoencv = trainAutoencoder(matrizv);
autoenca = trainAutoencoder(matriza);
autoencrotv = trainAutoencoder(matrizrotv);
autoencrota = trainAutoencoder(matrizrota);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% recorriendo todos los trayectos 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear V;
clear A;
rec=1;

for i=100:length(Velocidad)-100
if(atipicos(i)==1)
V(:,rec)=Velocidad(i-10:i+10);
A(:,rec)=Aceleracin(i-10:i+10);
rec=rec+1;
end;
end;





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% encontrando similitudes 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

zv = encode(autoencv,V);
xv = decode(autoencv,zv);

clear simv;
clear sima;

for i=1:length(V(1,:))
simv(i)=(V(:,i)-mean(V(:,i)))'*(xv(:,i)-mean(xv(:,i)))/(norm(V(:,i)-mean(V(:,i)))*norm(xv(:,i)-mean(xv(:,i))));
end;
% figure;
% plot(sim);

za = encode(autoenca,A);
xa = decode(autoenca,za);

for i=1:length(A(1,:))
sima(i)=(A(:,i)-mean(A(:,i)))'*(xa(:,i)-mean(xa(:,i)))/(norm(A(:,i)-mean(A(:,i)))*norm(xa(:,i)-mean(xa(:,i))));
end;

clear puntos;

thr=0.95;
j=1;
for i=1:length(simv)
if (simv(i)>thr && sima(i)>thr)
%%%puntos(11*(j-1)+1:11*(j),1:2)=[Latitud(i:i+10) Longitud(i:i+10)];
puntos(j,1:2)=coordat(i,1:2);
j=j+1;
end
end

clear simv;
clear sima;

zv = encode(autoencrotv,V);
xv = decode(autoencrotv,zv);

for i=1:length(V(1,:))
simv(i)=(V(:,i)-mean(V(:,i)))'*(xv(:,i)-mean(xv(:,i)))/(norm(V(:,i)-mean(V(:,i)))*norm(xv(:,i)-mean(xv(:,i))));
end;
% figure;
% plot(sim);

za = encode(autoencrota,A);
xa = decode(autoencrota,za);

for i=1:length(A(1,:))
sima(i)=(A(:,i)-mean(A(:,i)))'*(xa(:,i)-mean(xa(:,i)))/(norm(A(:,i)-mean(A(:,i)))*norm(xa(:,i)-mean(xa(:,i))));
end;

clear puntosrot;

thr=0.95;
j=1;
for i=1:length(simv)
if (simv(i)>thr && sima(i)>thr)
%%%puntos(11*(j-1)+1:11*(j),1:2)=[Latitud(i:i+10) Longitud(i:i+10)];
puntosrot(j,1:2)=coordat(i,1:2);
j=j+1;
end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% entrenando autocodificadores con 3 capas RBM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

autoencv1 = trainAutoencoder(matrizv,20);
zv1 = encode(autoencv1,matrizv);
autoencv2 = trainAutoencoder(zv1,10);
zv2 = encode(autoencv2,zv1);
autoencv3 = trainAutoencoder(zv2,3);
zv2 = encode(autoencv3,zv2);
autoenca1 = trainAutoencoder(matriza,20);
za1 = encode(autoenca1,matriza);
autoenca2 = trainAutoencoder(za1,10);
za2 = encode(autoenca2,za1);
autoenca3 = trainAutoencoder(za2,3);
za2 = encode(autoenca3,za2);
autoencrotv1 = trainAutoencoder(matrizrotv,20);
zrotv1 = encode(autoencrotv1,matrizrotv);
autoencrotv2 = trainAutoencoder(zrotv1,10);
zrotv2 = encode(autoencrotv2,zrotv1);
autoencrotv3 = trainAutoencoder(zrotv2,3);
zrotv2 = encode(autoencrotv3,zrotv2);
autoencrota1 = trainAutoencoder(matrizrota,20);
zrota1 = encode(autoencrota1,matrizrota);
autoencrota2 = trainAutoencoder(zrota1,10);
zrota2 = encode(autoencrota2,zrota1);
autoencrota3 = trainAutoencoder(zrota2,3);
zrota2 = encode(autoencrota3,zrota2);

autoenccruv1 = trainAutoencoder(matrizcruv,20);
zcruv1 = encode(autoenccruv1,matrizcruv);
autoenccruv2 = trainAutoencoder(zcruv1,10);
zcruv2 = encode(autoenccruv2,zcruv1);
autoenccruv3 = trainAutoencoder(zcruv2,3);
zcruv2 = encode(autoenccruv3,zcruv2);
autoenccrua1 = trainAutoencoder(matrizcrua,20);
zcrua1 = encode(autoenccrua1,matrizcrua);
autoenccrua2 = trainAutoencoder(zcrua1,10);
zcrua2 = encode(autoenccrua2,zcrua1);
autoenccrua3 = trainAutoencoder(zcrua2,3);
zcrua2 = encode(autoenccrua3,zcrua2);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% recorriendo todos los trayectos 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


semaforos=[40.337563 -3.756103;40.334738 -3.757735; 40.340432 -3.756063; 40.336111 -3.758840; 40.333432 -3.760971; 40.332698 -3.761443; 40.336860 -3.757826];

rotondas=[40.342853 -3.756519; 40.347349 -3.755183; 40.334508 -3.760247; 40.335050 -3.760346; 40.322612 -3.716328; 40.327216 -3.717999];

nullclass=[40.336514 -3.728332; 40.341487 -3.727938; 40.349327 -3.737709; 40.319441 -3.715132; 40.324322 -3.714737; 40.351432 -3.752578; 40.338980 -3.755624; 40.326491 -3.718470];

cruces=[40.324288 -3.714727; 40.322070 -3.716666; 40.332173 -3.762090; 40.326514 -3.718466];

clear V;
clear A;
clear C;
rec=1;

for i=100:length(Velocidad)-100
if(atipicos(i)==1)

for s=1:length(cruces)
lat1=cruces(s,1);
long1=cruces(s,2);
lat2=Latitud(i);
long2=Longitud(i);
distsem=acos(sin(lat1*pi/180)*sin(lat2*pi/180)+cos(lat1*pi/180)*cos(lat2*pi/180)*cos(long2*pi/180-long1*pi/180))*6371000;
if (distsem<50)
V(:,rec)=Velocidad(i-20:i+20);
A(:,rec)=Aceleracin(i-20:i+20);
C(rec)=0;
rec=rec+1;
end
end

for s=1:length(semaforos)
lat1=semaforos(s,1);
long1=semaforos(s,2);
lat2=Latitud(i);
long2=Longitud(i);
distsem=acos(sin(lat1*pi/180)*sin(lat2*pi/180)+cos(lat1*pi/180)*cos(lat2*pi/180)*cos(long2*pi/180-long1*pi/180))*6371000;
if (distsem<50)
V(:,rec)=Velocidad(i-20:i+20);
A(:,rec)=Aceleracin(i-20:i+20);
C(rec)=1;
rec=rec+1;
end
end

for s=1:length(rotondas)
lat1=rotondas(s,1);
long1=rotondas(s,2);
lat2=Latitud(i);
long2=Longitud(i);
distsem=acos(sin(lat1*pi/180)*sin(lat2*pi/180)+cos(lat1*pi/180)*cos(lat2*pi/180)*cos(long2*pi/180-long1*pi/180))*6371000;
if (distsem<50)
V(:,rec)=Velocidad(i-20:i+20);
A(:,rec)=Aceleracin(i-20:i+20);
C(rec)=2;
rec=rec+1;
end
end

end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% preparando clasificador 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

zv1 = encode(autoencv1,V);
zv2 = encode(autoencv2,zv1);
zv3 = encode(autoencv3,zv2);
zrotv1 = encode(autoencrotv1,V);
zrotv2 = encode(autoencrotv2,zrotv1);
zrotv3 = encode(autoencrotv3,zrotv2);
zcruv1 = encode(autoenccruv1,V);
zcruv2 = encode(autoenccruv2,zcruv1);
zcruv3 = encode(autoenccruv3,zcruv2);
za1 = encode(autoenca1,A);
za2 = encode(autoenca2,za1);
za3 = encode(autoenca3,za2);
zrota1 = encode(autoencrota1,A);
zrota2 = encode(autoencrota2,zrota1);
zrota3 = encode(autoencrota3,zrota2);
zcrua1 = encode(autoenccrua1,A);
zcrua2 = encode(autoenccrua2,zcrua1);
zcrua3 = encode(autoenccrua3,zcrua2);
samples=[zv3' za3' zrotv3' zrota3' zcruv3' zcrua3' C']; 
%%% samples=[zv1' za1' zrotv1' zrota1' C']; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% densidad de atipicos
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


semaforos=[40.337563 -3.756103;40.334738 -3.757735; 40.340432 -3.756063; 40.336111 -3.758840; 40.333432 -3.760971; 40.332698 -3.761443; 40.336860 -3.757826];

rotondas=[40.342853 -3.756519; 40.347349 -3.755183; 40.334508 -3.760247; 40.335050 -3.760346; 40.322612 -3.716328; 40.327216 -3.717999];

nullclass=[40.336514 -3.728332; 40.341487 -3.727938; 40.349327 -3.737709; 40.319441 -3.715132; 40.324322 -3.714737; 40.351432 -3.752578; 40.338980 -3.755624; 40.326491 -3.718470];

cruces=[40.324288 -3.714727; 40.322070 -3.716666; 40.332173 -3.762090; 40.326514 -3.718466];

clear C;
clear Lat;
clear Long;
rec=1;

for i=100:length(Velocidad)-100
if(atipicos(i)==1)

for s=1:length(cruces)
lat1=cruces(s,1);
long1=cruces(s,2);
lat2=Latitud(i);
long2=Longitud(i);
distsem=acos(sin(lat1*pi/180)*sin(lat2*pi/180)+cos(lat1*pi/180)*cos(lat2*pi/180)*cos(long2*pi/180-long1*pi/180))*6371000;
if (distsem<50)
Lat(rec)=Latitud(i);
Long(rec)=Longitud(i);
C(rec)=0;
rec=rec+1;
end
end

for s=1:length(semaforos)
lat1=semaforos(s,1);
long1=semaforos(s,2);
lat2=Latitud(i);
long2=Longitud(i);
distsem=acos(sin(lat1*pi/180)*sin(lat2*pi/180)+cos(lat1*pi/180)*cos(lat2*pi/180)*cos(long2*pi/180-long1*pi/180))*6371000;
if (distsem<50)
Lat(rec)=Latitud(i);
Long(rec)=Longitud(i);
C(rec)=1;
rec=rec+1;
end
end

for s=1:length(rotondas)
lat1=rotondas(s,1);
long1=rotondas(s,2);
lat2=Latitud(i);
long2=Longitud(i);
distsem=acos(sin(lat1*pi/180)*sin(lat2*pi/180)+cos(lat1*pi/180)*cos(lat2*pi/180)*cos(long2*pi/180-long1*pi/180))*6371000;
if (distsem<50)
Lat(rec)=Latitud(i);
Long(rec)=Longitud(i);
C(rec)=2;
rec=rec+1;
end
end

end;
end;

contadores=zeros(1,length(Lat));
for i=1:length(Lat)
for j=1:length(Lat)
lat1=Lat(i);
long1=Long(i);
lat2=Lat(j);
long2=Long(j);
dist=acos(sin(lat1*pi/180)*sin(lat2*pi/180)+cos(lat1*pi/180)*cos(lat2*pi/180)*cos(long2*pi/180-long1*pi/180))*6371000;
if (dist<20)
contadores(i)=contadores(i)+1;
end
end
end

%%% samples=[Lat' Long' contadores' C']; 
samples=[samples contadores']; 












%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Casos claros
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%
% deteccion si un punto es atipico
%%%%%%%


semaforos=[40.334736 -3.757773; 40.336874 -3.757851];

rotondas=[40.343041 -3.756421; 40.347382 -3.755148];

nullclass=[40.336514 -3.728332; 40.341487 -3.727938; 40.349327 -3.737709; 40.319441 -3.715132; 40.324322 -3.714737; 40.351432 -3.752578; 40.338980 -3.755624; 40.326491 -3.718470];

cruces=[40.324288 -3.714727; 40.332182 -3.762079];

atipicos=zeros(length(Velocidad),1);
clear v;
clear a;
clear matrizv;
clear matriza;
clear matrizrotv;
clear matrizrota;
clear matrizcruv;
clear matrizcrua;
thr=5;

clear coordat;
clear csem;
clear crot;
clear ccru;

sc=1;
rc=1;
cc=1;
zc=1;

for i=100:length(Velocidad)-100
v=Velocidad(i-10:i+10);
a=Aceleracin(i-10:i+10);
X=[v a];
Y=[Velocidad(i) Aceleracin(i)];
dm=0;
if (det(X'*X)>0)
dm=mahal(Y,X);
end
%%% dm=sqrt((Velocidad(i)-mean(v))^2/var(v)+(Aceleracin(i)-mean(a))^2/var(a));
%%% if (dm>thr && Velocidad(i)<10)
if (dm>thr)
atipicos(i)=1;
coordat(cc,:)=[Latitud(i) Longitud(i)];
cc=cc+1;
for s=1:length(semaforos)
lat1=semaforos(s,1);
long1=semaforos(s,2);
lat2=Latitud(i);
long2=Longitud(i);
distsem=acos(sin(lat1*pi/180)*sin(lat2*pi/180)+cos(lat1*pi/180)*cos(lat2*pi/180)*cos(long2*pi/180-long1*pi/180))*6371000;
if (distsem<50)
v=Velocidad(i-10:i+10);
a=Aceleracin(i-10:i+10);
matrizv(:,sc)=v;
matriza(:,sc)=a;
csem(1:3,sc)=[lat2; long2; i];
sc=sc+1;
end
end

for s=1:length(rotondas)
lat1=rotondas(s,1);
long1=rotondas(s,2);
lat2=Latitud(i);
long2=Longitud(i);
distsem=acos(sin(lat1*pi/180)*sin(lat2*pi/180)+cos(lat1*pi/180)*cos(lat2*pi/180)*cos(long2*pi/180-long1*pi/180))*6371000;
if (distsem<50)
v=Velocidad(i-10:i+10);
a=Aceleracin(i-10:i+10);
matrizrotv(:,rc)=v;
matrizrota(:,rc)=a;
crot(1:3,rc)=[lat2; long2; i];
rc=rc+1;
end
end

for s=1:length(cruces)
lat1=cruces(s,1);
long1=cruces(s,2);
lat2=Latitud(i);
long2=Longitud(i);
distsem=acos(sin(lat1*pi/180)*sin(lat2*pi/180)+cos(lat1*pi/180)*cos(lat2*pi/180)*cos(long2*pi/180-long1*pi/180))*6371000;
if (distsem<50)
v=Velocidad(i-10:i+10);
a=Aceleracin(i-10:i+10);
matrizcruv(:,zc)=v;
matrizcrua(:,zc)=a;
ccru(1:3,zc)=[lat2; long2; i];
zc=zc+1;
end
end

end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% entrenando autocodificadores con 3 capas RBM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

autoencv1 = trainAutoencoder(matrizv,10);
zv1 = encode(autoencv1,matrizv);
autoencv2 = trainAutoencoder(zv1,5);
zv2 = encode(autoencv2,zv1);
autoencv3 = trainAutoencoder(zv2,1);
zv2 = encode(autoencv3,zv2);
autoenca1 = trainAutoencoder(matriza,10);
za1 = encode(autoenca1,matriza);
autoenca2 = trainAutoencoder(za1,5);
za2 = encode(autoenca2,za1);
autoenca3 = trainAutoencoder(za2,1);
za2 = encode(autoenca3,za2);
autoencrotv1 = trainAutoencoder(matrizrotv,10);
zrotv1 = encode(autoencrotv1,matrizrotv);
autoencrotv2 = trainAutoencoder(zrotv1,5);
zrotv2 = encode(autoencrotv2,zrotv1);
autoencrotv3 = trainAutoencoder(zrotv2,1);
zrotv2 = encode(autoencrotv3,zrotv2);
autoencrota1 = trainAutoencoder(matrizrota,10);
zrota1 = encode(autoencrota1,matrizrota);
autoencrota2 = trainAutoencoder(zrota1,5);
zrota2 = encode(autoencrota2,zrota1);
autoencrota3 = trainAutoencoder(zrota2,1);
zrota2 = encode(autoencrota3,zrota2);

autoenccruv1 = trainAutoencoder(matrizcruv,10);
zcruv1 = encode(autoenccruv1,matrizcruv);
autoenccruv2 = trainAutoencoder(zcruv1,5);
zcruv2 = encode(autoenccruv2,zcruv1);
autoenccruv3 = trainAutoencoder(zcruv2,1);
zcruv2 = encode(autoenccruv3,zcruv2);
autoenccrua1 = trainAutoencoder(matrizcrua,10);
zcrua1 = encode(autoenccrua1,matrizcrua);
autoenccrua2 = trainAutoencoder(zcrua1,5);
zcrua2 = encode(autoenccrua2,zcrua1);
autoenccrua3 = trainAutoencoder(zcrua2,1);
zcrua2 = encode(autoenccrua3,zcrua2);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% recorriendo todos los trayectos 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


semaforos=[40.334736 -3.757773; 40.336874 -3.757851];

rotondas=[40.343041 -3.756421; 40.347382 -3.755148];

nullclass=[40.336514 -3.728332; 40.341487 -3.727938; 40.349327 -3.737709; 40.319441 -3.715132; 40.324322 -3.714737; 40.351432 -3.752578; 40.338980 -3.755624; 40.326491 -3.718470];

cruces=[40.324288 -3.714727; 40.332182 -3.762079];

clear V;
clear A;
clear C;
rec=1;

for i=100:length(Velocidad)-100
if(atipicos(i)==1)

for s=1:length(cruces)
lat1=cruces(s,1);
long1=cruces(s,2);
lat2=Latitud(i);
long2=Longitud(i);
distsem=acos(sin(lat1*pi/180)*sin(lat2*pi/180)+cos(lat1*pi/180)*cos(lat2*pi/180)*cos(long2*pi/180-long1*pi/180))*6371000;
if (distsem<50)
V(:,rec)=Velocidad(i-10:i+10);
A(:,rec)=Aceleracin(i-10:i+10);
C(rec)=0;
rec=rec+1;
end
end

for s=1:length(semaforos)
lat1=semaforos(s,1);
long1=semaforos(s,2);
lat2=Latitud(i);
long2=Longitud(i);
distsem=acos(sin(lat1*pi/180)*sin(lat2*pi/180)+cos(lat1*pi/180)*cos(lat2*pi/180)*cos(long2*pi/180-long1*pi/180))*6371000;
if (distsem<50)
V(:,rec)=Velocidad(i-10:i+10);
A(:,rec)=Aceleracin(i-10:i+10);
C(rec)=1;
rec=rec+1;
end
end

for s=1:length(rotondas)
lat1=rotondas(s,1);
long1=rotondas(s,2);
lat2=Latitud(i);
long2=Longitud(i);
distsem=acos(sin(lat1*pi/180)*sin(lat2*pi/180)+cos(lat1*pi/180)*cos(lat2*pi/180)*cos(long2*pi/180-long1*pi/180))*6371000;
if (distsem<50)
V(:,rec)=Velocidad(i-10:i+10);
A(:,rec)=Aceleracin(i-10:i+10);
C(rec)=2;
rec=rec+1;
end
end

end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% preparando clasificador 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

zv1 = encode(autoencv1,V);
zv2 = encode(autoencv2,zv1);
zv3 = encode(autoencv3,zv2);
zrotv1 = encode(autoencrotv1,V);
zrotv2 = encode(autoencrotv2,zrotv1);
zrotv3 = encode(autoencrotv3,zrotv2);
zcruv1 = encode(autoenccruv1,V);
zcruv2 = encode(autoenccruv2,zcruv1);
zcruv3 = encode(autoenccruv3,zcruv2);
za1 = encode(autoenca1,A);
za2 = encode(autoenca2,za1);
za3 = encode(autoenca3,za2);
zrota1 = encode(autoencrota1,A);
zrota2 = encode(autoencrota2,zrota1);
zrota3 = encode(autoencrota3,zrota2);
zcrua1 = encode(autoenccrua1,A);
zcrua2 = encode(autoenccrua2,zcrua1);
zcrua3 = encode(autoenccrua3,zcrua2);
samples=[zv3' za3' zrotv3' zrota3' zcruv3' zcrua3' C']; 
%%% samples=[zv2' za2' zrotv2' zrota2' zcruv2' zcrua2' C']; 
%%% samples=[zv1' za1' zrotv1' zrota1' zcruv1' zcrua1' C']; 
%%% samples=[zv1' za1' zrotv1' zrota1' C']; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% densidad de atipicos
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


semaforos=[40.334736 -3.757773; 40.336874 -3.757851];

rotondas=[40.343041 -3.756421; 40.347382 -3.755148];

nullclass=[40.336514 -3.728332; 40.341487 -3.727938; 40.349327 -3.737709; 40.319441 -3.715132; 40.324322 -3.714737; 40.351432 -3.752578; 40.338980 -3.755624; 40.326491 -3.718470];

cruces=[40.324288 -3.714727; 40.332182 -3.762079];

clear C;
clear Lat;
clear Long;
rec=1;

for i=100:length(Velocidad)-100
if(atipicos(i)==1)

for s=1:length(cruces)
lat1=cruces(s,1);
long1=cruces(s,2);
lat2=Latitud(i);
long2=Longitud(i);
distsem=acos(sin(lat1*pi/180)*sin(lat2*pi/180)+cos(lat1*pi/180)*cos(lat2*pi/180)*cos(long2*pi/180-long1*pi/180))*6371000;
if (distsem<50)
Lat(rec)=Latitud(i);
Long(rec)=Longitud(i);
C(rec)=0;
rec=rec+1;
end
end

for s=1:length(semaforos)
lat1=semaforos(s,1);
long1=semaforos(s,2);
lat2=Latitud(i);
long2=Longitud(i);
distsem=acos(sin(lat1*pi/180)*sin(lat2*pi/180)+cos(lat1*pi/180)*cos(lat2*pi/180)*cos(long2*pi/180-long1*pi/180))*6371000;
if (distsem<50)
Lat(rec)=Latitud(i);
Long(rec)=Longitud(i);
C(rec)=1;
rec=rec+1;
end
end

for s=1:length(rotondas)
lat1=rotondas(s,1);
long1=rotondas(s,2);
lat2=Latitud(i);
long2=Longitud(i);
distsem=acos(sin(lat1*pi/180)*sin(lat2*pi/180)+cos(lat1*pi/180)*cos(lat2*pi/180)*cos(long2*pi/180-long1*pi/180))*6371000;
if (distsem<50)
Lat(rec)=Latitud(i);
Long(rec)=Longitud(i);
C(rec)=2;
rec=rec+1;
end
end

end;
end;

contadores=zeros(1,length(Lat));
for i=1:length(Lat)
for j=1:length(Lat)
lat1=Lat(i);
long1=Long(i);
lat2=Lat(j);
long2=Long(j);
dist=acos(sin(lat1*pi/180)*sin(lat2*pi/180)+cos(lat1*pi/180)*cos(lat2*pi/180)*cos(long2*pi/180-long1*pi/180))*6371000;
if (dist<20)
contadores(i)=contadores(i)+1;
end
end
end

%%% samples=[Lat' Long' contadores' C']; 
samples=[samples contadores']; 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% entrenando autocodificadores 1 capa para clasificar 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

autoencv1 = trainAutoencoder(matrizv);
autoenca1 = trainAutoencoder(matriza);
autoencrotv1 = trainAutoencoder(matrizrotv);
autoencrota1 = trainAutoencoder(matrizrota);
autoenccruv1 = trainAutoencoder(matrizcruv);
autoenccrua1 = trainAutoencoder(matrizcrua);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% recorriendo todos los trayectos 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear V;
clear A;
rec=1;

for i=100:length(Velocidad)-100
if(atipicos(i)==1)
V(:,rec)=Velocidad(i-10:i+10);
A(:,rec)=Aceleracin(i-10:i+10);
rec=rec+1;
end;
end;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% encontrando similitudes 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

zv = encode(autoencv1,V);
xv = decode(autoencv1,zv);

clear simv;
clear sima;

for i=1:length(V(1,:))
simv(i)=(V(:,i)-mean(V(:,i)))'*(xv(:,i)-mean(xv(:,i)))/(norm(V(:,i)-mean(V(:,i)))*norm(xv(:,i)-mean(xv(:,i))));
end;
% figure;
% plot(sim);

za = encode(autoenca1,A);
xa = decode(autoenca1,za);

for i=1:length(A(1,:))
sima(i)=(A(:,i)-mean(A(:,i)))'*(xa(:,i)-mean(xa(:,i)))/(norm(A(:,i)-mean(A(:,i)))*norm(xa(:,i)-mean(xa(:,i))));
end;

clear simrotv;
clear simrota;

zv = encode(autoencrotv1,V);
xv = decode(autoencrotv1,zv);

for i=1:length(V(1,:))
simrotv(i)=(V(:,i)-mean(V(:,i)))'*(xv(:,i)-mean(xv(:,i)))/(norm(V(:,i)-mean(V(:,i)))*norm(xv(:,i)-mean(xv(:,i))));
end;
% figure;
% plot(sim);

za = encode(autoencrota1,A);
xa = decode(autoencrota1,za);

for i=1:length(A(1,:))
simrota(i)=(A(:,i)-mean(A(:,i)))'*(xa(:,i)-mean(xa(:,i)))/(norm(A(:,i)-mean(A(:,i)))*norm(xa(:,i)-mean(xa(:,i))));
end;

clear simcruv;
clear simcrua;

zv = encode(autoenccruv1,V);
xv = decode(autoenccruv1,zv);

for i=1:length(V(1,:))
simcruv(i)=(V(:,i)-mean(V(:,i)))'*(xv(:,i)-mean(xv(:,i)))/(norm(V(:,i)-mean(V(:,i)))*norm(xv(:,i)-mean(xv(:,i))));
end;
% figure;
% plot(sim);

za = encode(autoenccrua1,A);
xa = decode(autoenccrua1,za);

for i=1:length(A(1,:))
simcrua(i)=(A(:,i)-mean(A(:,i)))'*(xa(:,i)-mean(xa(:,i)))/(norm(A(:,i)-mean(A(:,i)))*norm(xa(:,i)-mean(xa(:,i))));
end;

clear puntos;
clear puntosrot;
clear puntoscru;

thr=0.96;
j=1;
k=1;
l=1;
for i=1:length(simv)
minsem=min(simv(i), sima(i));
minrot=min(simrotv(i), simrota(i));
mincru=min(simcruv(i), simcrua(i));
maxim=max(max(minsem, minrot), mincru);
if (maxim>thr && maxim==minsem)
puntos(j,1:2)=coordat(i,1:2);
j=j+1;
elseif (maxim>thr && maxim==minrot)
puntosrot(k,1:2)=coordat(i,1:2);
k=k+1;
elseif (maxim>thr && maxim==mincru)
puntoscru(l,1:2)=coordat(i,1:2);
l=l+1;
end
end

