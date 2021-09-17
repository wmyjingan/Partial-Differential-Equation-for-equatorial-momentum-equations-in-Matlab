clear
clc

load topo_025deg
[lon01,lat01]=meshgrid(lon,lat);

hh=30;

load uv_remove_14_40day_600_3kkm_smooth_std
ub=u;
vb=v;

poslon=lon>=175&lon<=275;
poslat=abs(lat-0)<=11;

lon1=lon(poslon);
lat1=lat(poslat);
ub=ub(:,poslat,poslon);
vb=vb(:,poslat,poslon);
tim1=tim;
date1=datestr(tim1);

clear u up v vp

load ssh_14_40day_600_3kkm_fft2 ssh lon lat tim

ssh(ssh==0)=nan;

sst=zeros(size(ssh));
sss=zeros(size(sst));

poslon=lon>=175&lon<=275;
poslat=abs(lat-0)<=11;
lon2=lon(poslon);
lat2=lat(poslat);
ssh=ssh(:,poslat,poslon);
lon=lon2;
lat=lat2;

gg=9.8;
fai=gg*(ssh);
ssh1=squeeze(fai(200,:,:));
tim2=tim;
date2=datestr(tim2);

%
timup=max(tim1(1),tim2(1));
timend=min(tim1(end),tim2(end));
postim1=tim1>=timup&tim1<=timend;
postim2=tim2>=timup&tim2<=timend;
ub=ub(postim1,:,:);
vb=vb(postim1,:,:);
fai=fai(postim2,:,:);
tim=tim1(postim1);
date=datestr(tim);
tim0=tim;
lon0=lon;
lat0=lat;

sst=sst(postim2,poslat,poslon);
sss=sss(postim2,poslat,poslon);

% fai, ub, lon, lat, tim
%%
nx=length(lon);
ny=length(lat);
nt=length(tim);
a=6371393;
deg=2*pi*a/360;
omega=2*pi/24/3600;
bb=2*omega/a;
yy=lat*deg;
ff=bb*yy;
% lon=lon';
% lat=lat';
lonu=[lon(1)-0.5*diff(lon(1:2)),lon(1:end-1)+0.5*diff(lon),...
    lon(end)+0.5*diff(lon(end-1:end))];
latv=[lat(1)-0.5*diff(lat(1:2)),lat(1:end-1)+0.5*diff(lat),...
    lat(end)+0.5*diff(lat(end-1:end))];
[x0,y0]=meshgrid(lon,lat);
[xu,yu]=meshgrid(lonu,lat);
[xv,yv]=meshgrid(lon,latv);

cc=2.8;
rr=sqrt(cc/2/bb);
aa=bb*rr;
dt=diff(tim(1:2))*24*3600;
dx=diff(lon(1:2))*deg;
dy=diff(lat(1:2))*deg;

nn_b=ny*(nx+1)+nx*(ny+1);
nn=(ny-2)*(nx-1)+(nx-2)*(ny-1);

% load tauxy_14_40day_600_3kkm_fft2
% postim=tim>=timup&tim<=timend;
% poslon=lon>=170&lon<=270; 
% poslat=abs(lat-0)<=11;
taux=zeros(size(ssh));
tauy=zeros(size(ssh));
rou0=1024;
hh=30;

%% boundary and initial conditions
% load uvg_14_40day_600_3kkm_fft2
% u(isnan(u))=0;
% v(isnan(v))=0;
% 
% [xr,yr]=meshgrid(lon,lat);
% postim=tim>=timup&tim<=timend;
% u=u(postim,:,:);
% v=v(postim,:,:);
% data=squeeze(v(45,abs(lat)<=6,lon>=175&lon<=265));

ubdn=single(zeros(nt,ny,nx+1));
vbdn=single(zeros(nt,ny+1,nx));
% for t=1:nt
%     data=squeeze(u(t,:,:));
%     ubdn(t,:,:)=interp2(xr,yr,data,xu,yu,'linear');
%     data=squeeze(v(t,:,:));
%     vbdn(t,:,:)=interp2(xr,yr,data,xv,yv,'linear');
% end   
ubdn(:,:,:)=0;
vbdn(:,:,:)=0;

ul=ubdn(:,:,1);
ur=ubdn(:,:,end);
uu=ubdn(:,end,:);
ud=ubdn(:,1,:);
vu=vbdn(:,end,:);
vd=vbdn(:,1,:);
vl=vbdn(:,:,1);
vr=vbdn(:,:,end);

u00=squeeze(ubdn(1,:,:));
v00=squeeze(vbdn(1,:,:));

u1=u00(2:end-1,2:end-1);
v1=v00(2:end-1,2:end-1);
U1=[reshape(u1',1,(ny-2)*(nx-1)),reshape(v1',1,(nx-2)*(ny-1))];

u0=zeros(nt,ny,nx+1);
v0=zeros(nt,ny+1,nx);

u0(1,:,:)=u00;
u0(:,:,1)=ul;
u0(:,:,end)=ur;
u0(:,1,:)=ud;
u0(:,end,:)=uu;
v0(1,:,:)=v00;
v0(:,1,:)=vd;
v0(:,end,:)=vu;
v0(:,:,1)=vl;
v0(:,:,end)=vr;
U0=zeros(nn_b,nt);
for t=1:nt
U0(:,t)=[reshape(squeeze(u0(t,:,:))',ny*(nx+1),1);...
    reshape(squeeze(v0(t,:,:))',nx*(ny+1),1)];
end

u=u0;
v=v0;

vb(:,:,:)=0;

lon=lon0;
lat=lat0;
tim=tim0;

%%

% ah=1e3;
ah=0;
aat=3e-4*1.;
bbs=7.4e-4*1.;

bsst=gg*hh/2*aat*sst;
bsss=-gg*hh/2*bbs*sss;

jj=5;
ii=10;
au=1/ii/24/3600;
av=1/jj/24/3600;
aau=1/dt+au/2;
ccu=1/dt-au/2;
aav=1/dt+av/2;
ccv=1/dt-av/2;

clear sst sss ssh 

topou=interp2(lon01,lat01,topo,xu,yu,'linear');
topov=interp2(lon01,lat01,topo,xv,yv,'linear');

taux(:,:,:)=0;
tauy(:,:,:)=0;

for t=1:nt-1
%     t=1;
fai1=squeeze(nanmean(fai(t:t+1,:,:),1));
bsst1=squeeze(nanmean(bsst(t:t+1,:,:),1));
bsss1=squeeze(nanmean(bsss(t:t+1,:,:),1));
fai1=fai1-bsst1-bsss1;

taux1=squeeze(nanmean(taux(t:t+1,:,:),1))/rou0/hh;
tauy1=squeeze(nanmean(tauy(t:t+1,:,:),1))/rou0/hh;

ub1=squeeze(nanmean(ub(t:t+1,:,:),1));
vb1=squeeze(nanmean(vb(t:t+1,:,:),1));
% for j=1:ny
%     ub1(j,:)=smooth(ub1(j,:),3);
%     vb1(j,:)=smooth(vb1(j,:),5);
% end

% fai1(isnan(topo))=nan;
% ub1(isnan(topo))=nan;
% vb1(isnan(topo))=nan;
% taux1(isnan(topo))=nan;
% tauy1(isnan(topo))=nan;

% fai1(isnan(fai1))=0;
% ub1(isnan(ub1))=0;
% vb1(isnan(vb1))=0;
% taux1(isnan(taux1))=0;
% tauy1(isnan(tauy1))=0;

%% F_b
A_b=single(zeros(nn,nn_b));
% u-momentum
for j=1:ny-2
for i=1:nx-1
A11=zeros(ny,nx+1);
A12=zeros(ny+1,nx);

ub11=0.5*(ub1(j+1,i)+ub1(j+1,i+1));
vb11=0.5*(vb1(j+1,i)+vb1(j+1,i+1));
uy1=(ub1(j+2,i)-ub1(j,i)+ub1(j+2,i+1)-ub1(j,i+1))/4/dy;
ux1=(ub1(j+1,i+1)-ub1(j+1,i))/dx;
ux1=0;

uy1=1.0*uy1;

ub11(isnan(ub11))=0;
vb11(isnan(vb11))=0;
uy1(isnan(uy1))=0;
ux1(isnan(ux1))=0;

A11(j,i+1)=-vb11/4/dy-ah/2/dy/dy;
A11(j+1,i)=-ub11/4/dx-ah/2/dx/dx;
A11(j+1,i+1)=aau+ux1/2+ah/dx/dx+ah/dy/dy;
A11(j+1,i+2)=ub11/4/dx-ah/2/dx/dx;
A11(j+2,i+1)=vb11/4/dy-ah/2/dy/dy;

A12(j+1,i)=-(bb*lat(j+1)*deg-uy1)/8;
A12(j+1,i+1)=-(bb*lat(j+1)*deg-uy1)/8;
A12(j+2,i)=-(bb*lat(j+1)*deg-uy1)/8;
A12(j+2,i+1)=-(bb*lat(j+1)*deg-uy1)/8;

A_b((j-1)*(nx-1)+i,:)=[reshape(A11',1,ny*(nx+1)),reshape(A12',1,nx*(ny+1))];
end
end
% v-momentum
for j=1:ny-1
for i=1:nx-2
A21=zeros(ny,nx+1);
A22=zeros(ny+1,nx);

ub11=0.5*(ub1(j,i+1)+ub1(j+1,i+1));
vb11=0.5*(vb1(j,i+1)+vb1(j+1,i+1));
vy1=(vb1(j+1,i+1)-vb1(j,i+1))/dy;
vx1=(vb1(j,i+2)-vb1(j,i)+vb1(j+1,i+2)-vb1(j+1,i))/4/dx;

ub11(isnan(ub11))=0;
vb11(isnan(vb11))=0;
vy1(isnan(vy1))=0;
vx1(isnan(vx1))=0;

A21(j,i+1)=(vx1+bb*latv(j+1)*deg)/8;
A21(j,i+2)=(vx1+bb*latv(j+1)*deg)/8;
A21(j+1,i+1)=(vx1+bb*latv(j+1)*deg)/8;
A21(j+1,i+2)=(vx1+bb*latv(j+1)*deg)/8;

A22(j,i+1)=-vb11/4/dy-ah/2/dy/dy;
A22(j+1,i)=-ub11/4/dx-ah/2/dx/dx;
A22(j+1,i+1)=aav+vy1/2+ah/dx/dx+ah/dy/dy;
A22(j+1,i+2)=ub11/4/dx-ah/2/dx/dx;
A22(j+2,i+1)=vb11/4/dy-ah/2/dy/dy;

A_b((ny-2)*(nx-1)+(nx-2)*(j-1)+i,:)=[reshape(A21',1,ny*(nx+1)),reshape(A22',1,nx*(ny+1))];

end
end

B_b=single(zeros(nn,nn_b));
% u-momentum
for j=1:ny-2
for i=1:nx-1
A11=zeros(ny,nx+1);
A12=zeros(ny+1,nx);

ub11=0.5*(ub1(j+1,i)+ub1(j+1,i+1));
vb11=0.5*(vb1(j+1,i)+vb1(j+1,i+1));
uy1=(ub1(j+2,i)-ub1(j,i)+ub1(j+2,i+1)-ub1(j,i+1))/4/dy;
ux1=(ub1(j+1,i+1)-ub1(j+1,i))/dx;
ux1=0;

uy1=1.0*uy1;

ub11(isnan(ub11))=0;
vb11(isnan(vb11))=0;
uy1(isnan(uy1))=0;
ux1(isnan(ux1))=0;

A11(j,i+1)=vb11/4/dy+ah/2/dy/dy;
A11(j+1,i)=ub11/4/dx+ah/2/dx/dx;
A11(j+1,i+1)=ccu-ux1/2-ah/dx/dx-ah/dy/dy;
A11(j+1,i+2)=-ub11/4/dx+ah/2/dx/dx;
A11(j+2,i+1)=-vb11/4/dy+ah/2/dy/dy;

A12(j+1,i)=(bb*lat(j+1)*deg-uy1)/8;
A12(j+1,i+1)=(bb*lat(j+1)*deg-uy1)/8;
A12(j+2,i)=(bb*lat(j+1)*deg-uy1)/8;
A12(j+2,i+1)=(bb*lat(j+1)*deg-uy1)/8;

B_b((j-1)*(nx-1)+i,:)=[reshape(A11',1,ny*(nx+1)),reshape(A12',1,nx*(ny+1))];
end
end
% v-momentum
for j=1:ny-1
for i=1:nx-2
A21=zeros(ny,nx+1);
A22=zeros(ny+1,nx);

ub11=0.5*(ub1(j,i+1)+ub1(j+1,i+1));
vb11=0.5*(vb1(j,i+1)+vb1(j+1,i+1));
vy1=(vb1(j+1,i+1)-vb1(j,i+1))/dy;
vx1=(vb1(j,i+2)-vb1(j,i)+vb1(j+1,i+2)-vb1(j+1,i))/4/dx;

ub11(isnan(ub11))=0;
vb11(isnan(vb11))=0;
vy1(isnan(vy1))=0;
vx1(isnan(vx1))=0;

A21(j,i+1)=-(vx1+bb*latv(j+1)*deg)/8;
A21(j,i+2)=-(vx1+bb*latv(j+1)*deg)/8;
A21(j+1,i+1)=-(vx1+bb*latv(j+1)*deg)/8;
A21(j+1,i+2)=-(vx1+bb*latv(j+1)*deg)/8;

A22(j,i+1)=vb11/4/dy+ah/2/dy/dy;
A22(j+1,i)=ub11/4/dx+ah/2/dx/dx;
A22(j+1,i+1)=ccv-vy1/2-ah/dx/dx-ah/dy/dy;
A22(j+1,i+2)=-ub11/4/dx+ah/2/dx/dx;
A22(j+2,i+1)=-vb11/4/dy+ah/2/dy/dy;

B_b((ny-2)*(nx-1)+(nx-2)*(j-1)+i,:)=[reshape(A21',1,ny*(nx+1)),reshape(A22',1,nx*(ny+1))];

end
end

F_b=zeros(nn,1);
F_b(:,1)=B_b*U0(:,t)-A_b*U0(:,t+1);

%% F
F=single(zeros(nn,1));
% u-momentum
for j=1:ny-2
for i=1:nx-1
F((j-1)*(nx-1)+i,1)=-(fai1(j+1,i+1)-fai1(j+1,i))/dx+...
    0.5*(taux1(j+1,i+1)+taux1(j+1,i));
end
end
% v-momentum
for j=1:ny-1
for i=1:nx-2
F((ny-2)*(nx-1)+(nx-2)*(j-1)+i,1)=-(fai1(j+1,i+1)-fai1(j,i+1))/dy+...
    0.5*(tauy1(j+1,i+1)+tauy1(j,i+1));
end
end
F(isnan(F))=0;

%% NLN
F_n=single(zeros(nn,1));
ut1=squeeze(u(t,:,:));
vt1=squeeze(v(t,:,:));

% u-momentum
for j=1:ny-2
for i=1:nx-1    
uux=ut1(j+1,i+1)*(ut1(j+1,i+2)-ut1(j+1,i))/2/dx;
vuy=(vt1(j+1,i)+vt1(j+1,i+1)+vt1(j+2,i)+vt1(j+2,i+1))/4 ...
    *(ut1(j+2,i+1)-ut1(j,i+1))/2/dy;
F_n((j-1)*(nx-1)+i,1)=-(uux+vuy);
end
end
% v-momentum
for j=1:ny-1
for i=1:nx-2
uvx=(ut1(j,i+1)+ut1(j,i+2)+ut1(j+1,i+1)+ut1(j+1,i+2))/4 ...
    *(vt1(j+1,i+2)-vt1(j+1,i))/2/dx;
vvy=vt1(j+1,i+1)*(vt1(j+2,i+1)-vt1(j,i+1))/2/dy;
F_n((ny-2)*(nx-1)+(nx-2)*(j-1)+i,1)=-(uvx+vvy);
end
end
F_n(isnan(F_n))=0;

%% A
A=single(zeros(nn,nn));
% u-momentum
for j=1:ny-2
for i=1:nx-1
A11=zeros(ny,nx+1);
A12=zeros(ny+1,nx);

ub11=0.5*(ub1(j+1,i)+ub1(j+1,i+1));
vb11=0.5*(vb1(j+1,i)+vb1(j+1,i+1));
uy1=(ub1(j+2,i)-ub1(j,i)+ub1(j+2,i+1)-ub1(j,i+1))/4/dy;
ux1=(ub1(j+1,i+1)-ub1(j+1,i))/dx;
ux1=0;

uy1=1.0*uy1;

ub11(isnan(ub11))=0;
vb11(isnan(vb11))=0;
uy1(isnan(uy1))=0;
ux1(isnan(ux1))=0;

A11(j,i+1)=-vb11/4/dy-ah/2/dy/dy;
A11(j+1,i)=-ub11/4/dx-ah/2/dx/dx;
A11(j+1,i+1)=aau+ux1/2+ah/dx/dx+ah/dy/dy;
A11(j+1,i+2)=ub11/4/dx-ah/2/dx/dx;
A11(j+2,i+1)=vb11/4/dy-ah/2/dy/dy;

A12(j+1,i)=-(bb*lat(j+1)*deg-uy1)/8;
A12(j+1,i+1)=-(bb*lat(j+1)*deg-uy1)/8;
A12(j+2,i)=-(bb*lat(j+1)*deg-uy1)/8;
A12(j+2,i+1)=-(bb*lat(j+1)*deg-uy1)/8;

A11(:,1)=[];A11(:,end)=[];
A11(1,:)=[];A11(end,:)=[];
A12(1,:)=[];A12(end,:)=[];
A12(:,1)=[];A12(:,end)=[];

A((j-1)*(nx-1)+i,:)=[reshape(A11',1,(ny-2)*(nx-1)),reshape(A12',1,(nx-2)*(ny-1))];
end
end
% v-momentum
for j=1:ny-1
for i=1:nx-2
A21=zeros(ny,nx+1);
A22=zeros(ny+1,nx);

ub11=0.5*(ub1(j,i+1)+ub1(j+1,i+1));
vb11=0.5*(vb1(j,i+1)+vb1(j+1,i+1));
vy1=(vb1(j+1,i+1)-vb1(j,i+1))/dy;
vx1=(vb1(j,i+2)-vb1(j,i)+vb1(j+1,i+2)-vb1(j+1,i))/4/dx;

ub11(isnan(ub11))=0;
vb11(isnan(vb11))=0;
vy1(isnan(vy1))=0;
vx1(isnan(vx1))=0;

A21(j,i+1)=(vx1+bb*latv(j+1)*deg)/8;
A21(j,i+2)=(vx1+bb*latv(j+1)*deg)/8;
A21(j+1,i+1)=(vx1+bb*latv(j+1)*deg)/8;
A21(j+1,i+2)=(vx1+bb*latv(j+1)*deg)/8;

A22(j,i+1)=-vb11/4/dy-ah/2/dy/dy;
A22(j+1,i)=-ub11/4/dx-ah/2/dx/dx;
A22(j+1,i+1)=aav+vy1/2+ah/dx/dx+ah/dy/dy;
A22(j+1,i+2)=ub11/4/dx-ah/2/dx/dx;
A22(j+2,i+1)=vb11/4/dy-ah/2/dy/dy;

A21(:,1)=[];A21(:,end)=[];
A21(1,:)=[];A21(end,:)=[];
A22(1,:)=[];A22(end,:)=[];
A22(:,1)=[];A22(:,end)=[];

A((ny-2)*(nx-1)+(nx-2)*(j-1)+i,:)=[reshape(A21',1,(ny-2)*(nx-1)),reshape(A22',1,(nx-2)*(ny-1))];

end
end


%% B

B=single(zeros(nn,nn));
% u-momentum
for j=1:ny-2
for i=1:nx-1
A11=zeros(ny,nx+1);
A12=zeros(ny+1,nx);

ub11=0.5*(ub1(j+1,i)+ub1(j+1,i+1));
vb11=0.5*(vb1(j+1,i)+vb1(j+1,i+1));
uy1=(ub1(j+2,i)-ub1(j,i)+ub1(j+2,i+1)-ub1(j,i+1))/4/dy;
ux1=(ub1(j+1,i+1)-ub1(j+1,i))/dx;
ux1=0;

uy1=1.0*uy1;

ub11(isnan(ub11))=0;
vb11(isnan(vb11))=0;
uy1(isnan(uy1))=0;
ux1(isnan(ux1))=0;

A11(j,i+1)=vb11/4/dy+ah/2/dy/dy;
A11(j+1,i)=ub11/4/dx+ah/2/dx/dx;
A11(j+1,i+1)=ccu-ux1/2-ah/dx/dx-ah/dy/dy;
A11(j+1,i+2)=-ub11/4/dx+ah/2/dx/dx;
A11(j+2,i+1)=-vb11/4/dy+ah/2/dy/dy;

A12(j+1,i)=(bb*lat(j+1)*deg-uy1)/8;
A12(j+1,i+1)=(bb*lat(j+1)*deg-uy1)/8;
A12(j+2,i)=(bb*lat(j+1)*deg-uy1)/8;
A12(j+2,i+1)=(bb*lat(j+1)*deg-uy1)/8;

A11(:,1)=[];A11(:,end)=[];
A11(1,:)=[];A11(end,:)=[];
A12(1,:)=[];A12(end,:)=[];
A12(:,1)=[];A12(:,end)=[];

B((j-1)*(nx-1)+i,:)=[reshape(A11',1,(ny-2)*(nx-1)),reshape(A12',1,(nx-2)*(ny-1))];
end
end
% v-momentum
for j=1:ny-1
for i=1:nx-2
A21=zeros(ny,nx+1);
A22=zeros(ny+1,nx);

ub11=0.5*(ub1(j,i+1)+ub1(j+1,i+1));
vb11=0.5*(vb1(j,i+1)+vb1(j+1,i+1));
vy1=(vb1(j+1,i+1)-vb1(j,i+1))/dy;
vx1=(vb1(j,i+2)-vb1(j,i)+vb1(j+1,i+2)-vb1(j+1,i))/4/dx;

ub11(isnan(ub11))=0;
vb11(isnan(vb11))=0;
vy1(isnan(vy1))=0;
vx1(isnan(vx1))=0;

A21(j,i+1)=-(vx1+bb*latv(j+1)*deg)/8;
A21(j,i+2)=-(vx1+bb*latv(j+1)*deg)/8;
A21(j+1,i+1)=-(vx1+bb*latv(j+1)*deg)/8;
A21(j+1,i+2)=-(vx1+bb*latv(j+1)*deg)/8;

A22(j,i+1)=vb11/4/dy+ah/2/dy/dy;
A22(j+1,i)=ub11/4/dx+ah/2/dx/dx;
A22(j+1,i+1)=ccv-vy1/2-ah/dx/dx-ah/dy/dy;
A22(j+1,i+2)=-ub11/4/dx+ah/2/dx/dx;
A22(j+2,i+1)=-vb11/4/dy+ah/2/dy/dy;

A21(:,1)=[];A21(:,end)=[];
A21(1,:)=[];A21(end,:)=[];
A22(1,:)=[];A22(end,:)=[];
A22(:,1)=[];A22(:,end)=[];

B((ny-2)*(nx-1)+(nx-2)*(j-1)+i,:)=[reshape(A21',1,(ny-2)*(nx-1)),reshape(A22',1,(nx-2)*(ny-1))];

end
end


%% solve
A1=A;
B1=B;
F1=F;
F_b1=F_b;
U11=U1;
%%%%% island
% Arow=mean(A,1);
% ocean=~isnan(Arow);
% iland=isnan(Arow);
% nnan=length(Arow(iland));
% 
% A1=A;
% A1(:,iland)=[];
% A1(iland,:)=[];
% 
% B1=B;
% B1(:,iland)=[];
% B1(iland,:)=[];
% 
% F1=F;
% F1(iland)=[];
% 
% F_b1=F_b;
% F_b1(iland)=[];
% 
% U11=U1;
% U11(iland)=[];
%%%%%%
A1(isnan(A1))=0;
U11(isnan(U11))=0;
F1(isnan(F1))=0;
F_b1(isnan(F_b1))=0;

lhs=A1;
rhs=B1*U11'+F1+F_b1;

solv=lhs\rhs;
%     solv(solv>=2)=2;
%     solv(solv<=-2)=-2;
% solv1=nan(nn,1);
% solv1(ocean)=solv;
solv1=solv;

    u2=reshape(solv1(1:(ny-2)*(nx-1)),nx-1,ny-2)';
    v2=reshape(solv1(1+(ny-2)*(nx-1):end),nx-2,ny-1)';
    
    u2(isnan(topou(2:end-1,2:end-1)))=nan;
    v2(isnan(topov(2:end-1,2:end-1)))=nan;
    
U2=[reshape(u2',1,(ny-2)*(nx-1)),reshape(v2',1,(ny-1)*(nx-2))];

u(t+1,2:end-1,2:end-1)=u2;
v(t+1,2:end-1,2:end-1)=v2;

% eke(t)=0.5*(nanmean(nanmean(u2.^2))+nanmean(nanmean(v2.^2)))*1e4;
% figure(1)
% plot(eke)

U1=U2;

t
end


um=u;
vm=v;


u=um(1:end,:,:);
v=vm(1:end,:,:);
u(u==0)=nan;
v(v==0)=nan;
tim=tim(1:end);
save('uv_c_U_aa_5_10_14_40day_600_3kkm_1993_2020_oscar_aviso.mat',...
    'u','v', 'lon', 'lonu','lat','latv','tim','-v7.3')

figure
data=squeeze(nanstd(v(:,:,:),1,1));
contourf(lon,latv,data,20,'linestyle','none')
colorbar
caxis([0 0.3])
figure
data=squeeze(nanstd(u(:,:,:),1,1));
contourf(lonu,lat,data,20,'linestyle','none')
colorbar
caxis([0 0.3])



