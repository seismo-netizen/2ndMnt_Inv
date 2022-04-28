clear; clc; close all;
% This is the demo code for 2nd moment inversion

%% change the path to your dirctory
addpath('/Users/qianyunlu/Documents/hmeng/2ndMnt_Italy/2ndMnt_Inv/2NDS_2018/src')

%% load stla stlo tauc for P and S waves
load Syn.mat

for i = 1:length(st)
    stla(i) = st{i}.stla;
    stlo(i) = st{i}.stlo;
    Tauc(i) = st{i}.tauc_P;
    phas(i) = 'P';
    stel(i) = 0;

    stla(i+length(st)) = st{i}.stla;
    stlo(i+length(st)) = st{i}.stlo;
    Tauc(i+length(st)) = st{i}.tauc_S;
    phas(i+length(st)) = 'S';
    stel(i+length(st)) = 0;
end

%% get partials for inversion
[G,az,toa] = getpartials_2d_generic(stla,stlo,stel,evla,evlo,evdp,...
    Vp',Vs',topl',phas,stk,dip);

d = Tauc'; 
for i=1:length(d)
    d(i)=(d(i)/2)^2;
end
[m2,~,~,vx,vy,~]=seconds_2d_cvx(G,d); % inversion using cvx tool

%% optimal solution with source parameters
m2 = m2(1:6);
X = [m2(4), m2(5); m2(5), m2(6);];
[U,Sig,V]=svd(X);
L_c=2*sqrt(max(max(Sig))); W_c=2*sqrt(Sig(2,2));
v0=m2(2:3)/m2(1); mv0=sqrt(sum(v0.^2));
v0_up = -v0(2)*sin(dip/180*pi);
v0_east = v0(1)*sin(stk/180*pi) + v0(2)*cos(dip/180*pi)*cos(stk/180*pi);
v0_north = v0(1)*cos(stk/180*pi) - v0(2)*cos(dip/180*pi)*sin(stk/180*pi);
tauc=2*sqrt(m2(1)); L0=tauc*mv0; dir=L0/L_c;

Tauc_pre = 2*sqrt(G*m2);

%% estimate upper and lower bounds of model induced uncertainties
NDF=length(d)-3;   %
misfit=sum(((G*m2(1:6)-d).^2));
sigmeas=sqrt(misfit/NDF);
chisq=sum(((G*m2(1:6)-d).^2)/sigmeas^2);
chi95=chi2inv(0.95,NDF);
MF95=chi95*sigmeas^2;

[~,L_cmax,W_cmax,vx_max,vy_max,tauc_max]=seconds_2d_cvx_maxA(G,d,MF95,1);
[~,L_cmin,W_cmin,vx_min,vy_min,tauc_min]=seconds_2d_cvx_maxA(G,d,MF95,4);
mv0_max = sqrt(vx_max^2 + vy_max^2);
mv0_min = sqrt(vx_min^2 + vy_min^2);
dir_max = mv0_max/L_cmax*tauc_max;
dir_min = mv0_min/L_cmin*tauc_min;
Vsr = (1+dir)/2*L_c/tauc;
Vsr_max = (1+dir_max)/2*L_cmax/tauc_max;
Vsr_min = (1+dir_min)/2*L_cmin/tauc_min;

%% estimate stress drop
totmom = 10^(1.5*mag + 9.1); % mag is in Mw
L_cb=L_c*1000; W_cb=W_c*1000; L_cmax=L_cmax*1000;
L_cmin=L_cmin*1000; W_cmax=W_cmax*1000; W_cmin=W_cmin*1000;

AR=W_cb/L_cb;
if(AR>=0.8)
    Abest=pi*L_cb*W_cb;
    SDBest=2.44*totmom*(1/1e6)*Abest^(-3/2);
    Amax=pi*L_cmax*W_cmax;
    SDMaxA=2.44*totmom*(1/1e6)*Amax^(-3/2);
    Amin=pi*L_cmin*W_cmin;
    SDMinA=2.44*totmom*(1/1e6)*Amin^(-3/2);

elseif(AR<=0.8)
    %ellipse formulas;   %convert stress drops to MPa
    [SDY]=getSDell(L_cb,W_cb,totmom);
    SDBest=SDY/1e6;
    [SDY]=getSDell(L_cmax,W_cmax,totmom);
    SDMaxA=SDY/1e6;
    [SDY]=getSDell(L_cmin,W_cmin,totmom);
    SDMinA=SDY/1e6;
end

%% display results
display(L_c); display(W_c); display(tauc); display(mv0); display(dir);
display(v0); display(v0_up); display(v0_east); display(v0_north);
disp(['Best Stress Drop is ',num2str(SDBest),' (MPa)'])
disp(['Upper bound Stress Drop is ',num2str(SDMinA),' (MPa)'])
disp(['Lower bound Stress Drop is ',num2str(SDMaxA),' (MPa)'])


%%
figure
subplot(2,2,1)
for i = 1:length(st)
    scatter(st{i}.stlo,st{i}.stla,100,st{i}.tauc_P,'filled'); hold on
end
plot(evlo,evla,'r+','MarkerSize',20,'linewidth',2)
colormap turbo; colorbar; box on; axis square
title('Observed \tau_c P (sec)')
ylabel('Lat')
xlabel('Lon')
set(gca,'fontsize',15)

subplot(2,2,3)
for i = 1:length(st)
    scatter(st{i}.stlo,st{i}.stla,100,st{i}.tauc_S,'filled'); hold on
end
plot(evlo,evla,'r+','MarkerSize',20,'linewidth',2)
colormap turbo; colorbar; box on; axis square
title('Observed \tau_c S (sec)')
ylabel('Lat')
xlabel('Lon')
set(gca,'fontsize',15)

subplot(2,2,2)
scatter(stlo(1:length(st)),stla(1:length(st)),100,Tauc_pre(1:length(st)),'filled'); 
hold on;
plot(evlo,evla,'r+','MarkerSize',20,'linewidth',2)
colormap turbo; colorbar; box on; axis square
title('Model \tau_c P (sec)')
ylabel('Lat')
xlabel('Lon')
set(gca,'fontsize',15)

subplot(2,2,4)
scatter(stlo(length(st)+1:end),stla(length(st)+1:end),100,Tauc_pre(length(st)+1:end),'filled');
hold on;
plot(evlo,evla,'r+','MarkerSize',20,'linewidth',2)
colormap turbo; colorbar; box on; axis square
title('Model \tau_c S (sec)')
ylabel('Lat')
xlabel('Lon')
set(gca,'fontsize',15)
