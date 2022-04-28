clear; clc; close all;

%% This section is to load a source model from SRCMOD HOME and plot slip distribution
% http://equake-rc.info/SRCMOD/searchmodels/viewmodel/S/
load('s2011LORCAx01LOPE.mat');

S = s2011LORCAx01LOPE;

% mesh and slip distribution
X = S.geoX;
Y = S.geoY;
Z = S.geoZ;
slip = S.slipSPL;

% hypocenter
src_hyp_X = 0;
src_hyp_Y = 0;
src_hyp_Z = S.evDPT;

% figure
f = figure;
plot3(src_hyp_X,src_hyp_Y,src_hyp_Z,'p','MarkerFaceColor','k','MarkerSize',15)
hold on; surf(X,Y,Z,slip); grid on;
axis equal
f.CurrentAxes.ZDir = 'Reverse';
xlabel('X==EW (km)')
ylabel('Y==NS (km)')
zlabel('Depth (km)');
colorbar;
zlim([0 max(max(Z))])
set(gca,'fontsize',15)
title(S.evTAG)
view(-352+180,45)

%% This section is to generate apparent source time function (ASTF) 
Fs = 100; % sampling frequency in Hz
T_d = 3; % duration in seconds
N_t = round(T_d*Fs)+1; % total samples
t = 0:N_t-1; t = t/Fs; % time

% number of cells along-strike and along-dip
N_s = S.invNzNx(2);
N_d = S.invNzNx(1);

% cell length along-strike and along-dip in km
len_s = S.invDzDx(2); 
len_d = S.invDzDx(1); 

% slip distribution
slip = S.slipSPL;

% strike and dip angle in degrees
stk = S.srcAStke; % strike
dip = S.srcDipAn; % dip 

% mesh
x_stk = (0:N_s-1)*len_s; x_stk = x_stk - mean(x_stk);
x_dip = (0:N_d-1)*len_d; x_dip = x_dip - mean(x_dip);

% hypocenter
x_stk_hyp = 0;% hypocenter
x_dip_hyp = 0;

% rupture velocity
vr = S.srcAveVr;
vs = 3.24; % S wave velocity in the source region km/s
vp = 5.95; % P wave velocity in the source region

% Compute rupture time of each grid
for i = 1:length(x_stk)
    for j = 1:length(x_dip)
           t_offset(j,i) = sqrt((x_stk(i)-x_stk_hyp)^2 + (x_dip(j)-x_dip_hyp)^2)/vr;
    end
end

%% plot slip distribution and rupture process
figure
set(gcf,'PaperPositionMode','auto');
set(gcf,'Units','inches');
afFigurePosition = [0 0 24 13.5];
set(gcf,'Position',afFigurePosition);
set(gcf,'PaperSize',[afFigurePosition(1)*2+afFigurePosition(3) afFigurePosition(2)*2+afFigurePosition(4)]);

subplot(2,2,1)
pcolor(x_stk,x_dip,slip); hold on;
plot(x_stk_hyp,x_dip_hyp,'p','MarkerFaceColor','k','MarkerSize',20)
axis equal; axis ij;
set(gca,'fontsize',20)
title(['Fault Slip (m) stk = ' num2str(stk) '^\circ dip = ' num2str(dip) '^\circ'])
xlabel('along strike (km)')
ylabel('along dip (km)')
colorbar

subplot(2,2,3)
pcolor(x_stk,x_dip,t_offset); hold on;
plot(x_stk_hyp,x_dip_hyp,'p','MarkerFaceColor','k','MarkerSize',20)
axis equal; axis ij; colormap(jet)
set(gca,'fontsize',20)
title('Rupture Time (sec)')
xlabel('along strike (km)')
ylabel('along dip (km)')
colorbar

%% Source Time function
tr = S.srcAveTr;
STF = zeros(size(t));
for i = 1:length(x_stk)
    for j = 1:length(x_dip)
        y = STF_cell(Fs,T_d,t_offset(j,i),tr);
        STF = STF + slip(j,i)*y;
    end
end
STF = STF/max(STF);
subplot(2,2,[2 4])
plot(t,STF,'r','linewidth',3); hold on;
set(gca,'fontsize',20)
xlabel('time (sec)')


%% Apparent Source Time Function
% slowness vector [s_ex s_ey s_ez]/vs or vp in the source region
% toa: take off angle in degrees
% az: azimuth in degrees

%% read in station list
evla = 44.789; evlo = 10.747; evdp = S.evDPT; % event location
mag = S.srcMwMoS(1); 
stk = S.srcAStke; % strike
dip = S.srcDipAn; % dip 
x = [0 5 10 15 20 25 30];
Vp = [5.85 5.95 6.21 6.30 6.37 6.66 7.74]; Vp = Vp';
Vs = [3.12 3.24 3.60 3.71 3.53 3.79 4.52]; Vs = Vs';
NL = length(x); topl = x';

for i = 1:length(x)-1
    if evdp > x(i) && evdp <= x(i+1)
        Vp0 = Vp(i); Vs0 = Vs(i);
    end
end


FID = fopen('stlist.txt');
stlist = textscan(FID,'%s %s %s %s %s');
fclose(FID);

for k = 1:length(stlist{1})
    st{k}.ntwk = stlist{1}{k};
    st{k}.stnm = stlist{2}{k};
    st{k}.stla = str2double(stlist{3}{k});
    st{k}.stlo = str2double(stlist{4}{k});
    st{k}.az = azimuth(evla,evlo,st{k}.stla,st{k}.stlo);
    st{k}.dist = str2double(stlist{5}{k});

    %% generate P ASTF
    ddum = st{k}.dist; 
    save -ascii toppinputs ddum evdp NL Vp topl
    [~, ~] = system(['export DYLD_LIBRARY_PATH=""; ' './topp']);

    fid=fopen('toppoutputs');
    tt=fscanf(fid,'%g',1);
    toa=fscanf(fid,'%g',1);
    st{k}.toa_P = toa;

    toa = st{k}.toa_P; az = st{k}.az;
    s_ex = sin(toa/180*pi)*sin(az/180*pi);
    s_ey = sin(toa/180*pi)*cos(az/180*pi);
    s_ez = -cos(toa/180*pi);

    s = [s_ex s_ey s_ez]/Vp0; % P wave slowness

    ASTF = zeros(size(t));
    for i = 1:length(x_stk)
        for j = 1:length(x_dip)
            x(j,i) = (x_stk(i)-x_stk_hyp)*sin(stk/180*pi) - (x_dip(j)-x_dip_hyp)*cos(dip/180*pi)*cos(stk/180*pi);
            y(j,i) = (x_stk(i)-x_stk_hyp)*cos(stk/180*pi) + (x_dip(j)-x_dip_hyp)*cos(dip/180*pi)*sin(stk/180*pi);
            z(j,i) = (x_dip(j)-x_dip_hyp)*sin(dip/180*pi);
            t_shift(j,i) = [x(j,i) y(j,i) z(j,i)]*s';
            temp = STF_cell(Fs,T_d,t_offset(j,i)-t_shift(j,i),tr);
            ASTF = ASTF + slip(j,i)*temp;
        end
    end
    ASTF = ASTF/max(ASTF);
    st{k}.ASTF_P = ASTF;
    st{k}.tauc_P = Tau_c(ASTF,1/Fs);
    %%
    ddum = st{k}.dist; 
    save -ascii toppinputs ddum evdp NL Vs topl
    [~, ~] = system(['export DYLD_LIBRARY_PATH=""; ' './topp']);

    fid=fopen('toppoutputs');
    tt=fscanf(fid,'%g',1);
    toa=fscanf(fid,'%g',1);
    st{k}.toa_S = toa;
    
    toa = st{k}.toa_S; az = st{k}.az;
    s_ex = sin(toa/180*pi)*sin(az/180*pi);
    s_ey = sin(toa/180*pi)*cos(az/180*pi);
    s_ez = -cos(toa/180*pi);

    s = [s_ex s_ey s_ez]/Vs0; % P wave slowness

    ASTF = zeros(size(t));
    for i = 1:length(x_stk)
        for j = 1:length(x_dip)
            x(j,i) = (x_stk(i)-x_stk_hyp)*sin(stk/180*pi) - (x_dip(j)-x_dip_hyp)*cos(dip/180*pi)*cos(stk/180*pi);
            y(j,i) = (x_stk(i)-x_stk_hyp)*cos(stk/180*pi) + (x_dip(j)-x_dip_hyp)*cos(dip/180*pi)*sin(stk/180*pi);
            z(j,i) = (x_dip(j)-x_dip_hyp)*sin(dip/180*pi);
            t_shift(j,i) = [x(j,i) y(j,i) z(j,i)]*s';
            temp = STF_cell(Fs,T_d,t_offset(j,i)-t_shift(j,i),tr);
            ASTF = ASTF + slip(j,i)*temp;
        end
    end
    ASTF = ASTF/max(ASTF);
    st{k}.ASTF_S = ASTF;
    st{k}.tauc_S = Tau_c(ASTF,1/Fs);
end

figure
subplot(3,2,1)
for i = 1:length(stlist{1})
    scatter(st{i}.stlo,st{i}.stla,100,st{i}.dist,'filled'); hold on
end
plot(evlo,evla,'r+','MarkerSize',20,'linewidth',2)
colormap turbo; colorbar; box on; axis square
title('Epi. Dist. (km)')
ylabel('Lat')
set(gca,'fontsize',15)

subplot(3,2,2)
for i = 1:length(stlist{1})
    scatter(st{i}.stlo,st{i}.stla,100,st{i}.az,'filled'); hold on
end
plot(evlo,evla,'r+','MarkerSize',20,'linewidth',2)
colormap turbo; colorbar; box on; axis square
title('Azimuth (degree)')
set(gca,'fontsize',15)

subplot(3,2,3)
for i = 1:length(stlist{1})
    scatter(st{i}.stlo,st{i}.stla,100,st{i}.toa_P,'filled'); hold on
end
plot(evlo,evla,'r+','MarkerSize',20,'linewidth',2)
colormap turbo; colorbar; box on; axis square
caxis([80 120])
title('toa P (degree)')
ylabel('Lat')
set(gca,'fontsize',15)

subplot(3,2,4)
for i = 1:length(stlist{1})
    scatter(st{i}.stlo,st{i}.stla,100,st{i}.toa_S,'filled'); hold on
end
plot(evlo,evla,'r+','MarkerSize',20,'linewidth',2)
colormap turbo; colorbar; box on; axis square
caxis([80 120])
title('toa S (degree)')
set(gca,'fontsize',15)


subplot(3,2,5)
for i = 1:length(stlist{1})
    scatter(st{i}.stlo,st{i}.stla,100,st{i}.tauc_P,'filled'); hold on
end
plot(evlo,evla,'r+','MarkerSize',20,'linewidth',2)
colormap turbo; colorbar; box on; axis square
title('\tau_c P (sec)')
ylabel('Lat')
xlabel('Lon')
set(gca,'fontsize',15)

subplot(3,2,6)
for i = 1:length(stlist{1})
    scatter(st{i}.stlo,st{i}.stla,100,st{i}.tauc_S,'filled'); hold on
end
plot(evlo,evla,'r+','MarkerSize',20,'linewidth',2)
colormap turbo; colorbar; box on; axis square
title('\tau_c S (sec)')
xlabel('Lon')
set(gca,'fontsize',15)


save('Syn.mat','evla','evlo','evdp','mag','stk','dip','x','Vp','Vs','st','topl')

%% uncomment this section to generate ASTF from any direction in the source region
% toa = 85; az = 120;
% s_ex = sin(toa/180*pi)*sin(az/180*pi);
% s_ey = sin(toa/180*pi)*cos(az/180*pi);
% s_ez = -cos(toa/180*pi);
% 
% s = [s_ex s_ey s_ez]/vs; % S wave slowness
% %s = [s_ex s_ey s_ez]/vp; % P wave slowness
%
% ASTF = zeros(size(t));
% for i = 1:length(x_stk)
%     for j = 1:length(x_dip)
%         x(j,i) = (x_stk(i)-x_stk_hyp)*sin(stk/180*pi) - (x_dip(j)-x_dip_hyp)*cos(dip/180*pi)*cos(stk/180*pi);
%         y(j,i) = (x_stk(i)-x_stk_hyp)*cos(stk/180*pi) + (x_dip(j)-x_dip_hyp)*cos(dip/180*pi)*sin(stk/180*pi);
%         z(j,i) = (x_dip(j)-x_dip_hyp)*sin(dip/180*pi);
%         t_shift(j,i) = [x(j,i) y(j,i) z(j,i)]*s';
%         temp = STF_cell(Fs,T_d,t_offset(j,i)-t_shift(j,i),tr);
%         ASTF = ASTF + slip(j,i)*temp;
%     end
% end
% ASTF = ASTF/max(ASTF);

%% 

% for k = 1:18
%     
% 
%     toa = 90; az = (k-1)*20;
%     s_ex = sin(toa/180*pi)*sin(az/180*pi);
%     s_ey = sin(toa/180*pi)*cos(az/180*pi);
%     s_ez = -cos(toa/180*pi);
%     
%     % S wave
%     s = [s_ex s_ey s_ez]/vs;
%     
%     % fault plane location in x y z
%     ASTF = zeros(size(t));
%     for i = 1:length(x_stk)
%         for j = 1:length(x_dip)
%             x(j,i) = (x_stk(i)-x_stk_hyp)*sin(stk/180*pi) - (x_dip(j)-x_dip_hyp)*cos(dip/180*pi)*cos(stk/180*pi);
%             y(j,i) = (x_stk(i)-x_stk_hyp)*cos(stk/180*pi) + (x_dip(j)-x_dip_hyp)*cos(dip/180*pi)*sin(stk/180*pi);
%             z(j,i) = (x_dip(j)-x_dip_hyp)*sin(dip/180*pi);
%             t_shift(j,i) = [x(j,i) y(j,i) z(j,i)]*s';
%             temp = STF_cell(Fs,T_d,t_offset(j,i)-t_shift(j,i),tr);
%             ASTF = ASTF + slip(j,i)*temp;
%         end
%     end
%     ASTF = ASTF/max(ASTF);
%     plot(t,ASTF-k,'k','linewidth',3); hold on;
% 
%     
%     % P wave
%     s = [s_ex s_ey s_ez]/vp;
%     
%     % fault plane location in x y z
%     ASTF = zeros(size(t));
%     for i = 1:length(x_stk)
%         for j = 1:length(x_dip)
%             x(j,i) = (x_stk(i)-x_stk_hyp)*sin(stk/180*pi) - (x_dip(j)-x_dip_hyp)*cos(dip/180*pi)*cos(stk/180*pi);
%             y(j,i) = (x_stk(i)-x_stk_hyp)*cos(stk/180*pi) + (x_dip(j)-x_dip_hyp)*cos(dip/180*pi)*sin(stk/180*pi);
%             z(j,i) = (x_dip(j)-x_dip_hyp)*sin(dip/180*pi);
%             t_shift(j,i) = [x(j,i) y(j,i) z(j,i)]*s';
%             temp = STF_cell(Fs,T_d,t_offset(j,i)-t_shift(j,i),tr);
%             ASTF = ASTF + slip(j,i)*temp;
%         end
%     end
%     ASTF = ASTF/max(ASTF);
%     plot(t,ASTF-k,'color',[0.5 0.5 0.5],'linewidth',3); hold on;
% 
% end
% 
% legend('STF','S wave ASTF','P wave ASTF')
% 
% y_tick = -18:1:0;
% yticks(y_tick)
% yticklabels({'340^\circ','320^\circ','300^\circ','280^\circ','260^\circ','240^\circ','220^\circ','200^\circ','180^\circ',...
%     '160^\circ','140^\circ','120^\circ','100^\circ','80^\circ','60^\circ','40^\circ','20^\circ','0^\circ','STF'})
% title(['Apparent Source Time Function; toa = ' num2str(toa) '^\circ'])


%% wavelet 0 to tr intergrate to unity
function y = STF_cell(Fs,T_d,t_offset,tr)
% Fs: sampling rate in Hz
% T_d: total duration in sec
% T_offset: ruptured time in sec
% tr: rise time in sec
% The wavelet is a parabolic function integrate to unity within the tr
% window

N_t = round(T_d*Fs)+1; % total samples
t = 0:N_t-1; t = t/Fs; % time
y = -6/tr^3*(t-t_offset).*(t-t_offset-tr);
for i = 1:length(t)
    if t(i) < t_offset
        y(i) = 0;
    elseif t(i) > t_offset+tr
        y(i) = 0;
    end
end

end