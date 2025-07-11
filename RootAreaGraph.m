
clearvars -except aalldata corrdata; clc ; close all ; 
addpath('./ClaraFunctions') ; 
addpath('./data') ; 
addpath('./mfiles') ; 
set(0,'defaultTextInterpreter','latex')
set(groot,'defaultAxesTickLabelInterpreter','latex') 
set(groot, 'defaultLegendInterpreter','latex')

load("keltyprofile.mat");

xx = transpose(kelty_At ); 
yy = kelty_z ; 

for ii = 1:length(kelty_z)
    mmn(ii) = mean(kelty_At(1:ii)) ; 
end

for jj = 1:29:length(kelty_z)
    if jj==1
        a2(1) = xx(jj) ; 
        z2(1) = yy(jj); 
    else 
        a2(end+1) = xx(jj) ; 
        z2(end+1) = yy(jj) ;
    end
end
for kk = 1:length(a2)
    mm2(kk) = mean(a2(1:kk)) ; 
end

figure(1)
plot(xx,yy, "blue", 'LineWidth', 1.5) ; hold on
plot(mmn, yy, 'red', 'LineWidth', 1.5, 'LineStyle','--'); hold on 
xlabel('Projected Area per Unit Height A (m)')
ylabel('Water Level at Vegetation hv (m)')
xlim([.1 .7])
ylim([0 3])
lgd = legend({'$A_t$','$\bar{A_v}$'}) ; 
fontsize(lgd,14,'points')

figure(2)
plot(a2,z2, 'blue', 'LineWidth', 1.5) ; hold on
plot(mm2, z2, 'red', 'LineWidth', 1.5, 'LineStyle','--'); hold on 
scatter(a2,z2, 36, 'filled', 'MarkerFaceColor', 'blue')
xlabel('Projected Area per Unit Height A (m)')
ylabel('Water Level at Vegetation hv (m)')
xlim([.1 .7])
ylim([0 3])
lgd = legend({'$A_t$','$\bar{A_v}$'}) ; 
fontsize(lgd,14,'points') 


currentfolder = pwd ; %MATLABONLINE
savfolderpath = join([currentfolder, '/ClaraFigures/RootGraphs/'], '') ; 
savfigname = "CrossCorr_CB_TpComp.png" ;
  %exportgraphics(gcf, fullfile(savfolderpath, savfigname), 'Resolution', 500) ;