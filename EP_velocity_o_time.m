%  clear ; clc ; close all ; 
addpath('./ClaraFunctions') ; 
addpath('./data') ; 
addpath('./mfiles') ; 
     load('aalldata_Mar282025.mat') ;
%  savfigname = 'NewADVVelocityPlotsAllLayoutsWide.png' ; 
savfolderpath = '/home/elizabeth/Desktop/cshorex-main/osu_mangrove/ClaraFigures/Velocity_over_Time_20250414/'; 
%        saveas(gcf, fullfile(savfolderpath, savfigname)) 
% set(0, 'DefaultAxesFontName', 'Nimbus Roman', 'DefaultTextFontName', 'Nimbus Roman')
 set(0,'defaultTextInterpreter','latex')
 set(groot,'defaultAxesTickLabelInterpreter','latex') 
 set(groot, 'defaultLegendInterpreter','latex')
 cnt = 1;
 clear fieldnames ; fieldnames = fieldnames(aalldata) ;
 ADVlabel = {'ADV 2, $d=1.404$m', 'ADV 3, $d=1.550$m','ADV 4, $d=1.720$m','ADV 5, $d=1.858$m'} ; 
zu = [1.4040, 1.5500, 1.7200, 1.8580] ; % ADV 2,3,4,5 (order is correct and checked)
for totalnum = 1:length(fieldnames)
    categoryname = fieldnames{totalnum} ; % 'HighDensity_h270_hv182_NoWall' ; %
    set_category_variables
%% Length Dep Constants
NumofTrials = length(t) ; 
%% Start of Trial Indp Sections
for num = [1,NumofTrials]
    if num == NumofTrials && NumofTrials <10 
        num4tt = join(['0',string(NumofTrials)],'') ; 
    elseif num ==1
        num4tt = '01' ;
    else
         num4tt = string(NumofTrials) ; 
    end
%     savfigname = join([strrep(categoryname,'_','-'),'-Trial',num4tt,'-OG.PNG']) ; 
    trialnum = num ; titlename_trial = join([strrep(categoryname, '_', '-'), ' Trial ', string(trialnum)], "") ;
    titlename = join([strrep(categoryname, '_', '-'),' Trial ',string(num)], "") ;
    uu = u{num} ; 
    tt = t{num} ; 

clear savfigname ; savfigname = 'VelocityProfilesZOOM2.png' ; 
savfigname = join([categoryname, 'Trial',num4tt,savfigname],'_') ; 
figure(cnt) ; cnt=cnt+1;
plot(tt,uu, '-o', 'LineWidth', 2)
legend(ADVlabel)
xlim([0, 5])
xlabel('$t[s]$')
ylabel('$u [m/s]$')
title(titlename)
set(gcf, 'Units','normalized','OuterPosition',[ 0 0 1 1]) ; 
saveas(gcf, fullfile(savfolderpath, savfigname)) ; close all
end 

end
