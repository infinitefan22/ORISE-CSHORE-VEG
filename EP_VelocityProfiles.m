%  clear ; clc ; close all ; 
addpath('./ClaraFunctions') ; 
addpath('./data') ; 
addpath('./mfiles') ; 
%     load('aalldata_Mar282025.mat') ;
 savfigname = 'NewADVVelocityProfiles.png' ; 
savfolderpath = '/home/elizabeth/Desktop/cshorex-main/osu_mangrove/ClaraFigures/EP_ADVPlotsApr102025/'; 
%        saveas(gcf, fullfile(savfolderpath, savfigname)) 
% set(0, 'DefaultAxesFontName', 'Nimbus Roman', 'DefaultTextFontName', 'Nimbus Roman')
 set(0,'defaultTextInterpreter','latex')
 set(groot,'defaultAxesTickLabelInterpreter','latex') 
 set(groot, 'defaultLegendInterpreter','latex')
 cnt =1;
 %%
clear fieldnames ; fieldnames = fieldnames(aalldata) ;
zu = [1.4040, 1.5500, 1.7200, 1.8580] ; % ADV 2,3,4,5 (order is correct and checked)
for totalnum = 1:length(fieldnames)
    categoryname = fieldnames{totalnum} ; %'HighDensity_h270_hv182_NoWall' ; %
    set_category_variables
%% Length Dep Constants
NumofTrials = length(t) ; 
%% Start of Trial Indp Sections
for num = 1:NumofTrials
    trialnum = num ; titlename_trial = join([strrep(categoryname, '_', '-'), ' Trial ', string(trialnum)], "") ;
    titlename = join([strrep(categoryname, '_', '-')], "") ;
    uu = u{num} ; 
    tt = t{num} ; 
%% Calculations
    vrms = zeros(1, width(uu)) ; % vrmstrial = zeros(num,width(uu)) ; 
    for ii = 1:width(uu)
        vrms(ii) = rms(uu(:,ii)) ;
    end
    vrmstrial(num, :) = vrms ; 
    vrmscategories.(categoryname) = rms(vrmstrial) ;

% plotting inv trials
trialnumbers(num) = join(['Trial ',string(num)], '') ; trialnumbers(end+1) = '$rms$ all trials' ; 
trialcolors = hsv(NumofTrials) ; 
clear savfigname ; savfigname = 'VelocityProfiles.png' ; 
savfigname = join([categoryname, 'Trial',string(num),savfigname],'_') ; 
end 
%% Plotting 
figure(cnt) ; cnt=cnt+1 ;
    for num = 1:NumofTrials
        scatter(vrmstrial(num,:), zu, 36, trialcolors(num,:), 'filled') ; hold on 
    end
    scatter(vrmscategories.(categoryname), zu, 100, 'black', 'filled')
    legend(trialnumbers)
    ylim([1, 3])
    xlabel('$u_{rms}$')
    ylabel('Depth')
    title(titlename)
    saveas(gcf, fullfile(savfolderpath, savfigname))
end
