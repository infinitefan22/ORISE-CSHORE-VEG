%  clear ; clc ; close all ; 
addpath('./ClaraFunctions') ; 
addpath('./data') ; 
addpath('./mfiles') ; 
%      load('aalldata_Mar282025.mat') ;
 savfigname = 'NewADVVelocityPlotsAllLayoutsWide.png' ; 
savfolderpath = '/home/elizabeth/Desktop/cshorex-main/osu_mangrove/ClaraFigures/EP_ADVPlotsApr092025/'; 
%        saveas(gcf, fullfile(savfolderpath, savfigname)) 
% set(0, 'DefaultAxesFontName', 'Nimbus Roman', 'DefaultTextFontName', 'Nimbus Roman')
 set(0,'defaultTextInterpreter','latex')
 set(groot,'defaultAxesTickLabelInterpreter','latex') 
 set(groot, 'defaultLegendInterpreter','latex')
 %%
clear fieldnames ; fieldnames = fieldnames(aalldata) ;
for totalnum = 1:length(fieldnames)
    categoryname = fieldnames{totalnum} ; %'HighDensity_h270_hv182_NoWall' ;
    set_category_variables
%% Length Dep Constants
NumofTrials = length(t) ; 
%% Start of Trial Indp Sections
for num = 1:NumofTrials
    trialnum = num ; titlename_trial = join([strrep(categoryname, '_', '-'), ' Trial ', string(trialnum)], "") ;
    titlename = join([strrep(categoryname, '_', '-')], "") ;
    savfig1name = join([categoryname,'_Trial',string(num),'_ADV.png'],'') ; 
    
    uu = u{num} ; 
    tt = t{num} ; 
%% Calculations
    vrms = zeros(1, width(uu)) ; % vrmstrial = zeros(num,width(uu)) ; 
    for ii = 1:width(uu)
        vrms(ii) = rms(uu(:,ii)) ;
    end
    vrmstrial(num, :) = vrms ; 
    vrmscategories.(categoryname) = rms(vrmstrial) ;
    stdvcategories.(categoryname) = std(vrmstrial) ; 
end
%% preplotting 
if contains(categoryname, 'h158')
    subp = 1 ; 
    subptitle = 'Water Height $h = 158$' ;
elseif contains(categoryname, 'h188')
    subp = 2 ; 
    subptitle = 'Water Height $h = 188$' ;
elseif contains(categoryname, 'h233')
    subp=3 ;
    subptitle = 'Water Height $h = 233$' ;
elseif contains(categoryname, 'h270') 
    subp=4 ; 
    subptitle = 'Water Height $h = 270$' ;
end
if contains(categoryname, 'Baseline')
    xtick = 1 ; 
elseif contains(categoryname, 'High')
    xtick = 2 ; 
elseif contains(categoryname,'Low')
    xtick = 3; 
end

%% plotting
figure(100) ;  %can't put cnt+1 here, so need at start of next fig (cause of subplots)
    sgtitle('ADV Velocity_{rms}') 
    subplot(4,1,subp)
    scatter(xtick*ones(1,1),vrmscategories.(categoryname), 36, 'filled') ; hold on
    errorbar(xtick,vrmscategories.(categoryname),stdvcategories.(categoryname)) ; hold on
    xlim([0, 4]) ; xticks(0:1:NumofTrials+1)
%     ylim([0,.8]) ; yticks(0:.2:.8) ; %make sure the max are the same for ylim and yticks
    xticklabels({' ','Baseline', 'High Density' , 'Low Density'}) %,'fontsize',8)
    ylabel('$u_{rms}$ [m/s]') %,'fontsize',8)
    title(subptitle)
end