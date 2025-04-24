clearvars -except aalldata ; clc ; close all ; 
addpath('./ClaraFunctions') ; 
addpath('./data') ; 
addpath('./mfiles') ; 
     if ~exist('aalldata', "var") ;  load('aalldata_20250421.mat') ; end

 set(0,'defaultTextInterpreter','latex')
 set(groot,'defaultAxesTickLabelInterpreter','latex') 
 set(groot, 'defaultLegendInterpreter','latex')

% savfolderpath = '/home/elizabeth/Desktop/cshorex-main/osu_mangrove/ClaraFigures/EP_ADVPlotsApr102025/'; 
currentfolder = pwd ; %MATLABONLINE
savfolderpath = join([currentfolder, '/ClaraFigures/20250424_VelocityProfiles/'], '') ; 

%        saveas(gcf, fullfile(savfolderpath, savfigname)) 
% set(0, 'DefaultAxesFontName', 'Nimbus Roman', 'DefaultTextFontName', 'Nimbus Roman')
%%
cnt =1;
clear fieldnames ; fieldnames = fieldnames(aalldata) ;
zu = [1.4040, 1.5500, 1.7200, 1.8580] ; % ADV 2,3,4,5 (order is correct and checked)
for totalnum = 1:length(fieldnames)
    categoryname = fieldnames{totalnum} ; %'HighDensity_h270_hv182_NoWall' ; %
    set_category_variables
%% Length Dep Constants
NumofTrials = length(t) ; 
clear trialnumbers vrms
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
    vrmstrial.(categoryname)(num, :) = vrms ; 
    vrmscategories.(categoryname) = rms(vrmstrial.(categoryname)) ;
    vmeancategories.(categoryname) = mean(vrmstrial.(categoryname)) ;

% plotting inv trials

trialnumbers(num) = join(['Trial ',string(num)], '') ; trialnumbers(end+1) = '$rms$ all trials' ; 
trialcolors = hsv(NumofTrials) ; 
clear savfigname ; savfigname = 'VelocityProfilesRMS_' ; 
savfigname = join([savfigname,categoryname, '.png'],'') ; % , 'Trial',string(num)
end 
%% Plotting 
figure(cnt) ; cnt=cnt+1 ;
    for num = 1:NumofTrials
        scatter(vrmstrial.(categoryname)(num,:), zu, 36, trialcolors(num,:), 'filled') ; hold on 
    end
    scatter(vrmscategories.(categoryname), zu, 100, 'black', 'filled') %RMS
    % scatter(vmeancategories.(categoryname), zu, 100, 'black', 'filled') %Mean
    legend(trialnumbers)
    ylim([1, 3])
    xlabel('$u_{rms}$ (m/s)')
    ylabel('Elevation (m)')
    title(join([titlename, ' $v_{rms}$'],''))
    set(gcf, 'Position', [10, 10, 600, 400]);
    saveas(gcf, fullfile(savfolderpath, savfigname))
end
