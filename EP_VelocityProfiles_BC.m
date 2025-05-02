% Velocity Profiles, Baseline Comparison
clearvars -except aalldata ; clc ; close all ; 
addpath('./ClaraFunctions') ; 
addpath('./data') ; 
addpath('./mfiles') ; 
     if ~exist('aalldata', "var") ;  load('aalldata_20250428.mat') ; end
 set(0,'defaultTextInterpreter','latex')
 set(groot,'defaultAxesTickLabelInterpreter','latex') 
 set(groot, 'defaultLegendInterpreter','latex')

currentfolder = pwd ; %MATLABONLINE
savfolderpath = join([currentfolder, '/ClaraFigures/EP_VelocityProfiles/20250501 plots/'], '') ; 

%        saveas(gcf, fullfile(savfolderpath, savfigname)) 
%%
cnt =1;
clear fieldnames ; fieldnames = fieldnames(aalldata) ;
zu = [1.4040, 1.5500, 1.7200, 1.8580] ; % ADV 2,3,4,5 (order is correct and checked)
catcolors = hsv(length(fieldnames)-4) ; 
for totalnum = 1:length(fieldnames)
    categoryname = fieldnames{totalnum} ; %'HighDensity_h270_hv182_NoWall' ; %
    set_category_variables
%% Length Dep Constants
NumofTrials = length(t) ; 
bmatch.(categoryname).match = ones(1,NumofTrials) ; 
clear trialnumbers vrms basetvrms
%% Start of Trial Indp Sections
for num = 1:NumofTrials
    trialnum = num ; titlename_trial = join([strrep(categoryname, '_', '-'), ' Trial ', string(trialnum)], "") ;
    titlename = join([strrep(categoryname, '_', '-')], "") ;
    tt = t{num} ; 
    uu = u{num} ; 
    [fbaseline, tnumb] = matchB_HDLDTrials(aalldata, categoryname, num, .5, .05) ; 
    if tnumb == 0 % will happen when there are no matches
        baset = zeros(height(uu),width(uu)) ; % 0000000s 
        bmatch.(categoryname).match(num) = 0 ; 
    else
        baset = aalldata.(fbaseline).u{tnumb} ; 
    end
%% Calculations
    vrms = zeros(1, width(uu)) ; basetvrms = zeros(1, width(uu)) ;% vrmstrial = zeros(num,width(uu)) ; 
    for ii = 1:width(uu)
        vrms(ii) = rms(uu(:,ii)) ;
        basetvrms(ii) = rms(baset(:,ii)) ; %has the same width as uu (all u has same width)
    end
    difvrmstrial1.(categoryname)(num, :) = vrms - basetvrms ; % * bmatch.(categoryname).match(num) 
    difvrmscategories1.(categoryname) = rms(difvrmstrial1.(categoryname)) ;
    difvmeancategories1.(categoryname) = mean(difvrmstrial1.(categoryname)) ;

    difvrmstrial2.(categoryname)(num, :) = vrms * bmatch.(categoryname).match(num) - basetvrms ; % 
    difvrmscategories2.(categoryname) = rms(difvrmstrial2.(categoryname)) ;
    difvmeancategories2.(categoryname) = mean(difvrmstrial2.(categoryname)) ;
   
    % vrmstrial.(categoryname)(num, :) = vrms ; 
    % vrmscategories.(categoryname) = rms(vrmstrial.(categoryname)) ;
    % vmeancategories.(categoryname) = mean(vrmstrial.(categoryname)) ;

trialnumbers(num) = join(['Trial ',string(num)], '') ; trialnumbers(end+1) = '$mean$ matched trials' ; 
trialcolors = hsv(NumofTrials) ; 
clear savfigname ; savfigname = join(['VelocityProfilesMatched_',categoryname, '.png'],'') ; % , 'Trial',string(num)
end 
%% Plotting 
% figure(cnt) ; cnt=cnt+1 ;
%     for num = 1:NumofTrials
%         if bmatch.(categoryname).match(num) == 0
%             scatter(difvrmstrial1.(categoryname)(num,:), zu, 55, trialcolors(num,:), 'x', 'filled', 'MarkerEdgeColor', trialcolors(num,:)) ; hold on 
%         else
%         scatter(difvrmstrial1.(categoryname)(num,:), zu, 36, trialcolors(num,:), 'o', 'filled', 'MarkerEdgeColor', "black") ; hold on 
%         end
%     end
%     % scatter(vrmscategories.(categoryname), zu, 100, 'black', 'filled') %RMS
%     scatter(difvmeancategories2.(categoryname), zu, 100, 'black', 'filled', 'MarkerEdgeColor', "black") %Mean not including x's
%     legend(trialnumbers)
%     ylim([1, 3])
%     xlabel('$u_{mean}$ (m/s)')
%     ylabel('Elevation (m)')
%     title(join([titlename, ' $v_{mean}$'],''))
%     set(gcf, 'Position', [10, 10, 600, 400]);
%     % saveas(gcf, fullfile(savfolderpath, savfigname))

end
%%
figure(100)
for totalnum = 5:length(fieldnames)
    categoryname = fieldnames{totalnum} ;
    set_category_variables
    NumofTrials = length(t) ; 
    bmatch.(categoryname).match = ones(1,NumofTrials) ; 
    % figure(100) ; 
    for num = 1:NumofTrials
        if bmatch.(categoryname).match(num) == 0
            scatter(difvrmstrial1.(categoryname)(num,:), zu, 55, catcolors(totalnum-4,:), 'x', 'filled', 'MarkerEdgeColor', white) ; hold on 
        else
        scatter(difvrmstrial1.(categoryname)(num,:), zu, 36, 'white', 'square', 'filled', 'MarkerEdgeColor', "white") ; hold on 
        end
    end
end

for totalnum = 5:length(fieldnames)
    categoryname = fieldnames{totalnum} ;
    set_category_variables
    NumofTrials = length(t) ; 
    bmatch.(categoryname).match = ones(1,NumofTrials) ; 
    gg1(totalnum-4) = scatter(difvmeancategories2.(categoryname), zu, 100, catcolors(totalnum-4,:), 'o', 'filled', 'MarkerEdgeColor', "black") ; %Mean not including x's
    ylim([1, 3])
    xlabel('$u_{mean}$ (m/s)')
    ylabel('Elevation (m)')
    title('HD/LD $v_{mean}$ of all trials')
    set(gcf, 'Position', [10, 10, 600, 400]);
    % saveas(gcf, fullfile(savfolderpath, savfigname))
end
legend(gg1, {'High Density h158 hv070', 'High Density h188 hv100', 'High Density h270 hv182', 'Low Density h158 hv073', 'Low Density h188 hv103', 'Low Density h233 hv148', 'Low Density h270 hv185'})
