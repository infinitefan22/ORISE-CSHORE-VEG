clearvars -except aalldata corrdata; clc ; close all ; 
addpath('./ClaraFunctions') ; 
addpath('./data') ; 
addpath('./mfiles') ; 
     if ~exist('aalldata', "var") ;  load('aalldatanewCdcalc_20250711.mat') ; end

set(0,'defaultTextInterpreter','latex')
set(groot,'defaultAxesTickLabelInterpreter','latex') 
set(groot, 'defaultLegendInterpreter','latex')

currentfolder = pwd ; %MATLABONLINE
savfolderpath = join([currentfolder, '/ClaraFigures/RootGraphs/'], '') ; 

%% Saving Data
savfig1name = 'CdCorrelationPlotsomezoom.png' ; 
%      saveas(gcf, fullfile(savfolderpath, savfig1name)) 
%% Setting Constants and Loading Data 
clear fieldnames ; fieldnames = fieldnames(aalldata) ;
clear graph legendnames
max=100 ; %x and y limits on graph
perfectcoor = linspace(0,max,10) ; 
graphcolors = lines(length(fieldnames)-4) ; 
legendnames = strrep(strrep(fieldnames, '_',' '), 'NoWall', '') ; 
legendnames{end+1} = 'Perfect Correlation' ; 
%% total category names loop
for totalnum = 5:length(fieldnames) %not including baseline sections

    categoryname = fieldnames{totalnum} ; %'HighDensity_h270_hv182_NoWall' ;
    set_category_variables
    if contains(categoryname, "HighDensity")|| contains(categoryname, "LowDensity")
        figure(1)
    %     for gnum = 1:length(graphcolors)
        for onum = 1:length(CdKelty)
        graph(totalnum-4) = scatter(CdKelty{onum}, Cdexact2{onum}, 100, graphcolors(totalnum-4,:), 'filled') ; hold on %graph(totalnum-4) = 
        end
    %     end
         xlim([0,20])
         ylim([0,50])
        xlabel('Kelty Cd','fontsize',16)
        ylabel('Johnson Cd','fontsize',16)
    %     legend(legendnames(5:end)) ;  
        title('Cd Calculation Correlation')
    end
end

figure(1)
graph(end+1) = plot(perfectcoor, perfectcoor, 'LineWidth', 2) ; hold on
    legend([graph], legendnames(5:end))
    