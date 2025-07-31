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
keltypaperCd
%      saveas(gcf, fullfile(savfolderpath, savfig1name)) 
%% Setting Constants and Loading Data 
clear fieldnames ; fieldnames = fieldnames(aalldata) ;
clear graph legendnames
max=100 ; %x and y limits on graph
perfectcoor = linspace(0,max,10) ; 
graphcolors = hsv(length(fieldnames)-4) ; 
legendnames = strrep(strrep(fieldnames, '_',' '), 'NoWall', '') ; 
legendnames{end+1} = 'Perfect Correlation' ; 
CdtableRand = readtable('Excel_SummaryTable_Random.xlsx') ; 
CdtableReg = readtable('Excel_SummaryTable_Regular.xlsx') ;
cnt =1; 
loopmin = 10 ; 
%% total category names loop
for totalnum = loopmin:length(fieldnames) %not including baseline sections
    clear keltyCdRand
    categoryname = fieldnames{totalnum} ; %'HighDensity_h270_hv182_NoWall' ;
    set_category_variables
    keltyCd = keltypaper.(categoryname).Cd ; 
    keltyerr = keltypaper.(categoryname).Cd ;
    disp("number " + categoryname + " kelty Cd length = " + length(keltyCd))
    disp("number " + categoryname + " Johnson Cd length = " + length(Cdexact2))
    disp(" ")
    cnt = cnt+1 ;
    if contains(categoryname, "HighDensity")|| contains(categoryname, "LowDensity")
        figure(1)
    %     for gnum = 1:length(graphcolors)
        for onum = 1:length(keltyCd) %length(CdKelty)
        graph(totalnum-(loopmin-1)) = scatter(keltyCd(onum), Cdexact2{onum}, 100, graphcolors(totalnum-(loopmin-1),:), 'filled') ; hold on %graph(totalnum-4) = 
        errorbar(keltyCd(onum), Cdexact2{onum}, keltyerr(onum), "horizontal", 'Color', graphcolors(totalnum-(loopmin-1),:))
        end
    %     end
         xlim([0,10])
         ylim([0,10])
        xlabel('Kelty Cd','fontsize',16)
        ylabel('Johnson Cd','fontsize',16)
    %     legend(legendnames(5:end)) ;  
        title('Cd Calculation Correlation')
    end
end

figure(1)
graph(end+1) = plot(perfectcoor, perfectcoor, 'LineWidth', 2, 'Color', '#c0c0c0') ; hold on
    legend([graph], legendnames(loopmin:end))
    