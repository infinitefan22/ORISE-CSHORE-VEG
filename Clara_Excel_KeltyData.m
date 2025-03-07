 clear ;  close all ;  clc ;
addpath('./data') ;
addpath('./ClaraFunctions') ;
%1 if you want to plot
plottingKelty = 0 ; plottingall = 0;
cnt = 1 ; %so I don't need to count number of figures (and fig numbers aren't definied)
%% Cdtable
CdtableRand = readtable('Excel_SummaryTable_Random.xlsx') ; 
CdtableReg = readtable('Excel_SummaryTable_Regular.xlsx') ;
[RandHD, yerrRandHDtop, yerrRandHDbot, RandLD, yerrRandLDtop, yerrRandLDbot] = definetableHDLD(CdtableRand, 'yes') ; 
[RegHD,  yerrRegHDtop,  yerrRegHDbot,  RegLD,  yerrRegLDtop,  yerrRegLDbot ] = definetableHDLD(CdtableReg, 'yes') ; 
CdAll = [CdtableRand.CD; CdtableReg.CD] ; 
ReAll = [CdtableRand.ReUde; CdtableReg.ReUde] ;
%% line of best fit equs
ReequRand = linspace(min(CdtableRand.ReUde)-1000, max(CdtableRand.ReUde)+1000, 1000) ; 
CdequRand = .61 + (10600./ReequRand).^1.26 ; 
ReequReg = linspace(min(CdtableReg.ReUde)-1000, max(CdtableReg.ReUde)+1000, 1000) ; 
CdequReg = lineofbestfit_fract_exp(CdtableReg.ReUde, CdtableReg.CD, length(ReequReg), 1000, 1000) ;
ReequAll = linspace(min(ReAll)-1000, max(ReAll)+1000, 1000) ; 
CdequAll = lineofbestfit_fract_exp(ReAll, CdAll, length(ReequAll), 1000, 1000) ;

%%
if plottingKelty == 1 
figure(cnt) % Regular and Random Wave Data Re
    cnt = cnt+1 ;
    scatter( RegHD.ReUde,RegHD.CD, 'filled','o', 'blue') ; hold on 
    scatter( RegLD.ReUde,RegLD.CD, 'filled','^','blue') ; hold on
    plot(ReequReg, CdequReg, 'Color', '#000080', 'LineWidth', 2) 
    scatter( RandHD.ReUde,RandHD.CD, 'filled','o', 'red') ; hold on 
    scatter( RandLD.ReUde,RandLD.CD, 'filled','^','red') ; hold on
    plot(ReequRand, CdequRand, 'Color', "#A2142F", 'LineWidth', 2)
    plot(ReequAll, CdequAll, 'Color', 'k', 'LineWidth', 2, 'LineStyle', '--')  
%     xlim([0,3e4]) ; ylim([0,6])
    title("Kelty Cd Plotted Against Reynold's Number Regular and Random Waves")
    xlabel("Reynold's Number (Re)") ; ylabel("Drag Coefficient (Cd)") ; 
    legend({'Regular Waves HD', 'Regular Waves LD', 'Line of Best Fit Regular Waves', ...
        'Random Waves HD', 'Random Waves LD', 'Line of Best Fit Random Waves', ...
        'Line of Best Fit All Data'}, 'Location', 'best') ; 
figure(cnt) % Regular Wave Data Re
    cnt = cnt+1 ;
    scatter( RegHD.ReUde,RegHD.CD, 'filled','o', 'b') ; hold on 
    errorbar(RegHD.ReUde,RegHD.CD, yerrRegHDbot, yerrRegHDtop,'LineStyle','none','Color', 'b', 'CapSize', 10)
    scatter( RegLD.ReUde,RegLD.CD, 'filled','^','b') ; hold on
    errorbar(RegLD.ReUde,RegLD.CD, yerrRegLDbot, yerrRegLDtop,'LineStyle','none','Color', 'b', 'CapSize', 10)
    plot(ReequReg, CdequReg, 'Color', '#000080', 'LineWidth', 2) 
%     xlim([0,3e4]) ; ylim([0,6])
    title("Kelty Cd Plotted Against Reynold's Number Regular Waves")
    xlabel("Reynold's Number (Re)") ; ylabel("Drag Coefficient (Cd)") ; 
    legend({'High Density','HD sd', 'Low Density', 'LD sd', 'Line of Best Fit'}, 'Location', 'best') ; 
figure(cnt) %Random Wave Data Re
    cnt = cnt+1 ;
    scatter( RandHD.ReUde,RandHD.CD, 'filled','o', 'r') ; hold on 
    errorbar(RandHD.ReUde,RandHD.CD, yerrRandHDbot, yerrRandHDtop,'LineStyle','none','Color', 'r', 'CapSize', 10)
    scatter( RandLD.ReUde,RandLD.CD, 'filled','^','r') ; hold on
    errorbar(RandLD.ReUde,RandLD.CD, yerrRandLDbot, yerrRandLDtop,'LineStyle','none','Color', 'r', 'CapSize', 10)
    plot(ReequRand, CdequRand, 'Color', "#A2142F", 'LineWidth', 2)
    xlim([0,3e4]) ; ylim([0,6])
    title("Kelty Cd Plotted Against Reynold's Number Random Waves")
    xlabel("Reynold's Number (Re)") ; ylabel("Drag Coefficient (Cd)") ; 
    legend({'High Density','HD sd', 'Low Density', 'LD sd', 'Line of Best Fit'}, 'Location', 'best') ; 
figure(cnt) %Random Wave Data KC
    cnt = cnt+1 ;
    scatter( RandHD.KCUde,RandHD.CD, 'filled','o', 'r') ; hold on 
%     errorbar(CdtableHD.KCUde,CdtableHD.CD, yerrHDbot, yerrHDtop,'LineStyle','none','Color', 'b', 'CapSize', 10)
    scatter( RandLD.KCUde,RandLD.CD, 'filled','^','r') ; hold on
%     errorbar(CdtableLD.KCUde,CdtableLD.CD, yerrLDbot, yerrLDtop,'LineStyle','none','Color', 'r', 'CapSize', 10)
    xlim([0,150]) ; ylim([0,6])
    title("Kelty Cd Plotted Against Keulegan-Carpenter Number Random Waves")
    xlabel("Keulegan-Carpenter Number(KC)") ; ylabel("Drag Coefficient (Cd)") ;
    legend({'High Density', 'Low Density'}, 'Location', 'best') ;
end

% %% kelty large table
% keltydataog = readtable('Experimental_Info.xlsx') ;
% %% Impt Locations
% Array1col = 22 ; % with 17 col
% Array2col = 121 ; %with 18 col
% 
% %% Processing Data, creating specific trial sections
% trialnames = unique(keltydataog.ExperimentFileName) ;
%     namecounts = histcounts(categorical(keltydataog{:,1}), categorical(trialnames));
%     trialnames = trialnames(namecounts > 1) ;
%     trialnames = trialnames(~cellfun('isempty',trialnames)) ; trialnames = trialnames(~cellfun(@(x)contains(x,'comment'),trialnames)) ;
% keltydata = struct();% Create a structure to store the tables for each trial
% 
% for i = 1:height(trialnames) % Loop over each unique trial name
%     trialdata = keltydataog(strcmp(keltydataog.ExperimentFileName, trialnames{i}), :);% Get the rows where the first column matches the current trial name
%     keltydata.(trialnames{i}) = trialdata;% Store the data in the structure with the trial name as the field
% %     disp(['Stored table for ', trialnames{i}]);% Optionally display a message for each trial
% end
% 
% trialnames = fieldnames(keltydata) ; %redifining trail names
% trialnamesBaseline = trialnames(~cellfun(@(x)contains(x,'Baseline'),trialnames)) ;
% trialnamesWall = trialnames(~cellfun(@(x)contains(x,'_Wall'),trialnames)) ;
% trialnamesNoWall = trialnames(~cellfun(@(x)contains(x,'_NoWall'),trialnames)) ;

