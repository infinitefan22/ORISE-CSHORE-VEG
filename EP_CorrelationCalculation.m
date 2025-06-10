% Pearson Correlation - good for linear relationships, if data goes up or down ONLY and steadly
% Spearman Correlation - ranks data to see if data sets go up and down together
% Kendall Correlation - same as Spearman, but used for smaller than n=30 data sets
% Cross Correlation - time lagged correlation 
% Coherance for Signal Processing - tells how well the frequency of signals
% match***

clearvars -except aalldata corrdata; clc ; close all ; 
addpath('./ClaraFunctions') ; 
addpath('./data') ; 
addpath('./mfiles') ; 
     if ~exist('aalldata', "var") ;  load('aalldata_20250428.mat') ; end
     if ~exist('corrdata', "var") ;  load('corrdata_20250529.mat') ; end

set(0,'defaultTextInterpreter','latex')
set(groot,'defaultAxesTickLabelInterpreter','latex') 
set(groot, 'defaultLegendInterpreter','latex')

currentfolder = pwd ; %MATLABONLINE
savfolderpath = join([currentfolder, '/ClaraFigures/CorrelationCalculation/20250609/'], '') ; 
cnt = 1;
clear fieldnames ; fieldnames = fieldnames(aalldata) ;
ADVlabel = {'ADV 2, $z=1.404$m', 'ADV 3, $z=1.550$m','ADV 4, $z=1.720$m','ADV 5, $z=1.858$m'} ;
graphcolors = hsv(7) ; 
zu = [1.4040, 1.5500, 1.7200, 1.8580] ; % ADV 2,3,4,5 (order is correct and checked) ; from the bottom of the wave flume, not the false floor
ffh = .85 ; %false floor height
% indvlayoutplot = 0 ; %plotting indv layout graphs. 1 for yes, anything else for no (0)
savfigs = 0;
 % errbars = 1 ; 
Qcosh = 1 ; %do you want to see eta v velocity (0) or eta*cosh v velocity (1)

for totalnum = 5:length(fieldnames) %CATEGORIES LOOP
    gcnum = totalnum-4 ; %should be changed with graphcolors to avoid error
    clear ms specADVvalms
    categoryname = fieldnames{totalnum} ; % 'HighDensity_h270_hv182_NoWall' ; %
    set_category_variables
    NumofTrials = length(t) ; 
%% Correlation calculation
    if Qcosh ==1 ; Qevc = 'etacosh' ; else ; Qevc = 'eta_pp' ; end
    
    if contains(categoryname, '158') || contains(categoryname, '188')
        totADV = 2 ; 
    else 
        totADV = 4 ; 
    end

    for num = 1:NumofTrials % Coherence ; TRIALS % milinx = ; %idx where y is .01 sec
        for ADVnum = 1:totADV
            clear x y idx
            [x,y] = mscohere(c.(Qevc){num}(:,ADVnum), u{num}(:,ADVnum)) ; % max y = pi
            [ms.etavu{num}(ADVnum), idx] = max(x(1:100)) ; 
            ms.etavu_y{num}(ADVnum) = y(idx) ; 
            % ch.etavu = mscohere(c.eta_pp{1}, u{1}(:,ADVnum)) ; % use this format for all the correlations that would vary with depth
            % % ch.etacoshvu = mscohere(c.etacosh(:,ADVnum), variable2(:,ADVnum)) ; 
        end
    end
    for num = 1:NumofTrials % Cross Correlation
        for ADVnum = 1:totADV
            clear x y idx
            [x, y] = xcorr(c.(Qevc){num}(:,ADVnum), u{num}(:,ADVnum),1000,'normalized') ; 
            [xc.etavu{num}(ADVnum), idx] = max(x) ; 
            xc.etavu_y{num}(ADVnum) = y(idx) ; 
        end
    end

    for ADVnnn = 1:totADV %Calculating Mean of Coherence & Cross Correlation
        for num = 1:NumofTrials
            specADVvalms(num) = ms.etavu{num}(:,ADVnnn) ; 
            ms.sdmean(ADVnnn) = std(specADVvalms) ; 
            specADVvalxc(num) = xc.etavu{num}(:,ADVnnn) ; 
            xc.sdmean(ADVnnn) = std(specADVvalxc) ; 
        end
        ms.mean{ADVnnn} = mean(specADVvalms) ; 
        xc.mean{ADVnnn} = mean(specADVvalxc) ; 
        clear specADVvalms specADVvalxc
    end
%% 

    %% Plotting 
     Coherence_Indv_Plots
     % Coherence_Plots
     % CrossCorrelation_Plots

end

% figure(100) ;  legend(graphes, legendentry) %fieldnames(5:end))
savfigname = "AllCorrelations_Cross_EtavU.png" ;
  %exportgraphics(gcf, fullfile(savfolderpath, savfigname), 'Resolution', 300) ; %close all


%% Plotting
% % Making sure correlation occurs near wave maker
% figure(cnt) ; cnt = cnt+1 ;
% for nnn = 1:NumofTrials
% scatter(ch.etavu_y{nnn}, ch.etavu{nnn}, 'filled') ; hold on
% end

%% EP_Cosh Correlation 
% coshratioPcorr = corr(c.coshratio(:),c.ucoshratio(:),'Type','Pearson')  
% coshratioScorr = corr(c.coshratio(:),c.ucoshratio(:),'Type','Spearman')  
% [coshratioXcorr, lags] = xcorr(c.coshratio,c.ucoshratio)  

% figure(200)
% stem(xy,xx)