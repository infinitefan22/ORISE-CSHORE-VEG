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

currentfolder = pwd ; %MATLABONLINE
savfolderpath = join([currentfolder, '/data/'], '') ; 
cnt = 1;
clear fieldnames ; fieldnames = fieldnames(aalldata) ;
ADVlabel = {'ADV 2, $z=1.404$m', 'ADV 3, $z=1.550$m','ADV 4, $z=1.720$m','ADV 5, $z=1.858$m'} ; 
zu = [1.4040, 1.5500, 1.7200, 1.8580] ; % ADV 2,3,4,5 (order is correct and checked) ; from the bottom of the wave flume, not the false floor
ffh = .85 ; %false floor height

for totalnum = 7 %1:length(fieldnames) %CATEGORIES LOOP
    categoryname = fieldnames{totalnum} ; % 'HighDensity_h270_hv182_NoWall' ; %
    set_category_variables
    NumofTrials = length(t) ; 

    for num = 1:NumofTrials % TRIALS LOOP
        %% Correlation calculation
        
        for ADVnum = 1:width(u{1})
            clear x y 
            [x,y] = mscohere(c.eta_pp{num}, u{num}(:,ADVnum)) ; % max y = pi
            [ch.etavu{num}(ADVnum), idx] = max(x(1:50)) ; 
            ch.etavu_y{num}(ADVnum) = y(idx) ; 
        
        
            % ch.etavu = mscohere(c.eta_pp{1}, u{1}(:,ADVnum)) ; % use this format for all the correlations that would vary with depth
            % % ch.etacoshvu = mscohere(c.etacosh(:,ADVnum), variable2(:,ADVnum)) ; 
        
        
        end
        
        
    end
        figure(1)
        for nnn = 1:NumofTrials
        scatter(ch.etavu_y{nnn}, ch.etavu{nnn}, 'filled') ; hold on
        end
end




%% Plotting


%% EP_Cosh Correlation 
% coshratioPcorr = corr(c.coshratio(:),c.ucoshratio(:),'Type','Pearson')  
% coshratioScorr = corr(c.coshratio(:),c.ucoshratio(:),'Type','Spearman')  
% [coshratioXcorr, lags] = xcorr(c.coshratio,c.ucoshratio)  


