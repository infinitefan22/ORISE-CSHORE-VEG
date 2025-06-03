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
savfolderpath = join([currentfolder, '/ClaraFigures/CorrelationCalculation/20250603/'], '') ; 
cnt = 1;
clear fieldnames ; fieldnames = fieldnames(aalldata) ;
ADVlabel = {'ADV 2, $z=1.404$m', 'ADV 3, $z=1.550$m','ADV 4, $z=1.720$m','ADV 5, $z=1.858$m'} ;
graphcolors = hsv(7) ; 
zu = [1.4040, 1.5500, 1.7200, 1.8580] ; % ADV 2,3,4,5 (order is correct and checked) ; from the bottom of the wave flume, not the false floor
ffh = .85 ; %false floor height
indvlayoutplot = 1 ; %plotting indv layout graphs. 1 for yes, anything else for no (0)

for totalnum = 5:length(fieldnames) %CATEGORIES LOOP
    clear ch specADVval
    categoryname = fieldnames{totalnum} ; % 'HighDensity_h270_hv182_NoWall' ; %
    set_category_variables
    NumofTrials = length(t) ; 
%% Correlation calculation
    for num = 1:NumofTrials % TRIALS 
        % milinx = ; %idx where y is .01 sec
        for ADVnum = 1:width(u{1})
            clear x y 
            [x,y] = mscohere(c.eta_pp{num}, u{num}(:,ADVnum)) ; % max y = pi
            [ch.etavu{num}(ADVnum), idx] = max(x(1:100)) ; 
            ch.etavu_y{num}(ADVnum) = y(idx) ; 
            % ch.etavu = mscohere(c.eta_pp{1}, u{1}(:,ADVnum)) ; % use this format for all the correlations that would vary with depth
            % % ch.etacoshvu = mscohere(c.etacosh(:,ADVnum), variable2(:,ADVnum)) ; 
        end
    end

    for ADVnnn = 1:4
        for num = 1:NumofTrials
            specADVval(num) = ch.etavu{num}(ADVnnn) ; 
        end
        ch.mean{ADVnnn} = mean(specADVval) ; 
        % if contains(categoryname, '158') || contains(categoryname, '188')
        %     ch.mean{3} = NaN ; 
        %     ch.mean{4} = NaN ; 
        % end
    end
    %% Plotting prep
    if contains(categoryname, "High")
        shapetype = 'o' ; 
        shorttit = 'High Density ' ; 
    elseif contains(categoryname, "Low")
        shapetype = 'square' ; 
        shorttit = 'Low Density ' ; 
    else 
        shapetype = '^' ; 
        shorttit = '' ;
    end

 %% Plotting Indv Layouts (like EP_VelocityProfiles)
    if indvlayoutplot == 1  
        figure(cnt) ; cnt = cnt+1 ;
        trialcolors = hsv(NumofTrials) ; 
        clear trialnumbers
        for num = 1:NumofTrials
            trialnumbers(num) = join(['Trial ',string(num)], '') ; trialnumbers(end+1) = '$mean$ all trials' ; 
            scatter(ch.etavu{num}, zu, 36, trialcolors(num,:), 'o', 'filled', 'MarkerEdgeColor', 'black') ; hold on 
        end
        % scatter(vrmscategories.(categoryname), zu, 100, 'black', 'filled') %RMS
        scatter([ch.mean{:}], zu, 100, 'black', 'filled', 'MarkerEdgeColor', "black") ; hold on %Mean
        legend(trialnumbers)
        xlim([0 1.01])
        xlabel('Coherence')
        ylim([1.2 2])
        yticks(zu) ; yticklabels(string(zu)) 
        ylabel('$z$ [m]')
        title(strrep(strrep(categoryname, '_NoWall', ''), '_', '-') + " Correlation $eta$ vs. $u$")
        set(gcf, 'Position', [10, 10, 600, 400]);
        % saveas(gcf, fullfile(savfolderpath, savfigname))
    end
%% All Combined
    % figure(100) ;
    % for nnn = 1:NumofTrials
    % scatter(ch.etavu{nnn}, zu, 36, graphcolors(totalnum-4,:), 'filled', shapetype, 'MarkerEdgeColor', 'k') ; hold on
    % end
    % gg1 = scatter([ch.mean{:}], zu, 100, graphcolors(totalnum-4,:), 'filled',shapetype, 'MarkerEdgeColor', "black") ; hold on %Mean
    % if ~exist('graphes','var') ; graphes(1) = gg1 ; else ; graphes(end+1) = gg1 ; end
    % if ~exist('legendentry','var') ; legendentry{1} = strrep(strrep(categoryname, '_NoWall', ''), '_', '-') ; else;legendentry{end+1} = strrep(strrep(categoryname, '_NoWall', ''), '_', '-') ; end 
    % title("All Layout Correlation $eta$ vs. $u$")
    % xlim([0 1.01])
    % xlabel('Coherence')
    % ylim([1.2 2])
    % yticks(zu) ; yticklabels(string(zu)) 
    % ylabel('$z$ [m]')
    % % legend(graphes, fieldnames(5:end))
    %  %see EP_Cosh for reference on legend for multiple trials but one shape being graphed 
end

% figure(100) ;  legend(graphes, legendentry) %fieldnames(5:end))
savfigname = "AllCorrelations_EtavU.png" ;
  % exportgraphics(gcf, fullfile(savfolderpath, savfigname), 'Resolution', 300) ; %close all


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