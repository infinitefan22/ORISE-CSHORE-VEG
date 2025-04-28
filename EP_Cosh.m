clearvars -except aalldata ; clc ; close all ; 
addpath('./ClaraFunctions') ; 
addpath('./ClaraFunctions/LOESS regression smoothing-2.1.0.0') ;
addpath('./data') ; 
addpath('./mfiles') ; 
     if ~exist('aalldata', "var") ;  load('aalldata_20250428.mat') ; end

currentfolder = pwd ; %MATLABONLINE
savfolderpath = join([currentfolder, '/ClaraFigures/EP_Cosh/20250428/'], '') ; 

set(0,'defaultTextInterpreter','latex')
set(groot,'defaultAxesTickLabelInterpreter','latex') 
set(groot, 'defaultLegendInterpreter','latex')
%% Constants
cnt = 1;
clear fieldnames ; fieldnames = fieldnames(aalldata) ;
ADVlabel = {'ADV 2, $z=1.404$m', 'ADV 3, $z=1.550$m','ADV 4, $z=1.720$m','ADV 5, $z=1.858$m'} ; 
zu = [1.4040, 1.5500, 1.7200, 1.8580] ; % ADV 2,3,4,5 (order is correct and checked)
catcolors = hsv(length(fieldnames)) ; 

for totalnum = 1:length(fieldnames)
    categoryname = fieldnames{totalnum} ; % 'HighDensity_h270_hv182_NoWall' ; %
    set_category_variables
    
%% Length Dep Constants
NumofTrials = length(t) ; 
%% Start of Trial Indp Sections
for num = 1:NumofTrials  %
    clear hh ucosh1 ucosh2 ucosh1rums ucosh2rms
    if num == NumofTrials && NumofTrials <10 
        num4tt = join(['0',string(NumofTrials)],'') ; 
    elseif num == 1
        num4tt = '01' ;
    else
         num4tt = string(NumofTrials) ; 
    end

    trialnum = num ; titlename_trial = join([strrep(categoryname, '_', '-'), ' Trial ', string(trialnum)], "") ;
    titlename = join([strrep(categoryname, '_', '-'),' Trial ',string(num)], "") ;
    uu = u{num} ; 
    tt = t{num} ; 

    hh = hv{num} ; % still water height
    cosh1 = cosh(k{num}.*(hh-zu(4))) ; 
    cosh2 = cosh(k{num}.*(hh-zu(1))) ; 
    ucosh1 = uu(:,4) .* cosh1 ;
    ucosh2 = uu(:,1) .* cosh2 ; 
if contains(categoryname, 'h158')
    uu(:,3) = zeros(height(uu),1) ;  
    uu(:,4) = zeros(height(uu),1) ; 
    ucosh1 = uu(:,2) .* cosh(k{num}.*(hh-zu(2))) ;
end

ucosh1rms = rms(ucosh1) ; 
ucosh2rms = rms(ucosh2) ; 
ucoshratio = ucosh1rms/ucosh2rms ; 
coshratio = cosh1./cosh2 ; 

% if contains(categoryname, 'Basename')
%     ucrbase(totalnum, num) = ucoshratio ; 
%     xcrbase(totalnum, num) = coshratio ;
% elseif contains(categoryname, 'HighDensity')
%     ucrhigh(totalnum-4, num) = ucoshratio ; 
%     xcrhigh(totalnum-4, num) = coshratio ;
% elseif contains(categoryname, 'LowDensity')
%     ucrlow(totalnum-7, num) = ucoshratio ; 
%     xcrlow(totalnum-7, num) = coshratio ;
% end

clear savfigname ; savfigname = 'CoshDependence.png' ; savfigname = join([categoryname, 'Trial',num4tt,savfigname],'_') ; 
%% plotting
figure(1) ;
if contains(categoryname, "Base")
    g1 = scatter(coshratio,ucoshratio, 45, catcolors(totalnum,:), "^", 'filled') ; hold on
elseif contains(categoryname, "High")
    g2 = scatter(coshratio,ucoshratio, 45, catcolors(totalnum,:), "o", 'filled') ; hold on
elseif contains(categoryname, "Low")
    g3 = scatter(coshratio,ucoshratio, 45, catcolors(totalnum,:), "square", 'filled') ; hold on
end
% if contains(categoryname, "Base")
%     g1 = scatter(xcrbase,ucrbase, 36, catcolors(totalnum,:), "^", 'filled') ; hold on
% elseif contains(categoryname, "High")
%     g2 = scatter(xcrhigh,ucrhigh, 36, catcolors(totalnum,:), "o", 'filled') ; hold on
% elseif contains(categoryname, "Low")
%     g3 = scatter(xcrlow,ucrlow, 36, catcolors(totalnum,:), "square", 'filled') ; hold on
% end


xlim([0,3])
xlabel('$cosh_1/cosh_2$')
ylim([0,4.5])
ylabel('$u_{1,rms}/u_{2,rms}$')
% title("Cosh Dependency")
% set(gcf, 'Units','normalized','OuterPosition',[ 0 0 1 1]) ; 
set(gcf, 'Position', [100, 100, 1200, 800]);

 % saveas(gcf, fullfile(savfolderpath, savfigname)) ; %close all
end 

end

xxx = linspace(0,10,10) ; 
yyy = linspace(0,10,10) ; 

figure(1)
g4 = plot(xxx, yyy) ; hold on
legend([g4, g1, g2, g3] , {"1:1", "Baseline", "High Density", "Low Density"})
