% this code uses the top most ADV that is collecting data and the 2nd from
% the bottom ADV. If there is no ADV available higher than the 2nd from
% bottom (h158), the bottom ADV will be compared to the 2nd bottom
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
catcolors = hsv(3) ; 
khcolors = winter(300) ; 
for totalnum = 1:length(fieldnames)
    categoryname = fieldnames{totalnum} ; % 'HighDensity_h270_hv182_NoWall' ; %
    set_category_variables
    kvalues = [aalldata.(categoryname).k{1:end}] ; %yes you need these lines. it's annoying
    hvvalues = [aalldata.(categoryname).hv{:}] ; 
    kh(totalnum,1) = max(kvalues)*max(hvvalues) ;
    kh(totalnum,2) = min(kvalues)*min(hvvalues) ; 
end ; clear totalnum
kmapmax = max(kh(:,1))  ;
kmapmin = min(kh(:,2))  ;
kfactor = (kmapmax)/length(khcolors) ; %-kmapmin

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
    cosh2 = cosh(k{num}.*(hh-zu(2))) ; 
    ucosh1 = uu(:,4) .* cosh1 ; %top most
    ucosh2 = uu(:,2) .* cosh2 ; %2nd from bottom
if contains(categoryname, 'h158')
    uu(:,3) = zeros(height(uu),1) ;  
    uu(:,4) = zeros(height(uu),1) ; 
    ucosh1 = uu(:,2) .* cosh(k{num}.*(hh-zu(2))) ;
    ucosh2 = uu(:,1) .* cosh2 ;
end

ucosh1rms = rms(ucosh1) ; 
ucosh2rms = rms(ucosh2) ; 
ucoshratio = ucosh1rms/ucosh2rms ; 
coshratio = cosh1./cosh2 ; 

kkhh = round( (k{num} .* hh) / kfactor) ; % value for color map

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
% if contains(categoryname, "Base")
%     g1 = scatter(coshratio,ucoshratio, 55, khcolors(kkhh,:), "^", 'filled', 'MarkerEdgeColor', "black") ; hold on
% elseif contains(categoryname, "High")
%     g2 = scatter(coshratio,ucoshratio, 55, khcolors(kkhh,:), "o", 'filled', 'MarkerEdgeColor', "black") ; hold on
% elseif contains(categoryname, "Low")
%     g3 = scatter(coshratio,ucoshratio, 55, khcolors(kkhh,:), "square", 'filled', 'MarkerEdgeColor', "black") ; hold on
% end
if contains(categoryname, "Base")
    g1 = scatter(coshratio,ucoshratio, 55, catcolors(1,:), "^", 'filled', 'MarkerEdgeColor', "black") ; hold on
elseif contains(categoryname, "High")
    g2 = scatter(coshratio,ucoshratio, 55, catcolors(2,:), "o", 'filled', 'MarkerEdgeColor', "black") ; hold on
elseif contains(categoryname, "Low")
    g3 = scatter(coshratio,ucoshratio, 55, catcolors(3,:), "square", 'filled', 'MarkerEdgeColor', "black") ; hold on
end


xlim([0,3])
xlabel('$cosh_{ADV top}/cosh_{ADV 2nd bottom}$')
ylim([0,4.5])
ylabel('$u_{ADV top,rms}/u_{ADV 2nd bottom,rms}$')
% colormap(gca, winter)
% colorbar = colorbar; 
% colorbar.Label.String = '$k*hv$';
% colorbar.Label.FontSize = 12 ;
% colorbar.Label.Interpreter = 'latex' ;
% clim([kmapmin kmapmax])
title("Cosh Dependency")
% set(gcf, 'Units','normalized','OuterPosition',[ 0 0 1 1]) ; 
set(gcf, 'Position', [100, 100, 1200, 800]);

end 

end

xxx = linspace(0,10,10) ; 
yyy = linspace(0,10,10) ; 

figure(1)
g4 = plot(xxx, yyy, 'black') ; hold on
% legend([g4] , {"1:1"})
legend([g4, g1, g2, g3] , {"1:1", "Baseline", "High Density", "Low Density"})

savfigname = "CoshDependency2ndBotCompColor2-ZOOM.png" ;
  saveas(gcf, fullfile(savfolderpath, savfigname)) ; %close all