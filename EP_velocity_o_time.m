%% Explaination of Some Things
% I'm using the fLOESS function from MATLAB FileExchange. A data set will
% be considered noisy if there is a velocity change of >.5m/s within a .25s period. 
clearvars -except aalldata ; clc ; close all ; 
addpath('./ClaraFunctions') ; 
addpath('./ClaraFunctions/LOESS regression smoothing-2.1.0.0') ;
addpath('./data') ; 
addpath('./mfiles') ; 
     if ~exist('aalldata', "var") ;  load('aalldata_20250421.mat') ; end
     
%  savfigname = 'NewADVVelocityPlotsAllLayoutsWide.png' ; 
% savfolderpath = '/home/elizabeth/Desktop/cshorex-main/osu_mangrove/ClaraFigures/Velocity_over_Time_20250414/'; 
currentfolder = pwd ; %MATLABONLINE
savfolderpath = join([currentfolder, '/ClaraFigures/20250424_velocity-o-time/'], '') ; 
%        saveas(gcf, fullfile(savfolderpath, savfigname)) 
% set(0, 'DefaultAxesFontName', 'Nimbus Roman', 'DefaultTextFontName', 'Nimbus Roman')
 set(0,'defaultTextInterpreter','latex')
 set(groot,'defaultAxesTickLabelInterpreter','latex') 
 set(groot, 'defaultLegendInterpreter','latex')
 %% Constants
 fLOESStime = .25 ; %time range to check for noisiness (within t time r change occurs
 fLOESSrange = 'sd' ; % put the exact number or 'sd' for sd
 cnt = 1;
 clear fieldnames ; fieldnames = fieldnames(aalldata) ;
 ADVlabel = {'ADV 2, $z=1.404$m', 'ADV 3, $z=1.550$m','ADV 4, $z=1.720$m','ADV 5, $z=1.858$m'} ; 
zu = [1.4040, 1.5500, 1.7200, 1.8580] ; % ADV 2,3,4,5 (order is correct and checked)
smoothvtrials = {} ; 
for totalnum = 5:length(fieldnames)
    categoryname = fieldnames{totalnum} ; % 'HighDensity_h270_hv182_NoWall' ; %
    set_category_variables
%% Length Dep Constants
NumofTrials = length(t) ; 
%% Start of Trial Indp Sections
for num = [1] %1:NumofTrials  %
    if num == NumofTrials && NumofTrials <10 
        num4tt = join(['0',string(NumofTrials)],'') ; 
    elseif num == 1
        num4tt = '01' ;
    else
         num4tt = string(NumofTrials) ; 
    end
%     savfigname = join([strrep(categoryname,'_','-'),'-Trial',num4tt,'-OG.PNG']) ; 
    trialnum = num ; titlename_trial = join([strrep(categoryname, '_', '-'), ' Trial ', string(trialnum)], "") ;
    titlename = join([strrep(categoryname, '_', '-'),' Trial ',string(num)], "") ;
    uu = u{num} ; 
    tt = t{num} ; 
if contains(categoryname, 'h158')
    uu(:,3) = zeros(height(uu),1) ;  
    uu(:,4) = zeros(height(uu),1) ; 
end

clear savfigname ; savfigname = 'vot smoothed ZOOM.png' ; savfigname = join([categoryname, 'Trial',num4tt,savfigname],'_') ; 

%% fLOESS Function (from MATLAB FileExchange)

fLtime = fLOESStime/(tt(2)-tt(1)) ; 
if contains(fLOESSrange, 'sd')
    fLrange = std(uu) ; 
else 
    fLrange = fLOESSrange  ; 
end
noisecheck = isNoisy(uu, fLrange*2, fLtime) ; 



for wnum = 1:4
    if noisecheck(wnum) ==1 
        lowpass = median(uu(:,wnum)) - std(uu(:,wnum))*1.25 ; 
        highpass = median(uu(:,wnum)) + std(uu(:,wnum))*1.25 ; 
        uuidx = uu(:,wnum) < highpass & uu(:,wnum) > lowpass; %corresponding indices are set to 1

        for jj = 1:length(uuidx) ;
            if uuidx(jj) == 0 ; 
                uu(jj,wnum) = NaN ; 
            end
        end
        uu(:, wnum) = fLOESS([transpose(tt),uu(:, wnum)], .03) ; %64/length(uu(:, wnum))
        wwnum = wnum+1 ; disp(titlename + " ADV " + wwnum + " has been smoothed")
        if ~isempty(smoothvtrials); smoothvtrials{end+1} = titlename; 
        else ; smoothvtrials{1} = titlename ; 
        end
    end
end

%% plotting
figure(cnt) ; cnt=cnt+1;
plot(tt,uu, 'LineWidth', 2)
legend(ADVlabel)
xlim([5,10])
xlabel('$t[s]$')
ylabel('$u [m/s]$')
title(titlename)
% set(gcf, 'Units','normalized','OuterPosition',[ 0 0 1 1]) ; 
set(gcf, 'Position', [100, 100, 1200, 800]);

 saveas(gcf, fullfile(savfolderpath, savfigname)) ; %close all
end 



end
