clearvars -except aalldata ; clc ; close all ; 
addpath('./ClaraFunctions') ; 
addpath('./ClaraFunctions/LOESS regression smoothing-2.1.0.0') ;
addpath('./data') ; 
addpath('./mfiles') ; 
     if ~exist('aalldata', "var") ;  load('aalldata_20250428.mat') ; end

currentfolder = pwd ; %MATLABONLINE
savfolderpath = join([currentfolder, '/ClaraFigures/EP_eta_v_u/20250506/Zoom 1 (5<x<35)/'], '') ; 

set(0,'defaultTextInterpreter','latex')
set(groot,'defaultAxesTickLabelInterpreter','latex') 
set(groot, 'defaultLegendInterpreter','latex')
savefig = 1 ; 
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

for totalnum = 5:length(fieldnames)
    categoryname = fieldnames{totalnum} ; % 'HighDensity_h270_hv182_NoWall' ; %
    set_category_variables
%% Length Dep Constants
NumofTrials = length(t) ; 

%% Start of Trial Indp Sections
for num = 1:NumofTrials  %
        clear uu tt hh eta_pp contsh etaplot
    if num == NumofTrials && NumofTrials <10 
        num4tt = join(['0',string(NumofTrials)],'') ; 
    elseif num == 1
        num4tt = '01' ;
    else
         num4tt = string(NumofTrials) ; 
    end
    trialnum = num ; titlename_trial = join([strrep(categoryname, '_', '-'), ' Trial ', string(trialnum)], "") ;
    titlename = join([strrep(categoryname, '_', '-'),' Trial ',string(num)], "") ;
    hh = hv{num} ; % still water height
    uu = u{num} ; 
    tt = t{num} ; 
    eta_pp = eta_p{num} ; 


%% Plotting
    if contains(categoryname, 'h158') 
        for ADVnum = 1:2 %ignoring the ADVs that do not collect v consistently at that water height
        contsh = cosh(k{num}*(hh+zu(ADVnum)) / sinh(k{num}*hh)) ; 
        etaplot = eta_pp .* contsh ; 
    
        figure(cnt)
            sgtitle(join([titlename, ' $\eta$ vs. $u$'],''))
            subplot(4,1, ADVnum)
            plot(tt, uu(:, ADVnum)) ; hold on
            plot(tt, etaplot(:, ADVnum)) ; hold on
            xlim([5, 35])
            title(join(["ADV ", string(ADVnum)], ''))
            legend({'Velocity', '$\eta$'})
    
        end ; cnt = cnt+1 ; 
        set(gcf, 'Position', [100, 100, 1200, 800]);
        savfigname = join(['EtaUComp ',titlename, 'zoom.png'],'') ;
        if savefig == 1 ; saveas(gcf, fullfile(savfolderpath, savfigname)) ; end %close all
    else
        for ADVnum = 1:4
        contsh = cosh(k{num}*(hh+zu(ADVnum)) / sinh(k{num}*hh)) ; 
        etaplot = eta_pp .* contsh ; 
    
        figure(cnt)
            sgtitle(join([titlename, ' $\eta$ vs. $u$'],''))
            subplot(4,1, ADVnum)
            plot(tt, uu(:, ADVnum)) ; hold on
            plot(tt, etaplot(:, ADVnum)) ; hold on
            xlim([5, 35])
            title(join(["ADV ", string(ADVnum)], ''))
            legend({'Velocity', '$\eta$'})

        end ; cnt = cnt+1 ; 
        set(gcf, 'Position', [100, 100, 1200, 800]);
        savfigname = join(['EtaUComp ',titlename, 'zoom.png'],'') ;
        if savefig == 1 ; saveas(gcf, fullfile(savfolderpath, savfigname)) ; end %close all
    end


end ; close all ; 

end