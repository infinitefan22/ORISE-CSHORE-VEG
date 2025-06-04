clearvars -except aalldata corrvalues; clc ; close all ;
addpath('./ClaraFunctions') ; 
addpath('./ClaraFunctions/LOESS regression smoothing-2.1.0.0') ;
addpath('./data') ; 
addpath('./mfiles') ; 
     if ~exist('aalldata', "var") ;  load('aalldata_20250428.mat') ; end

currentfolder = pwd ; %MATLABONLINE
savfolderpath = join([currentfolder, '/ClaraFigures/EP_eta_v_u/20250522/EtaCosh v U zoom/'], '') ; 

set(0,'defaultTextInterpreter','latex')
set(groot,'defaultAxesTickLabelInterpreter','latex') 
set(groot, 'defaultLegendInterpreter','latex')
savefig = 0 ; 
%% Constants
cnt = 1; 
clear fieldnames ; fieldnames = fieldnames(aalldata) ;
ADVlabel = {'ADV 2, $z=1.404$m', 'ADV 3, $z=1.550$m','ADV 4, $z=1.720$m','ADV 5, $z=1.858$m'} ; 
zu = [1.4040, 1.5500, 1.7200, 1.8580] ; % ADV 2,3,4,5 (order is correct and checked) ; from the bottom of the wave flume, not the false floor
ffh = .85 ; %false floor height
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

for totalnum = 1:length(fieldnames)
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
    if contains(categoryname, "Low") 
        titlenameshort = 'Low Density' ; 
    elseif contains(categoryname, "High")
        titlenameshort = 'High Density' ; 
    end

    hh = hv{num} ; % still water height at vegetation (in relation to false floor)
    uu = u{num} ; 
    tt = t{num} ; 
    Hrmsii = rms(Hrmsi{num}) ; 
    eta_pp = eta_p{num} ; % we are only interested in pressure at ADVS (press3 data)
    kh = k{num}*hh ; 
    hhv = Hrmsii/hh ; 
    corrvalues = structure_variables(corrvalues, categoryname, 'u',u) ; % variables/values stored for correlation analysis
    corrvalues = structure_variables(corrvalues, categoryname, 't',t) ; 
    corrvalues = structure_variables(corrvalues, categoryname, 'uu',uu) ; 
    corrvalues = structure_variables(corrvalues, categoryname, 'tt',tt) ; 
    corrvalues = structure_variables(corrvalues, categoryname, 'Hrmsii',Hrmsii) ; 
    etastore = eta_pp(:,3) ; 
    corrvalues = structure_variables(corrvalues, categoryname, 'eta_pp',etastore) ; %rn eta_pp is for 4 adv, we only need the 3, at the location of the ADV stack
    corrvalues = structure_variables(corrvalues, categoryname, 'kh',kh) ; 
    corrvalues = structure_variables(corrvalues, categoryname, 'hhv',hhv) ; 
%% Plotting
    subplotnum2 = [4,3,2,1] ; 
    if contains(categoryname, 'h158') 
        subpnum = 2 ; 
    elseif contains(categoryname, 'h188')
        subpnum = 3 ; 
    else 
        subpnum = 4 ; 
    end

    for ADVnum = 1:subpnum
        contsh = cosh(k{num}*(hh+zu(ADVnum)) / sinh(k{num}*hh)) ; 
        % contsh = cosh((k{num}*(zu(ADVnum)-ffh)) / sinh(k{num}*hh)) ;
        etaplot = eta_pp(:,3) .* contsh ; 

        corrvalues = structure_variables(corrvalues, categoryname, 'eta_pCosh',etaplot) ; %storing etaplot

        figure(cnt)
            sgtitle(join([titlenameshort, ', $kh$=', num2str(kh), ' $H_{rmsi}/hv$=' , num2str(hhv)],''))
            subplot(4,1, subplotnum2(ADVnum))
            plot(tt, uu(:, ADVnum)) ; hold on
            plot(tt, etaplot) ; hold on
            xlim([5, 35])
            ylim([-2,2])
            if ADVnum ==1 ; xlabel("time [s]") ; end
            title(join(["ADV ", string(ADVnum), ', z=', zu(ADVnum), 'm'], ''))
            legend({'$u$ [m/s]', '$\eta*cosh$'}) % u = velocity, eta = water level elevation

    end ; cnt = cnt+1 ; 

        set(gcf, 'Position', [100, 100, 1800, 1200]);
        savfigname = join(['EtaU ',titlename, 'zoom.png'],'') ;
        % pause(3)
        if savefig == 1 ; exportgraphics(gcf, fullfile(savfolderpath, savfigname), 'Resolution', 300); end 
        %if savefig == 1 ; saveas(gcf, fullfile(savfolderpath, savfigname)) ; end %close all


end   

end 