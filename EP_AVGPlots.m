%  clear ; clc ; close all ; 
addpath('./ClaraFunctions') ; 
addpath('./data') ; 
addpath('./mfiles') ; 
%  load('aalldata_Mar282025.mat') ;
 %% Design elements and ease of use
 savefigures = 0 ; % the figs must be plotted in order to save them
 plotfig1 = 0 ; %indv trial velocity vs time
 plotfig2 = 0 ; % indv layout ave velocity/trial
 plotfig205 = 0 ; % subplot ver of fig 2
 plotfig3 = 0 ; % all layouts ave velocity/layout
 plotfig4 = 0 ; % indv layout ave Re/trial
 plotfig405 = 1 ; % subplot ver of fig 4
 plotfig5 = 1 ; % all layouts ave Re/layout
 
 savfig205name = 'AllTrials_ADVmeanvelocity.png' ; 
 savfig3name = 'AllLayouts_ADVmeanvelocity.png' ; 
 savefig405name = 'AllTrials_ADVrmsRe.png' ;
 savefig5name = 'AllLayouts_ADVrmsRe.png' ;
 savfolderpath = '/home/elizabeth/Desktop/cshorex-main/osu_mangrove/ClaraFigures/EP_ADVPlotsMar312025/'; 

% set(0, 'DefaultAxesFontName', 'Nimbus Roman', 'DefaultTextFontName', 'Nimbus Roman')
 set(0,'defaultTextInterpreter','latex')
 set(groot,'defaultAxesTickLabelInterpreter','latex') 
 set(groot, 'defaultLegendInterpreter','latex')
%% Setting Constants and Loading Data 
cnt = 1 ; % for graphs 1-100
cnt100 = 1 ; % for graphs 100+
yw = [-1.425, -1.428, -1.428, -1.433] ; %location of AVGs along width of wave flume (w2-5)
graphcolors = lines(4) ; %number is the number of ADV gauges
ADVlabel = {'ADV 2, $d=1.404$m', 'ADV 3, $d=1.550$m','ADV 4, $d=1.720$m','ADV 5, $d=1.858$m'} ; 
Daverage = 0.038 ;
clear fieldnames ; fieldnames = fieldnames(aalldata) ;
%% total category names loop
for totalnum = 1:length(fieldnames)
%% Setting Variables
categoryname = fieldnames{totalnum} ; %'HighDensity_h270_hv182_NoWall' ;
   alpha = aalldata.(categoryname).alpha ; 
   Cdexact2 = aalldata.(categoryname).Cdexact2 ; 
   CdKelty = aalldata.(categoryname).CdKelty ; 
   d = aalldata.(categoryname).d ; 
   datEf = aalldata.(categoryname).datEf ; 
   dateta = aalldata.(categoryname).dateta ;
   datHrms = aalldata.(categoryname).datHrms ; 
   eta = aalldata.(categoryname).eta ; 
   eta_init = aalldata.(categoryname).eta_init ; 
   eta_p = aalldata.(categoryname).eta_p ; 
   eta0a = aalldata.(categoryname).eta0a ; 
   eta0b = aalldata.(categoryname).eta0b ; 
   F2 = aalldata.(categoryname).F2 ; 
   F2overCd = aalldata.(categoryname).F2overCd ; 
   Hrmsi = aalldata.(categoryname).Hrmsi ; 
   hv = aalldata.(categoryname).hv ; 
   KC = aalldata.(categoryname).KC ; 
   modeleta = aalldata.(categoryname).modeleta ;
   modelHrms = aalldata.(categoryname).modelHrms ;
   p = aalldata.(categoryname).p ; 
   p_init = aalldata.(categoryname).p_init ; 
%    Re = aalldata.(categoryname).Re ; 
   sav = aalldata.(categoryname).sav ; 
   stats = aalldata.(categoryname).stats ; 
   t = aalldata.(categoryname).t ; 
%    w = aalldata.(categoryname).w ; 
   u = aalldata.(categoryname).u ; 
   udum = aalldata.(categoryname).udum ; 
   waveperiod = aalldata.(categoryname).waveperiod ; 
   xi = aalldata.(categoryname).xi ; 
   xp = aalldata.(categoryname).xp ; 
   xwg = aalldata.(categoryname).xwg ; 
%    zw = aalldata.(categoryname).zw ; 
%% Length Dep Constants
NumofTrials = length(t) ; 
savfig2name = join([categoryname,'_ADVmeanvelocity.png'],'') ;
savfig4name = join([categoryname,'_ADVrmsRe.png'],'') ;

%% Start of Trial Indp Sections
for num = 1:NumofTrials
    trialnum = num ; titlename_trial = join([strrep(categoryname, '_', '-'), ' Trial ', string(trialnum)], "") ;
    titlename = join([strrep(categoryname, '_', '-')], "") ;
    savfig1name = join([categoryname,'_Trial',string(num),'_ADV.png'],'') ; 
    
    uu = u{num} ; 
    tt = t{num} ; 
%% Calculations
    vrms = zeros(1, width(uu)) ; % vrmstrial = zeros(num,width(uu)) ; 
    for ii = 1:width(uu)
        vrms(ii) = rms(uu(:,ii)) ;
    end
    vrmstrial(num, :) = vrms ; 
    vrmscategories.(categoryname) = rms(vrmstrial) ;

    Re = zeros(1,width(vrms)) ; % Retrial = zeros(num,length(vrms)) ; 
    for jj = 1:width(vrms)
        Re(jj) = vrms(jj) .* Daverage / 1002e-6 ; 
    end
    Retrial(num, :) = Re ; 
    Recategories.(categoryname) = rms(Retrial) ; 
%% Plotting
if plotfig1 == 1 % ADV veloccity over time / Trial
    figure(cnt) ; cnt = cnt+1 ; 
        for jnum=1:width(uu)
            plot(tt,uu(:,jnum),'LineWidth',2) ; hold on 
        end
    xlabel('$t[s]$','fontsize',16)
    ylabel('$u$','fontsize',16)
     ylim([min(min(uu))-.2, max(max(uu))+.2])
    title(titlename_trial,'interpreter','latex');
    legend(ADVlabel)
    if savefigures==1 ; saveas(gcf, fullfile(savfolderpath, savfig1name)) ; end
end

end %specific trials loop
%% Outside of Trial Loop 
%%Velocity Plots
if plotfig2 ==1 % ave velocity/trial/layout
    figure(cnt) ;cnt = cnt+1 ; 
    for num = 1:NumofTrials
    for onum = 1:height(graphcolors)
        scatter(num*ones(1,1),vrmstrial(num,onum), 100, graphcolors(onum,:), 'filled') ; hold on 
        for lnum =1:4 % for 4 ADVs
            text(num, vrmstrial(num,onum),string(round(vrmstrial(num,onum),3)), 'VerticalAlignment','bottom', 'HorizontalAlignment','right', 'fontsize', 10)
        end
    end
    end
    legend(vrmstrial(1,:), ADVlabel, 'location','northeast') ; %lgd.Position = [6.5, 3.5,.2, .08] ;%, 'location','northeast')
    xlim([0, NumofTrials+1]) ; xticks(0:1:NumofTrials+1)
    xlabel('Trial Number','fontsize',16)
    ylabel('$u_{rms}$ [m/s]','fontsize',16)
    title(join([titlename, ' ADV Velocity'], ""))
    if savefigures==1 ; saveas(gcf, fullfile(savfolderpath, savfig2name)) ; end
end

if plotfig205 ==1 % subplot ver of fig 2
    figure(100) ;  %can't put cnt+1 here, so need at start of next fig (cause of subplots)
    sgtitle('ADV Velocity') 
    cnum = 0 ; if contains(categoryname, 'LowDensity') ; cnum =1 ;end
    subplot(3,4,totalnum+cnum) ;
    for num = 1:NumofTrials
    for onum = 1:height(graphcolors)
        scatter(num*ones(1,1),vrmstrial(num,onum), 36, graphcolors(onum,:), 'filled') ; hold on 
    end
    end
    if totalnum==length(fieldnames) ; legend(vrmstrial(1,:), ADVlabel, 'location','bestoutside') ; end %lgd.Position = [6.5, 3.5,.2, .08] ;%, 'location','northeast')
    xlim([0, NumofTrials+1]) ; xticks(0:1:NumofTrials+1)
    ylim([0,.8]) ; yticks(0:.2:.8) ; %make sure the max are the same for ylim and yticks
    xlabel('Trial Number') %,'fontsize',8)
    ylabel('$u_{rms}$ [m/s]') %,'fontsize',8)
    title(titlename)
    if savefigures==1 ; saveas(gcf, fullfile(savfolderpath, savfig205name)) ; end
end

if plotfig3 ==1 % ave velocity of all layouts 
    figure(cnt) ; cnt=cnt+1;
    for onum = 1:height(graphcolors)
        scatter(totalnum*ones(1,1),vrmscategories.(categoryname)(onum), 100, graphcolors(onum,:), 'filled') ; hold on 
    end
    legend(vrmstrial(1,:), ADVlabel, 'location','northeast') ; %lgd.Position = [6.5, 3.5,.2, .08] ;%, 'location','northeast')
    xlim([0, length(fieldnames)+1]) ; xticks(0:1:length(fieldnames)+1)
    xticklabels(fieldnames)
    xlabel('Wave Flume Layout','fontsize',16)
    ylabel('$u_{rms}$ [m/s]','fontsize',16)
    title('ADV Mean Velocity per Layout')
    if savefigures==1 ; saveas(gcf, fullfile(savfolderpath, savfig3name)) ; end
end

if plotfig4 ==1 % ave velocity/trial/layout
    figure(cnt) ;cnt = cnt+1 ; 
    for num = 1:NumofTrials
    for onum = 1:height(graphcolors)
        scatter(num*ones(1,1),Retrial(num,onum), 100, graphcolors(onum,:), 'filled') ; hold on 
        for lnum =1:4 % for 4 ADVs
            text(num, Retrial(num,onum),string(round(Retrial(num,onum),3)), 'VerticalAlignment','bottom', 'HorizontalAlignment','right', 'fontsize', 10)
        end
    end
    end
    legend(Retrial(1,:), ADVlabel, 'location','northeast') ; %lgd.Position = [6.5, 3.5,.2, .08] ;%, 'location','northeast')
    xlim([0, NumofTrials+1]) ; xticks(0:1:NumofTrials+1)
    xlabel('Trial Number','fontsize',16)
    ylabel('$Re_{rms}$','fontsize',16)
    title(join([titlename, 'Reynolds Number'], ""))
    if savefigures==1 ; saveas(gcf, fullfile(savfolderpath, savfig4name)) ; end
end

if plotfig405 ==1 % subplot ver of fig 2
    figure(101) ; %can't put cnt+1 here, so need at start of next fig
    sgtitle('Reynolds Number') 
    cnum = 0 ; if contains(categoryname, 'LowDensity') ; cnum =1 ;end
    subplot(3,4,totalnum+cnum) ;
    for num = 1:NumofTrials
    for onum = 1:height(graphcolors)
        scatter(num*ones(1,1),Retrial(num,onum), 36, graphcolors(onum,:), 'filled') ; hold on 
    end
    end
    if totalnum==length(fieldnames) ; legend(Retrial(1,:), ADVlabel, 'location','bestoutside') ; end %lgd.Position = [6.5, 3.5,.2, .08] ;%, 'location','northeast')
    xlim([0, NumofTrials+1]) ; xticks(0:1:NumofTrials+1)
    ylim([0,.8]) ; yticks(0:.2:.8) ; %make sure the max are the same for ylim and yticks
    xlabel('Trial Number') %,'fontsize',8)
    ylabel('$Re_{rms}$') %,'fontsize',8)
    title(titlename)
    if savefigures==1 ; saveas(gcf, fullfile(savfolderpath, savfig405name)) ; end
end

if plotfig5 ==1 % ave velocity of all layouts 
    figure(102) ; 
    for onum = 1:height(graphcolors)
        scatter(totalnum*ones(1,1),Recategories.(categoryname)(onum), 100, graphcolors(onum,:), 'filled') ; hold on 
    end
    legend(Retrial(1,:), ADVlabel, 'location','northeast') ; %lgd.Position = [6.5, 3.5,.2, .08] ;%, 'location','northeast')
    xlim([0, length(fieldnames)+1]) ; xticks(0:1:length(fieldnames)+1)
    xticklabels(fieldnames)
    xlabel('Wave Flume Layout','fontsize',16)
    ylabel('$Re_{rms}$','fontsize',16)
    title('$Re_{rms}$ per Layout')
    if savefigures==1 ; saveas(gcf, fullfile(savfolderpath, savfig5name)) ; end
end

end %total category names loop

