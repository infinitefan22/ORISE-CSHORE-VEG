% clear ; clc ; close all ; 
addpath('./ClaraFunctions') ; 
addpath('./data') ; 
addpath('./mfiles') ; 
 %load('aalldata_Mar242025.mat') ;
 %% Design elements and ease of use
 plotfig1 = 1 ;
 plotfig2 = 1 ; 
 savefigures = 1 ; % the figs must be plotted in order to save them
 savfolderpath = '/home/elizabeth/Desktop/cshorex-main/osu_mangrove/ClaraFigures/EP_AVDPlotsMar272025/'; 

% set(0, 'DefaultAxesFontName', 'Nimbus Roman', 'DefaultTextFontName', 'Nimbus Roman')
 set(0,'defaultTextInterpreter','latex')
 set(groot,'defaultAxesTickLabelInterpreter','latex') 
 set(groot, 'defaultLegendInterpreter','latex')
%% Setting Constants and Loading Data 
cnt = 1 ; 
yw = [-1.425, -1.428, -1.428, -1.433] ; %location of AVGs along width of wave flume (w2-5)
graphcolors = lines(4) ; %number is the number of AVD gauges
AVDlabel = {'AVD 2, $d=1.404$m', 'AVD 3, $d=1.550$m','AVD 4, $d=1.720$m','AVD 5, $d=1.858$m'} ; 
% load('aalldata_Mar062025DELETE.mat') ;

%% Setttting Variables
categoryname = 'HighDensity_h270_hv182_NoWall' ;
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
   Re = aalldata.(categoryname).Re ; 
   sav = aalldata.(categoryname).sav ; 
   stats = aalldata.(categoryname).stats ; 
   t = aalldata.(categoryname).t ; 
   w = aalldata.(categoryname).w ; 
   udum = aalldata.(categoryname).udum ; 
   waveperiod = aalldata.(categoryname).waveperiod ; 
   xi = aalldata.(categoryname).xi ; 
   xp = aalldata.(categoryname).xp ; 
   xwg = aalldata.(categoryname).xwg ; 
   zw = aalldata.(categoryname).zw ; 
%% Length Dep Constants
NumofTrials = length(t) ; 
savfig2name = join([categoryname,'_AVDmeanvelocity.png'],'') ;
%% Start of Trial Indp Sections
for num = 1:NumofTrials
    trialnum = num ; titlename = join([strrep(categoryname, '_', '-'), ' Trial ', string(trialnum)], "") ;
    savfig1name = join([categoryname,'_Trial',string(num),'_ADV.png'],'') ; 
    
    ww = w{num} ; 
    tt = t{num} ; 
%% Calculations
    vrms = zeros(1, width(ww)) ; 
    for ii = 1:width(ww)
        vrms(ii) = mean(sqrt(ww(:,ii).^2)) ;
    end
    vrmstrial(num, :) = vrms ; 
%% Plotting
if plotfig1 == 1 
    figure(cnt) ; cnt = cnt+1 ; 
        for jnum=1:width(ww)
            plot(tt,ww(:,jnum),'LineWidth',2) ; hold on 
        end
    xlabel('$t[s]$','fontsize',16)
    ylabel('$u$','fontsize',16)
     ylim([min(min(ww))-.2, max(max(ww))+.2])
    title(titlename,'interpreter','latex');
    legend(AVDlabel)
    if savefigures==1 ; saveas(gcf, fullfile(savfolderpath, savfig1name)) ; end
end

end
%% Outside of Trial Loop
if plotfig2 ==1
    figure(cnt) ; cnt=cnt+1 ; 
    for num = 1:NumofTrials
    for onum = 1:height(graphcolors)
        scatter(num*ones(1,1),vrmstrial(num,onum), 100, graphcolors(onum,:), 'filled') ; hold on 
        for lnum =1:4 % for 4 AVDs
            text(num, vrmstrial(num,onum),string(round(vrmstrial(num,onum),3)), 'VerticalAlignment','bottom', 'HorizontalAlignment','right', 'fontsize', 10)
        end
    end
    end
    legend(vrmstrial(1,:), AVDlabel, 'location','northeast') ; %lgd.Position = [6.5, 3.5,.2, .08] ;%, 'location','northeast')
    xlim([0, NumofTrials+1]) ; xticks(0:1:NumofTrials+1)
    xlabel('Trial Number','fontsize',16)
    ylabel('$u$ [m/s]','fontsize',16)
    title(join([titlename, ' AVD Velocity'], ""))
    if savefigures==1 ; saveas(gcf, fullfile(savfolderpath, savfig2name)) ; end

end
