%  clear ; clc ; close all ; 
addpath('./ClaraFunctions') ; 
addpath('./data') ; 
addpath('./mfiles') ; 
    load('aalldata_Mar282025.mat') ;
%% Saving Data
savfig1name = 'CdCorrelationPlotsomezoom.png' ; 
savfolderpath = '/home/elizabeth/Desktop/cshorex-main/osu_mangrove/ClaraFigures/EP_CdPlotsApr082025/'; 
%      saveas(gcf, fullfile(savfolderpath, savfig1name)) 
%% Setting Constants and Loading Data 
clear fieldnames ; fieldnames = fieldnames(aalldata) ;
clear graph legendnames
max=100 ; %x and y limits on graph
perfectcoor = linspace(0,max,10) ; 
graphcolors = lines(length(fieldnames)-4) ; 
legendnames = strrep(strrep(fieldnames, '_',' '), 'NoWall', '') ; 
legendnames{end+1} = 'Perfect Correlation' ; 
%% total category names loop
for totalnum = 5:length(fieldnames) %not including baseline sections
%% Setting Variables
categoryname = fieldnames{totalnum} ; %'HighDensity_h270_hv182_NoWall' ;
%    alpha = aalldata.(categoryname).alpha ; 
   Cdexact2 = aalldata.(categoryname).Cdexact2 ; 
   CdKelty = aalldata.(categoryname).CdKelty ; 
%    d = aalldata.(categoryname).d ; 
%    datEf = aalldata.(categoryname).datEf ; 
%    dateta = aalldata.(categoryname).dateta ;
%    datHrms = aalldata.(categoryname).datHrms ; 
%    eta = aalldata.(categoryname).eta ; 
%    eta_init = aalldata.(categoryname).eta_init ; 
%    eta_p = aalldata.(categoryname).eta_p ; 
%    eta0a = aalldata.(categoryname).eta0a ; 
%    eta0b = aalldata.(categoryname).eta0b ; 
%    F2 = aalldata.(categoryname).F2 ; 
%    F2overCd = aalldata.(categoryname).F2overCd ; 
%    Hrmsi = aalldata.(categoryname).Hrmsi ; 
%    hv = aalldata.(categoryname).hv ; 
%    KC = aalldata.(categoryname).KC ; 
%    modeleta = aalldata.(categoryname).modeleta ;
%    modelHrms = aalldata.(categoryname).modelHrms ;
%    p = aalldata.(categoryname).p ; 
%    p_init = aalldata.(categoryname).p_init ; 
%    Re = aalldata.(categoryname).Re ; 
%    sav = aalldata.(categoryname).sav ; 
%    stats = aalldata.(categoryname).stats ; 
%    t = aalldata.(categoryname).t ; 
%    w = aalldata.(categoryname).w ; 
%    u = aalldata.(categoryname).u ; 
%    udum = aalldata.(categoryname).udum ; 
%    waveperiod = aalldata.(categoryname).waveperiod ; 
%    xi = aalldata.(categoryname).xi ; 
%    xp = aalldata.(categoryname).xp ; 
%    xwg = aalldata.(categoryname).xwg ; 
%    zw = aalldata.(categoryname).zw ; 
if contains(categoryname, "HighDensity")|| contains(categoryname, "LowDensity")
    figure(1)
%     for gnum = 1:length(graphcolors)
    for onum = 1:length(CdKelty)
    graph(totalnum-4) = scatter(CdKelty{onum}, Cdexact2{onum}, 100, graphcolors(totalnum-4,:), 'filled') ; hold on %graph(totalnum-4) = 
    end
%     end
    xlim([0,20])
    ylim([0,50])
    xlabel('Kelty Cd','fontsize',16)
    ylabel('Johnson Cd','fontsize',16)
%     legend(legendnames(5:end)) ;  
    title('Cd Calculation Correlation')
end
end

figure(1)
graph(end+1) = plot(perfectcoor, perfectcoor, 'LineWidth', 2) ; hold on
    legend([graph], legendnames(5:end))
    