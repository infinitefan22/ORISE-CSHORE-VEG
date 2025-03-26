clear ; clc ; close all
addpath('./ClaraFunctions') ; 
addpath('./data') ; 
addpath('./mfiles') ; 
 %% Design elements and ease of use
% plottingfig1 = 0 ;
% plottingfig2 = 0 ; 
% savefigures = 0 ; 
% savfolderpath = '/home/elizabeth/Desktop/cshorex-main/osu_mangrove/ClaraFigures/EP_KeltyPlotsMar192025/'; 

% set(0, 'DefaultAxesFontName', 'Nimbus Roman', 'DefaultTextFontName', 'Nimbus Roman')
 set(0,'defaultTextInterpreter','latex')
 set(groot,'defaultAxesTickLabelInterpreter','latex') 
 set(groot, 'defaultLegendInterpreter','latex')
%% Setting Constants and Loading Data 
cnt = 1 ; 
yw = [-1.425, -1.428, -1.428, -1.433] ; %location of AVGs along width of wave flume (w2-5)
% load('aalldata_Mar062025DELETE.mat') ;
load('aalldata_Mar242025.mat') ;
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
%% Start of Trial Indp Sections
for num = 1:NumofTrials
    trialnum = num ; titlename = join([strrep(categoryname, '_', '-'), ' Trial ', string(trialnum)], "") ;
    w = w{num} ; 
%% Calculations
    for ii = 1:width(w)
        vrms(ii) = mean(sqrt(w(:,ii)^2)) ;
%% Plotting
if plotfig1 == 1 
    plot(w{num}(1))
end

end