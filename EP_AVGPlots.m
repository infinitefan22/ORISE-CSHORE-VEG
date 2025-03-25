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
% load('aalldata_Mar062025DELETE.mat') ;
load('aalldata_Mar242025.mat') ;
%% Setttting Variables
categoryname = 'HighDensity_h270_hv182_NoWall' ;
   Cdexact2 = aalldata.(categoryname).Cdexact2 ; 
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
   KC = aalldata.(categoryname).KC ; 
   modeleta = aalldata.(categoryname).modeleta ;
   modelHrms = aalldata.(categoryname).modelHrms ;
   p = aalldata.(categoryname).p ; 
   p_init = aalldata.(categoryname).p_init ; 
   Re = aalldata.(categoryname).Re ; 
   sav = aalldata.(categoryname).sav ; 
   t = aalldata.(categoryname).t ; 
   udum = aalldata.(categoryname).udum ; 
   waveperiod = aalldata.(categoryname).waveperiod ; 
   xi = aalldata.(categoryname).xi ; 
   xp = aalldata.(categoryname).xp ; 
   xwg = aalldata.(categoryname).xwg ; 