clearvars -except aalldata corrdata; clc ; close all ; 
addpath('./ClaraFunctions') ; 
addpath('./data') ; 
addpath('./mfiles') ; 
     if ~exist('aalldata', "var") ;  load('aalldata_20250428.mat') ; end
     if ~exist('corrdata', "var") ;  load('corrdata_20250529.mat') ; end

set(0,'defaultTextInterpreter','latex')
set(groot,'defaultAxesTickLabelInterpreter','latex') 
set(groot, 'defaultLegendInterpreter','latex')

currentfolder = pwd ; %MATLABONLINE
savfolderpath = join([currentfolder, '/ClaraFigures/CorrelationCalculation/20250617/'], '') ; 
cnt = 1;
clear fieldnames ; fieldnames = fieldnames(aalldata) ;
ADVlabel = {'ADV 2, $z=1.404$m', 'ADV 3, $z=1.550$m','ADV 4, $z=1.720$m','ADV 5, $z=1.858$m'} ; 
zu = [1.4040, 1.5500, 1.7200, 1.8580] ; % ADV 2,3,4,5 (order is correct and checked) ; from the bottom of the wave flume, not the false floor
ffh = .85 ; %false floor height

%% Code Spec Variables
Tprange = .5 ; 
Hrmsirange = .05 ; 

% %% Loop
% for totalnum = 5:length(fieldnames) %CATEGORIES LOOP
%     categoryname = fieldnames{totalnum} ; % 'HighDensity_h270_hv182_NoWall' ; %
%     set_category_variables
%     NumofTrials = length(t) ; 
% 
%     for num = 1:NumofTrials 

waterhvtit = {'158','188','270'} ; % 233 is skipped since Hd doesn't have 233 trials
for numw = 1:3
    if numw == 1
        BL = 'Baseline_h158_hv073_NoWall' ;
        HD = 'HighDensity_h158_hv070_NoWall' ; 
        LD = 'LowDensity_h158_hv073_NoWall' ;
    elseif numw == 2
        BL = 'Baseline_h188_hv103_NoWall' ;
        HD = 'HighDensity_h188_hv100_NoWall' ; 
        LD = 'LowDensity_h188_hv103_NoWall' ;
    elseif numw == 3
        BL = 'Baseline_h270_hv185_NoWall' ;
        HD = 'HighDensity_h270_hv182_NoWall' ; 
        LD = 'LowDensity_h270_hv185_NoWall' ;
    end

    comptab = [] ; 
    for ii = 1:length(aalldata.(BL).t)
        for jj = 1:length(aalldata.(HD).t)
            for kk = 1:length(aalldata.(LD).t)
                dataBL = aalldata.(BL).Tp{ii} ;
                dataHD = aalldata.(HD).Tp{jj} ;
                dataLD = aalldata.(LD).Tp{kk} ;
                if abs(dataBL-dataHD) < Tprange && abs(dataLD-dataHD) < Tprange && abs(dataBL-dataLD) < Tprange
                    if abs(dataBL-dataHD) < Hrmsirange && abs(dataLD-dataHD) < Hrmsirange && abs(dataBL-dataLD) < Hrmsirange
                        comptab = [comptab; ii,jj,kk] ; 
                    end
                end
            end
        end
    end

    if isempty(comptab)
        disp(['No Trials Contain matches for Baseline, High Density, and Low Density Trials. ' ...
            'Tprange='] + string(Tprange)+' ; Hrmsirange=' + string(Hrmsirange)) ; 
    else 
        ComparisonTab = array2table(comptab,'VariableNames', ...
            {join(['Baseline ',waterhvtit{numw}],''), join(['High Density ',waterhvtit{numw}],''),join(['Low Density ',waterhvtit{numw}],'')} );
        disp('Matching trials across all 3 arrays:');
        disp(ComparisonTab);
    end

end



    % titmatch = contains(fieldnames, waterhvtit) ; 
    % if sum(titmatch) == 3 
    % 
    % else 
    % end
