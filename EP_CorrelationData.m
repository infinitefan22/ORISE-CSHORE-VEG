clearvars -except aalldata; clc ; close all ; 
addpath('./ClaraFunctions') ; 
addpath('./data') ; 
addpath('./mfiles') ; 
     if ~exist('aalldata', "var") ;  load('aalldata_20250428.mat') ; end

currentfolder = pwd ; %MATLABONLINE
savfolderpath = join([currentfolder, '/data/'], '') ; 

corrdata = struct() ; 

clear fieldnames ; fieldnames = fieldnames(aalldata) ;

zu = [1.4040, 1.5500, 1.7200, 1.8580] ; % ADV 2,3,4,5 (order is correct and checked) ; from the bottom of the wave flume, not the false floor
ffh = .85 ; %false floor height

for totalnum = 1:length(fieldnames)
    categoryname = fieldnames{totalnum} ; % 'HighDensity_h270_hv182_NoWall' ; %
    set_category_variables
    NumofTrials = length(t) ; 
    
%% EP_Cosh Data
    for num = 1:NumofTrials  %
        clear hh ucosh1 ucosh2 ucosh1rums ucosh2rms 
        uu = u{num} ; 
        tt = t{num} ; 
        hh = hv{num} ; % still water height
        % cosh1 = cosh(k{num}.*(hh-zu(4))) ; 
        % cosh2 = cosh(k{num}.*(hh-zu(1))) ; 
        cosh1 = cosh(k{num}.*(zu(4)-ffh)) ; 
        cosh2 = cosh(k{num}.*(zu(1)-ffh)) ; 
        ucosh1 = uu(:,4) .* cosh1 ;
        ucosh2 = uu(:,1) .* cosh2 ; 
    if contains(categoryname, 'h158')
        uu(:,3) = zeros(height(uu),1) ;  
        uu(:,4) = zeros(height(uu),1) ; 
        ucosh1 = uu(:,2) .* cosh(k{num}.*(hh-zu(2))) ;
    end
    
    ucosh1rms = rms(ucosh1) ; 
    ucosh2rms = rms(ucosh2) ; 
    ucoshratio = ucosh1rms/ucosh2rms ; 
    coshratio = cosh1./cosh2 ; 
    
    corrdata = structure_variables(corrdata, categoryname, 'ucoshratio',ucoshratio) ; % variables/values stored for correlation analysis
    corrdata = structure_variables(corrdata, categoryname, 'coshratio',coshratio) ; 
    
    end 

%% EP_eta_v_u Data
    for num = 1:NumofTrials  %
        clear uu tt hh eta_pp contsh etaplot
    
        hh = hv{num} ; % still water height at vegetation (in relation to false floor)
        uu = u{num} ; 
        tt = t{num} ; 
        Hrmsii = rms(Hrmsi{num}) ; 
        eta_pp = eta_p{num} ; % we are only interested in pressure at ADVS (press3 data)
        kh = k{num}*hh ; 
        hhv = Hrmsii/hh ; 
     
        etastore = eta_pp(:,3) ; 
        
        if contains(categoryname, 'h158') 
            subpnum = 2 ; 
        elseif contains(categoryname, 'h188')
            subpnum = 3 ; 
        else 
            subpnum = 4 ; 
        end
    
        for ADVnum = 1:subpnum
            contsh(ADVnum) = cosh(k{num}*(hh+zu(ADVnum))) / sinh(k{num}*hh) ; 
            etaplot(:, ADVnum) = etastore * contsh(ADVnum) ; 
        end

        corrdata = structure_variables(corrdata, categoryname, 'u',uu) ; 
        corrdata = structure_variables(corrdata, categoryname, 't',tt) ; 
        corrdata = structure_variables(corrdata, categoryname, 'Hrmsii',Hrmsii) ; 
        corrdata = structure_variables(corrdata, categoryname, 'eta_pp',etastore) ; %rn eta_pp is for 4 adv, we only need the 3, at the location of the ADV stack
        corrdata = structure_variables(corrdata, categoryname, 'contsh',contsh) ;
        corrdata = structure_variables(corrdata, categoryname, 'etacosh',etaplot) ;
        corrdata = structure_variables(corrdata, categoryname, 'kh',kh) ; 
    end 

end   

%save([savfolderpath, 'corrdata_20250529'], 'corrdata')