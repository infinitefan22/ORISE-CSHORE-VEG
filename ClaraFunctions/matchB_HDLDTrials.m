function [fbaseline, tnumb] = matchB_HDLDTrials(aalldata, categoryname, num, Tprange, Hrmsirange) 
clear tnumb fbaseline ogwavetype ogTp ogHrmsi fbwavetype fbTp fbHrmsi 
% Tprange = .5
% Hrmsirange = .05
if contains(categoryname, 'h158')
    fbaseline = 'Baseline_h158_hv073_NoWall' ; 
elseif contains(categoryname, 'h188')
    fbaseline = 'Baseline_h188_hv103_NoWall' ; 
elseif contains(categoryname, 'h233')
    fbaseline = 'Baseline_h233_hv148_NoWall' ; 
elseif contains(categoryname, 'h270')
    fbaseline = 'Baseline_h270_hv185_NoWall' ; 
else
    disp("invalid categoryname")
end

tnumb = [] ; 

ogwavetype = aalldata.(categoryname).wavetype{num} ; 
ogTp = aalldata.(categoryname).Tp{num} ; 
ogHrmsi = rms(aalldata.(categoryname).Hrmsi{num}) ; 

fbwavetype = aalldata.(fbaseline).wavetype ; 
fbTp = aalldata.(fbaseline).Tp ; 
        for ii = 1:length(aalldata.(fbaseline).Hrmsi)
fbHrmsi{ii} = rms(aalldata.(fbaseline).Hrmsi{ii}) ; %this is inside for loop
        end


for jj = 1:length(aalldata.(fbaseline).Tp)
if contains(fbwavetype{jj},ogwavetype)
    if ogTp < fbTp{jj} + Tprange && ogTp > fbTp{jj} - Tprange
        if ogHrmsi < fbHrmsi{jj} + Hrmsirange && ogHrmsi > fbHrmsi{jj} - Hrmsirange
            if isempty(tnumb)
                tnumb(1) = jj ; 
            else 
                tnumb(end+1) = jj ; 
            end
        end
    end
end
end 

% if isempty(tnumb)
%     disp("no Baseline trial matches could be found for " + categoryname + " Trial " + num + ". Factors considered: wavetype, Tp, Hrmsi")
% end
 clear weighteddif ; weighteddif = [] ; 

if length(tnumb) > 1
    for kk = 1:length(tnumb)
        difTp = abs(ogTp - fbTp{tnumb(kk)}) ; 
        difHrmsi = abs(ogHrmsi - fbHrmsi{tnumb(kk)}) ; 
        if isempty(weighteddif)
            weighteddif(1) = difTp/Tprange + difHrmsi/Hrmsirange ; 
        else 
            weighteddif(end+1) = difTp/Tprange + difHrmsi/Hrmsirange ;  
        end
    end

    [~, idxn] = min(weighteddif) ;
    tnumb = tnumb(idxn) ; 
elseif isscalar(tnumb)
    idxn = 1 ; 
    tnumb = tnumb(idxn) ; 
else 
     disp("no Baseline trial matches could be found for " + categoryname + " Trial " + num + ". Factors considered: wavetype, Tp, Hrmsi")
     tnumb = 0 ; 
end


end %idk if this will end up being a function, but we'll see