function [CdvalueK, sdvalueK, strJ, strK, Khrmsi, Jhrmsi, Ktp, Jtp] = match2KeltyTable(structog, struct, categoryname, trialnum)

Tprange = .5 ; 
Hrmsirange = .05 ; 
%%
ogwavetype = structog.(categoryname).wavetype{trialnum} ; 

if contains(ogwavetype, 'andom')
    ttable = 'Random' ;
elseif contains(ogwavetype, 'egular')
    ttable ='Regular' ;
end

if contains(categoryname, "High")
    mlayout = 'HD' ; 
elseif contains(categoryname, "Low")
    mlayout = 'LD' ;
elseif contains(categoryname, "Base")
    disp("match2KeltyTable cannot be run for Baseline Trials, please refilter to only include HD or LD trials")
end

if contains(categoryname, "h270")
    mtrial = 'h1' ;
elseif contains(categoryname,"h233")
    mtrial = 'h2' ;
elseif contains(categoryname,"h188")
    mtrial = 'h3' ;
elseif contains(categoryname,"h158")
    mtrial = 'h4' ;
end

tnumb = [] ; 
lll = length([struct.(ttable).Trial]) ;
for fnum = 3:lll
    if contains(struct.(ttable).Layout{fnum}, mlayout) && contains(struct.(ttable).Trial{fnum}, mtrial)
        if isempty(tnumb)
            tnumb(1) = fnum ; 
        else 
            tnumb(end+1) = fnum ; 
        end
    end
end

% this is how far I am editing, everything code wise above should be fine
% with spec exceptions
clear weighteddif ; weighteddif = [] ; 

ogTp = structog.(categoryname).Tp{trialnum} ; 
ogHrmsi = rms(structog.(categoryname).Hsjohnson{trialnum}) ; 
fbTp(:) = struct.(ttable).Tp(tnumb(:)) ; 
fbHrmsi(:) = struct.(ttable).Hm0(tnumb(:)) ; 

newt =[] ;
for jj = 1:length(fbTp)
    if ogTp < fbTp(jj) + Tprange && ogTp > fbTp(jj) - Tprange
        if ogHrmsi < fbHrmsi(jj) + Hrmsirange && ogHrmsi > fbHrmsi(jj) - Hrmsirange
            if isempty(newt)
                newt(1) = jj ; 
            else 
                newt(end+1) = jj ; 
            end
        end
    end
end 

if length(newt) > 1
    disp(categoryname +" Trial "+ trialnum + " has multiple matches. Kelty Table " + ogwavetype + ". Matched Rows: "+tnumb(newt(1)) +", " + tnumb(newt(2)))
    for kk = 1:length(newt)
        difTp = abs(ogTp - fbTp(newt(kk))) ; 
        difHrmsi = abs(ogHrmsi - fbHrmsi(newt(kk))) ; 
        if isempty(weighteddif)
            weighteddif(1) = difTp/Tprange + difHrmsi/Hrmsirange ; 
        else 
            weighteddif(end+1) = difTp/Tprange + difHrmsi/Hrmsirange ;  
        end
    end
    [~, idxn] = min(weighteddif) ;
    finalt = tnumb(newt(idxn)) ; 
    Khrmsi = fbHrmsi(newt(idxn)) ; 
    Ktp = fbTp(newt(idxn)) ;

elseif isscalar(newt)
    idxn = 1 ; 
    finalt = tnumb(newt(idxn)) ; 
    Khrmsi = fbHrmsi(newt(idxn)) ; 
    Ktp = fbTp(newt(idxn)) ; 
else 
     disp("no Kelty Table trial matches could be found for " + categoryname + " Trial " + trialnum + ". Factors considered: wavetype, Tp, Hrmsi")
     finalt = 1 ; 
     Khrmsi = NaN ; 
     Ktp = NaN ; 
end

CdvalueK = struct.(ttable).CD(finalt) ;
sdvalueK = struct.(ttable).sd(finalt) ;
strJ = join([categoryname, " Trial ", trialnum], "") ; 
strK = struct.(ttable).Trial(finalt) ; 
Jhrmsi = ogHrmsi ; 
Jtp = ogTp ; 

end %idk if this will end up being a function, but we'll see