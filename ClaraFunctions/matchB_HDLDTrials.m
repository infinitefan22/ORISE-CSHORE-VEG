function matchB_HDLDTrials(aalldata, categoryname, num)

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

ogTp = aalldata.(categoryname).Tp{num} ; 
ogHrmsi = aalldata.(categoryname).Hrmsi{num} ; 
ogwavetype = aalldata.(categoryname).wavetype{num} ; 

fbTp = aalldata.(fbaseline).Tp ; 
fbHrmsi = aalldata.(fbaseline).Hrmsi ; 
fbwavetype = aalldata.(fbaseline).wavetype ; 

for jj = 1:length(aalldata.(fbaseline).Tp)
if aalldata.(fbaseline).wavetype{jj} == aalldata.(categoryname).wavetype{num} ; 

if idxwavetype == 0
    disp("wavetype match not found")
end

idxTp =  == ogTp ; 


end %idk if this will end up being a function, but we'll see