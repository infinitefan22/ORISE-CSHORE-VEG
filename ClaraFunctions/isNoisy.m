% if there is a change greater than the sd of the function within a quarter
% second, then fLOESS will be applied
% data = the data set you want checked for noisiness
% range = the min change you want
% tiem = the min time the change must occur for the data set to be
% considered noisy
function res = isNoisy(data, range, time)

res = zeros(1, width(data)) ; 
check = zeros(1,length(data) - time) ; 

for i = 1:width(data)

for t = 1:length(data) - time
    mxx = max(data(t:t+time, i)) ; 
    mnn = min(data(t:t+time, i)) ; 
    if abs(mxx-mnn) <= range(i)
         check(t) = false ;
    else 
        check(t) = true ;
        t = length(data)-time ; 
    end
end

if any(check == 1 )
    res(i) = 1 ;
else 
    res(i) = 0 ; 
end


end



end % of function