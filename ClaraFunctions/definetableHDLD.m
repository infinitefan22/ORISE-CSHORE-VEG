function [resHD, resHDhigh, resHDlow, resLD, resLDhigh, resLDlow] = definetableHDLD(T, err) 
% T is the table with proper headings
% err is yes/no. If yes, error values will be taken from the table
disp(fieldnames(T)) % so you can see if the titles match the code
CdtableHD = T(strcmp(T.Layout,'HD'),:) ; CdtableLD = T(strcmp(T.Layout,'LD'),:) ; %separating the HD and LD
    yerrHDtop = zeros(height(CdtableHD), 1);   yerrHDbot = zeros(height(CdtableHD), 1);   yerrLDtop = zeros(height(CdtableLD), 1); yerrLDbot = zeros(height(CdtableLD), 1); %defining error array size

if strcmp(err, 'yes') == 1 %

    for i = 1:height(CdtableHD)
        if CdtableHD.Cdhigh(i) == 0 ; yerrHDtop(i) = 0 ;
        else ; yerrHDtop(i) = abs(CdtableHD.CD(i) - CdtableHD.Cdhigh(i)); end
        if CdtableHD.Cdlow(i) == 0 ; yerrHDbot(i) = 0 ;
        else ; yerrHDbot(i) = abs(CdtableHD.CD(i) - CdtableHD.Cdlow(i)); end
    end
    
    for i = 1:height(CdtableLD)
        if CdtableLD.Cdhigh(i) == 0 ; yerrLDtop(i) = 0 ;
        else ; yerrLDtop(i) = abs(CdtableLD.CD(i) - CdtableLD.Cdhigh(i)); end
        if CdtableLD.Cdlow(i) == 0 ; yerrLDbot(i) = 0 ;
        else ; yerrLDbot(i) = abs(CdtableLD.CD(i) - CdtableLD.Cdlow(i)); end
    end

else 
    disp('error values of T were not calculated')
    yerrHDtop = 0 ;
    yerrHDbot = 0 ;
    yerrLDtop = 0 ;
    yerrLdbot = 0 ; 
end
resHD = CdtableHD ;
resHDhigh = yerrHDtop ;
resHDlow = yerrHDbot ;
resLD = CdtableLD ;
resLDhigh = yerrLDtop ;
resLDlow = yerrLDbot ;
