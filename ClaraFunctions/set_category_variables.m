clear alpha Cdexct2 CdKelty d datEf datHrms eta eta_init eta_p eta0a eta0b F2 F2overCd Hrmsi hv 
clear k KC modeleta modelHrmsp p_init sav stats t Tp u udum waveperiod xi xp xwg
  
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
   % F2overCd = aalldata.(categoryname).F2overCd ; 
   Hrmsi = aalldata.(categoryname).Hrmsi ; 
   hv = aalldata.(categoryname).hv ; 
   k = aalldata.(categoryname).k ; 
   KC = aalldata.(categoryname).KC ; 
   modeleta = aalldata.(categoryname).modeleta ;
   modelHrms = aalldata.(categoryname).modelHrms ;
   p = aalldata.(categoryname).p ; 
   p_init = aalldata.(categoryname).p_init ; 
%    Re = aalldata.(categoryname).Re ; 
   sav = aalldata.(categoryname).sav ; 
   stats = aalldata.(categoryname).stats ; 
   t = aalldata.(categoryname).t ; 
   Tp = aalldata.(categoryname).Tp ;
%    w = aalldata.(categoryname).w ; 
   u = aalldata.(categoryname).u ; 
   % udum = aalldata.(categoryname).udum ; 
   waveperiod = aalldata.(categoryname).waveperiod ; 
   xi = aalldata.(categoryname).xi ; 
   xp = aalldata.(categoryname).xp ; 
   xwg = aalldata.(categoryname).xwg ; 
%    zw = aalldata.(categoryname).zw ;  

if exist('corrdata', "var") 
    clear c ucoshratio coshratio Hrmsi eta_pp contsh etacosh kh
    c.ucoshratio = corrdata.(categoryname).ucoshratio ;
    c.coshratio = corrdata.(categoryname).coshratio ;
    c.u = corrdata.(categoryname).u ;
    c.t = corrdata.(categoryname).t ;
    c.Hrmsii = corrdata.(categoryname).Hrmsii ;
    c.eta_pp = corrdata.(categoryname).eta_pp ;
    c.contsh = corrdata.(categoryname).contsh ;
    c.etacosh = corrdata.(categoryname).etacosh ;
    c.kh = corrdata.(categoryname).kh ;
end


% if exist('corrdata', "var") 
%     clear c ucoshratio coshratio Hrmsi eta_pp contsh etacosh kh
%     c.ucoshratio = cell2mat(corrdata.(categoryname).ucoshratio) ;
%     c.coshratio = cell2mat(corrdata.(categoryname).coshratio) ;
%     c.Hrmsii = cell2mat(corrdata.(categoryname).Hrmsii) ;
%     c.contsh = cell2mat(corrdata.(categoryname).contsh) ;
%     c.kh = cell2mat(corrdata.(categoryname).kh) ;
%     c.u = corrdata.(categoryname).u ;
%     c.t = corrdata.(categoryname).t ;
%     c.eta_pp = corrdata.(categoryname).eta_pp ;
%     c.etacosh = corrdata.(categoryname).etacosh ;
% end

% c.ucoshratio = cell2mat(corrdata.(categoryname).ucoshratio) ;
%     c.coshratio = cell2mat(corrdata.(categoryname).coshratio) ;
%     c.u = corrdata.(categoryname).u ;
%     c.t = corrdata.(categoryname).t ;
%     c.Hrmsii = corrdata.(categoryname).Hrmsii ;
%     c.eta_pp = corrdata.(categoryname).eta_pp ;
%     c.contsh = corrdata.(categoryname).contsh ;
%     c.etacosh = corrdata.(categoryname).etacosh ;
%     c.kh = corrdata.(categoryname).kh ;
   % c.ucoshratio = cell2mat(corrdata.(categoryname).ucoshratio) ;
   %  c.coshratio = cell2mat(corrdata.(categoryname).coshratio) ;
   %  % c.u = cell2mat(corrdata.(categoryname).u) ;
   %  % c.t = cell2mat(corrdata.(categoryname).t) ;
   %  c.Hrmsii = cell2mat(corrdata.(categoryname).Hrmsii) ;
   %   c.eta_pp = cell2mat(corrdata.(categoryname).eta_pp) ;
   %  c.contsh = cell2mat(corrdata.(categoryname).contsh) ;
   %  % c.etacosh = cell2mat(corrdata.(categoryname).etacosh) ;
   %  c.kh = cell2mat(corrdata.(categoryname).kh) ;