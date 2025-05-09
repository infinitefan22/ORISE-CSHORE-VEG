clear ; clc ; close all
addpath('./ClaraFunctions') ; 
addpath('./data') ; 
addpath('./mfiles') ; 
 %% Design elements and ease of use
plottingfig1 = 0 ;
plottingfig2 = 0 ; 
savefigures = 0 ; 
savfolderpath = '/home/elizabeth/Desktop/cshorex-main/osu_mangrove/ClaraFigures/EP_KeltyPlotsMar192025/'; 

% set(0, 'DefaultAxesFontName', 'Nimbus Roman', 'DefaultTextFontName', 'Nimbus Roman')
 set(0,'defaultTextInterpreter','latex')
 set(groot,'defaultAxesTickLabelInterpreter','latex') 
 set(groot, 'defaultLegendInterpreter','latex')
%% Setting Constants and Loading Data 
CNT = 1 ; 
% load('aalldata_Mar062025DELETE.mat') ;
load('aalldata_Mar242025.mat') ;
offset = [3.0 3.9 2.8 3.5 3.8 3.8 ]; 

%% variables (within aalldata)
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

%% length dependant constants
offset = [offset, zeros(1, length(xi) - length(offset))] ;
numofTrials = length(xi) ; 

% eta = ...
% etc/. 
%% for loops 

for num = 1:length(xi)
    %%
    trialnum = num ; titlename = join([strrep(categoryname, '_', '-'), ' Trial ', string(trialnum), ' $H_{rms}$=', sprintf('%2.2f',Hrmsi{num}(1))], "") ;
    savfig1name = join([categoryname,'_Trial',string(num),'_WaveVeg.png'],'') ; 
    savfig2name = join([categoryname,'_Trial',string(num),'_ModvsMeas.png'],'') ; 
    [j1, j2] = min(abs(xi{num}-43));
    j = num ; 

    if plottingfig1 == 1
    figure(CNT);clf;clear hh hl;
         CNT=CNT+1 ; % CLARA
         subplot(2,1,1)
        %plot(xp,sqrt(8)*std(p)/9810,'bs-');hold on
         hh1(1) = plot(xp{num}(1:6),sqrt(8)*std(eta_p{num}(:,1:6)),'bs-','markerfacecolor','b');hold on
        %plot(xuswg,sqrt(8)*std(etaus),'rv','markerfacecolor','r');hold on
         hh1(2) = plot(xwg{num},sqrt(8)*std(eta{num}),'rs-','markerfacecolor','r');hold on
         hh1(3) = plot(xi{num},modelHrms{num},'k-','linewidth',2);hold on ; 
         hh1(3) = plot(xi{num},datHrms{num},'k--','linewidth',2);hold on      %dname2=dname;dname2(dname=='_')='-';
         title(titlename);
%          title(categoryname+'_Trial '+string(trialnum),' ',sprintf('%2.2f',Hrmsi{num}(1)),'interpreter','latex','fontsize',12) ; %,' ',sprintf('%2.2f',stats.Tp)]
        %xlabel('$x[m]$','interpreter','latex','fontsize',16)
         ylabel('$H_{rms}[m]$','fontsize',16)
         a = axis;
         axis manual
         hf=fill([36 54 54 36],[a(3) a(3) a(4) a(4)],[.2 .5 .2]);set(hf,'facealpha',.2) ;
         text(37,a(3)+.1*(a(4)-a(3)),'Veg. Section','fontsize',14) ;
         text(15,a(3)+.1*(a(4)-a(3)),['$C_d = ',sprintf('%2.2f',Cdexact2{num}),'$'],'fontsize',14) ;
    
          subplot(2,1,2);
          plot(xp{num}(1:6),1000*(mean(p{num}(:,1:6))-p_init{num}(1:6))/9810,'bs-','markerfacecolor','b');hold on
          plot(xwg{num},1000*(mean(eta{num})-eta_init{num}),'rs-','markerfacecolor','r');hold on
          %plot(xuswg,1000*(mean(etaus)-etaus_init),'rv','markerfacecolor','r');hold on
          hh(1) = plot(xi{num},1000*modeleta{num},'k-','linewidth',2);
          hh(2) = plot(xi{num},1000*dateta{num},'k--','linewidth',2);
          hh(3) = plot(xi{num},1000*eta0a{num},'g-','linewidth',2);
          %  hh(2) = plot(xi,1000*eta0b,'g--','linewidth',2);
          %hh(2) = plot(xi,1000*etaiwavesbotshearveg,'k-.','linewidth',2);
          %plot(xi,1000*etaiwavesbotshearwrs,'m--','linewidth',2);
          %plot(xi,1000*etaiwavesbotshearveg2p,'m-','linewidth',2);
          %hh(2) = plot(xi,1000*etaiwavesbotshearveg2,'k-','linewidth',3);
          legend(hh,'Modeled  $\eta$, $u$ ','Measured $\eta, u$','{\bf F}= 0','AutoUpdate','off');
%           hl = legend(hh,'Modeled $ \eta, u $','Measued $\eta, u $','${\bf F}= 0$','AutoUpdate','off');
          %hl = legend(hh,'$S_{xx}, \tau_b$','$S_{xx}, \tau_b, {\bf F}$','AutoUpdate','off');
%           set(hl,'interpreter','latex','fontsize',14,'location','northwest')
          dumyl = ylim;
          if dumyl(2)-dumyl(1)<100
            ylim([mean(ylim)-40 mean(ylim)+40])
          end
          xlabel('$x[m]$','fontsize',16)
          ylabel('$\overline{\eta}[mm]$','fontsize',14)
          

%       dname2=categoryname;dname2(categoryname=='_'|categoryname=='/')='-';
%     %   print('-dpng',[dname2(21:end-1),'.png']) 
%       figname = [categoryname, '_Trial ' , trialnum] ;   % should give you  Baseline_h158_hv073_NoWall_Trial01'
%       print('-dpng',[figname,'.png']);
%     %   a = axis;
%     %   hf=fill([36 54 54 36],[a(3) a(3) a(4) a(4)],[.2 .5 .2]);set(hf,'facealpha',.2)
    %   axis(a)
    if savefigures==1 ; saveas(gcf, fullfile(savfolderpath, savfig1name)) ; end
    end
%%
if plottingfig2 ==1
    %%
  figure(CNT);clf
      CNT=CNT+1 ; % CLARA
      subplot(3,1,1)
      plot(t{num},udum{num},'r','linewidth',2);hold all
        plot(t{num},0*udum{num},'k','linewidth',1);hold all
      plot(sav{num}(j2).t+offset(j),mean(sav{num}(j2).u),'b','linewidth',2)
      ylabel('$u$','interpreter','latex','fontsize',16)
      %xlabel('$t[s]$','interpreter','latex','fontsize',16)
      xlim([0 20]);
      title(titlename,'interpreter','latex');
%       title([strrep(dname(21:end),'_','-'),' ',sprintf('%2.2f',Hrmsi(1)),' ',sprintf('%2.2f',stats.Tp)],'interpreter','latex','fontsize',12)
      %print('-dpng',[dname2(21:end-1),'u.png'])
    
      subplot(3,1,2)
      plot(t{num},-d{num},'r','linewidth',2);hold all
      plot(t{num},0*t{num}+mean(-d{num}),'r--','linewidth',1);hold all
      plot(sav{num}(j2).t+offset(j),-sav{num}(j2).Dx,'b','linewidth',2)
      plot(sav{num}(j2).t+offset(j),ones(size(sav{num}(j2).Dx))*mean(-sav{num}(j2).Dx),'b--','linewidth',1)
      ylabel('$-D$','interpreter','latex','fontsize',16)
      %xlabel('$t[s]$','interpreter','latex','fontsize',16)
      xlim([0 20]);
      %title([strrep(dname(21:end),'_','-'),' ',sprintf('%2.2f',Hrmsi(1)),' ',sprintf('%2.2f',stats.Tp)],'interpreter','latex','fontsize',12)
      %print('-dpng',[dname2(21:end-1),'u.png'])
      
      subplot(3,1,3)
      plot(t{num},0*F2{num},'k');hold all
      hf(1) = plot(t{num},F2{num},'r','linewidth',2);hold all
      plot(t{num},mean(F2{num})+0*t{num},'r--');hold all
      hf(2) = plot(sav{num}(j2).t+offset(j),sav{num}(j2).Fx,'b','linewidth',2) ; 
      plot(sav{num}(j2).t+offset(j),sav{num}(j2).meanFx+0*sav{num}(j2).t,'b--','linewidth',1)
      ylabel('$F_x $','interpreter','latex','fontsize',16)
      xlabel('$t[s]$','interpreter','latex','fontsize',16)
      xlim([0 20]);
      %title([strrep(dname(21:end),'_','-'),' ',sprintf('%2.2f',Hrmsi(1)),' ',sprintf('%2.2f',stats.Tp)],'interpreter','latex','fontsize',12)
      hl = legend(hf,'Using Measured Data','Using Modeled Data','AutoUpdate','off');
      set(hl,'interpreter','latex','fontsize',14,'location','northeast') ;
    %  print('-dpng',[dname2(21:end-1),'_u_d_F.png']) EDITED CLARA
     if savefigures==1 ; saveas(gcf, fullfile(savfolderpath, savfig2name)) ; end
end

end
