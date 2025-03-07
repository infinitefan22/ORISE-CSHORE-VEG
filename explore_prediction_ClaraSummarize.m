% Clara copied explore_prediction_KeltyData and modified it to run for all
% files created through ClaraSummarize 
% last edit: 2025/02/19
clear ; clc ; close all
plottingog = 0 ; plottingyn = 0 ; plottingpost = 0 ; %1 means plots, 0 means no plots
%og section  % plotting error    %plotting stored data
cnt = 0; %figure numbers
addpath('./ClaraFunctions') ; 
addpath('./data') ; 
%% constants
cf = .05;
icorrect =1;
Dtrunk = 0.1143 ; Droot = 0.0286 ; Daverage = 0.038 ;%both in meters, diameter of trunk and roots of model from Kelty full paper, average is ave of effective diameters p52
fluiddensity = 1000 ; fluidviscosity = .001 ; %kg/m3, kg/ms
offset = [3.0 3.9 2.8 3.5 3.8 3.8 ]; %this line is why j can only go to 6
%% accessing files
% basenameog = {'HighDensity_h270_hv182_NoWall' ,'HighDensity_h270_hv182_Wall'} ; % original basename
basename = '/home/elizabeth/Desktop/cshorex-main/osu_mangrove/data/SummaryFiles/' ; 
% dnames = dir(['./data/',basename,'/T*']);
dnames = dir([basename, '*Summary.mat']) ; 
bnames = dir([basename, 'Baseline*']) ; 
hdnames = dir([basename, 'HighDensity*']) ; 
ldnames = dir([basename, 'LowDensity*']) ; 
%% data processing 1
for j = 1:length(dnames)
  cnt = cnt+1;
  offset = [offset, zeros(1, length(dnames) - length(offset))] ; % making offset the same length as dnames
  clear p u eta eta_p ubp wbp
  dname = ['/home/elizabeth/Desktop/cshorex-main/osu_mangrove/data/SummaryFiles/', dnames(1).name] ; 
  dname = [basename,dnames(j).name,'/']; % idk why I can't load this
%   dname = ['./data/',basename,'/',dnames(j).name,'/']; %makes a string that forms the path to a specific subdirectory in the main data directory of ~mangrove
  load(dnames(j).name)
%   load([dname,'summary.mat']) %accessing the specific files (1-6)
  hv = str2num(dname(strfind(dname,'hv')+2:strfind(dname,'hv')+4))/100+.03; %finding hv in flies, calculating it, changing to number from string
  p = [dat.press.press];%gets press field of data, makes into set list
  u = [dat.u.u];u = u(:,2:5);%same as line 22 but for velocity
  xu = [dat.u.x];xu = xu(:,2:5);%velocity in x dir
  zu = [dat.u.z];zu = zu(:,2:5);%velocity in z dir
  w = [dat.w.w];w = w(:,2:5);
  zw = [dat.w.z];zw = zw(:,2:5);
  xp = [dat.press.x];%access x field of press field of the data

  eta = [dat.wg.eta]; % water level elevation (where the water line is)
  xwg = [dat.wg.x];
  etaus =[dat.uswg.eta];
  xuswg = [dat.uswg.x];

  activity = mean(std(eta),2); %average of standard deviation of each row
  startinds = find(mean(abs(eta(:,find(xwg>50))),2)>activity);startind = startinds(10);
  endinds = find(mean(abs(eta(:,find(xwg<50))),2)>activity);endind = endinds(end-10);
  p_init = mean(p(1:round(startind/2),:));
  eta_init = mean(eta(1:round(startind/2),:)); %average of eta rows 1 to startind/2 rounded to nearest integer and all columns
  etaus_init = mean(etaus(1:round(startind/2),:));%see line 39
  eta = eta(startind:endind,:);
  t = [0:size(eta,1)-1]./100;%element wise division of array of 1st col of eta
  stats = find_stats(t,eta(:,2),4);%function find_stats
  [k,n,c] = dispersion (2*pi/stats.Tp,hv);
  etaus = etaus(startind:endind,:);
  p = p(startind:endind,:);
  u = u(startind:endind,:);
  w = w(startind:endind,:);
  for jj =1:size(u,2)
    ubp(:,jj) = bandpass(u(:,jj),1/100,.5*1/stats.Tp,2*1/stats.Tp);%filters signals so only those within the specified range can pass
    wbp(:,jj) = bandpass(w(:,jj),1/100,.5*1/stats.Tp,2*1/stats.Tp);
  end

  for i = 1:size(p,2)
    if dat.press(i).z-dat.press(i).swd<0;
      eta_p(:,i) = p2eta([p(:,i)-mean(p(:,i))],1/100,hv,dat.press(i).z-dat.press(i).swd);
    else
      eta_p(:,i) = dat.press(i).z+p(:,i)/9810;
    end
  end
  dx = .1;
  xi = dat.press(1).x:dx:dat.press(6).x;
  Hrmsi = interp1(xwg,sqrt(8)*std(eta),xi);
  %Hrmsi = interp1(xp,sqrt(8)*std(eta_p),xi);
  [Sxy Sxx Syy] = radstress (Hrmsi, 0*ones(size(xi)), n*ones(size(xi)), c*ones(size(xi)));
  ddxSxx = cdiff(dx,Sxx);
  udum = u(:,3);
  %udum = udum-mean(udum);
  taub = 1000*cf*udum.*abs(udum);
  dissb = mean(-taub.*udum);
  De = .041; 
  Cd = 1;
  B = 3.66;L = 18;
  if contains(basename,'Base')
    N = 0;
  elseif contains(basename,'High')
    N= 50*8/(L*B); %number of plants ( and roots) per unit area
  end
  [j1 j2] = min(abs(xi-43));
  %% Cd Calculations (and some things for plots)
  % get exact value of Cd
  F2overCd = N*1000/2*De*(hv+eta_p(:,3)).*udum.*abs(udum);% uses p(3) for eta
  dissvegoverCd2 = mean(-F2overCd.*udum);
  Cdexact2 =((9810*n*c/8)*(Hrmsi(end)^2-Hrmsi(1)^2) - dissb*L)/(dissvegoverCd2*L); 
  if N==0;Cdexact = 1;Cdexact2 = 1;end
  %Cd = Cdexact;
  F2 = N*1000*Cdexact2/2*De*(hv+eta_p(:,3)).*udum.*abs(udum); %fluid force on vegetation
  d = -(abs(F2).*abs(udum)); %total time averaged force and dissipation
  dissveg2 = -mean(abs(F2).*abs(udum));
  
  modelEf(1) = (9810*n*c/8)*Hrmsi(1)^2;
  modelHrms(1) = sqrt(8*modelEf(1)/(9810*n*c)); 
  datEf(1) = (9810*n*c/8)*Hrmsi(1)^2;
  datHrms(1) = sqrt(8*modelEf(1)/(9810*n*c));
  in.cd = 1*Cdexact2;    in.d = hv;
  in.T = stats.Tp;  in.N = N;  in.D = De;  in.alpha = 0;  in.vegheight = Inf;
  %in.U = undertow_linear(in.d,in.Hrms,in.T);
  in.U = mean(udum);
  in.Hrms = 1.*modelHrms(1);
  in.inonlin = 1;
  [out] = veg_stress_dissipation(in);%function
  sav(1) = out;
  modeldissveg(1) = out.meanDx;
  modeleta(1) = (mean(p(:,1))-p_init(1))/9810;
  modelh(1) = hv+modeleta(1);
  dateta(1) = (mean(p(:,1))-p_init(1))/9810;
  eta0a(1) = (mean(p(:,1))-p_init(1))/9810;
  eta0b(1) = (mean(p(:,1))-p_init(1))/9810;
  [Sxy Sxx(1) Syy] = radstress(modelHrms(1),0, n, c);
  for k = 2:length(xi)
    datEf(k) = datEf(k-1)+dx*(dissb+dissveg2);
    modelEf(k) = modelEf(k-1)+dx*(dissb+modeldissveg(k-1));
    modelHrms(k) = sqrt(8*modelEf(k)/(9810*n*c));
    datHrms(k) = sqrt(8*datEf(k)/(9810*n*c));
    in.Hrms = modelHrms(k);
    [out] = veg_stress_dissipation(in);
    sav(k) = out;
    modeldissveg(k) = out.meanDx;
    modelF(k) = out.meanFx;
    [Sxy Sxx(k) Syy] = radstress(modelHrms(k),0, n, c);
    modeleta(k)= modeleta(k-1)-1*(Sxx(k)-Sxx(k-1))/(9810*hv)  - (dx/(9810*hv))*(mean(taub)+modelF(k)); 
    dateta(k)= dateta(k-1)  - (dx/(9810*hv))*(ddxSxx(k)+mean(taub)+mean(F2)); 
    eta0a(k)= eta0a(k-1)-1*(Sxx(k)-Sxx(k-1))/(9810*hv)  - (dx/(9810*hv))*(mean(taub)+0*modelF(k)); 
    eta0b(k)= eta0b(k-1)  - (dx/(9810*hv))*(ddxSxx(k)+mean(taub)+0*mean(F2)); 
  end
  
  Cdtable(j) = Cdexact2 ;  
%% Re Calculations
    Re(j) = fluiddensity*mean(abs(udum))*Daverage/fluidviscosity ; 
%% KC Calculations
   waveperiod(j) = range(sav(j2).t+offset(j)) ; %range(xi) 7.45 from Kelty paper
   KC(j) = mean(abs(udum)) * waveperiod(j) / Daverage ; 
  %% Plotting (original)
  if plottingog ==1    
  figure(cnt);clf;clear hh hl;
     subplot(2,1,1)
    %plot(xp,sqrt(8)*std(p)/9810,'bs-');hold on
     hh1(1) = plot(xp(1:6),sqrt(8)*std(eta_p(:,1:6)),'bs-','markerfacecolor','b');hold on
    %plot(xuswg,sqrt(8)*std(etaus),'rv','markerfacecolor','r');hold on
     hh1(2) = plot(xwg,sqrt(8)*std(eta),'rs-','markerfacecolor','r');hold on
     hh1(3) = plot(xi,modelHrms,'k-','linewidth',2);hold on
     hh1(3) = plot(xi,datHrms,'k--','linewidth',2);hold on      %dname2=dname;dname2(dname=='_')='-';
     title([strrep(dname(21:end),'_','-'),' ',sprintf('%2.2f',Hrmsi(1)),' ',sprintf('%2.2f',stats.Tp)],'interpreter','latex','fontsize',12)
    %xlabel('$x[m]$','interpreter','latex','fontsize',16)
     ylabel('$H_{rms}[m]$','interpreter','latex','fontsize',16)
     a = axis;
     axis manual
     hf=fill([36 54 54 36],[a(3) a(3) a(4) a(4)],[.2 .5 .2]);set(hf,'facealpha',.2) ;
     text(37,a(3)+.1*(a(4)-a(3)),'Veg. Section','interpreter','latex','fontsize',14) ;
     text(15,a(3)+.1*(a(4)-a(3)),['$C_d = ',sprintf('%2.2f',Cdexact2),'$'],'interpreter','latex','fontsize',14) ;

      subplot(2,1,2);
      plot(xp(1:6),1000*(mean(p(:,1:6))-p_init(1:6))/9810,'bs-','markerfacecolor','b');hold on
      plot(xwg,1000*(mean(eta)-eta_init),'rs-','markerfacecolor','r');hold on
      %plot(xuswg,1000*(mean(etaus)-etaus_init),'rv','markerfacecolor','r');hold on
      hh(1) = plot(xi,1000*modeleta,'k-','linewidth',2);
      hh(2) = plot(xi,1000*dateta,'k--','linewidth',2);
      hh(3) = plot(xi,1000*eta0a,'g-','linewidth',2);
      %  hh(2) = plot(xi,1000*eta0b,'g--','linewidth',2);
      %hh(2) = plot(xi,1000*etaiwavesbotshearveg,'k-.','linewidth',2);
      %plot(xi,1000*etaiwavesbotshearwrs,'m--','linewidth',2);
      %plot(xi,1000*etaiwavesbotshearveg2p,'m-','linewidth',2);
      %hh(2) = plot(xi,1000*etaiwavesbotshearveg2,'k-','linewidth',3);
      hl = legend(hh,'Modeled $ \eta, u $','Measued $\eta, u $','${\bf F}= 0$','AutoUpdate','off');
      %hl = legend(hh,'$S_{xx}, \tau_b$','$S_{xx}, \tau_b, {\bf F}$','AutoUpdate','off');
      set(hl,'interpreter','latex','fontsize',14,'location','northwest')
      dumyl = ylim;
      if dumyl(2)-dumyl(1)<100
        ylim([mean(ylim)-40 mean(ylim)+40])
      end
      xlabel('$x[m]$','interpreter','latex','fontsize',16)
      ylabel('$\overline{\eta}[mm]$','interpreter','latex','fontsize',14)
      set(gca,'TickLabelInterpreter','latex') ;

  dname2=dname;dname2(dname=='_'|dname=='/')='-';
%   print('-dpng',[dname2(21:end-1),'.png']) 
  figname = [basename, '_' ,  dnames(jj).name]    % should give you  Baseline_h158_hv073_NoWall_Trial01'
  print('-dpng',[figname,'.png']); end
%   a = axis;
%   hf=fill([36 54 54 36],[a(3) a(3) a(4) a(4)],[.2 .5 .2]);set(hf,'facealpha',.2)
%   axis(a)
  
if plottingog ==1
  figure(2);clf
  subplot(3,1,1)
  plot(t,udum,'r','linewidth',2);hold all
    plot(t,0*udum,'k','linewidth',1);hold all
  plot(sav(j2).t+offset(j),mean(sav(j2).u),'b','linewidth',2)
  ylabel('$u$','interpreter','latex','fontsize',16)
  %xlabel('$t[s]$','interpreter','latex','fontsize',16)
  xlim([0 20]);
  title([strrep(dname(21:end),'_','-'),' ',sprintf('%2.2f',Hrmsi(1)),' ',sprintf('%2.2f',stats.Tp)],'interpreter','latex','fontsize',12)
  %print('-dpng',[dname2(21:end-1),'u.png'])

  subplot(3,1,2)
  plot(t,-d,'r','linewidth',2);hold all
  plot(t,0*t+mean(-d),'r--','linewidth',1);hold all
  plot(sav(j2).t+offset(j),-sav(j2).Dx,'b','linewidth',2)
  plot(sav(j2).t+offset(j),ones(size(sav(j2).Dx))*mean(-sav(j2).Dx),'b--','linewidth',1)
  ylabel('$-D$','interpreter','latex','fontsize',16)
  %xlabel('$t[s]$','interpreter','latex','fontsize',16)
  xlim([0 20]);
  %title([strrep(dname(21:end),'_','-'),' ',sprintf('%2.2f',Hrmsi(1)),' ',sprintf('%2.2f',stats.Tp)],'interpreter','latex','fontsize',12)
  %print('-dpng',[dname2(21:end-1),'u.png'])
  
  subplot(3,1,3)
  plot(t,0*F2,'k');hold all
  hf(1) = plot(t,F2,'r','linewidth',2);hold all
  plot(t,mean(F2)+0*t,'r--');hold all
  hf(2) = plot(sav(j2).t+offset(j),sav(j2).Fx,'b','linewidth',2)
  plot(sav(j2).t+offset(j),sav(j2).meanFx+0*sav(j2).t,'b--','linewidth',1)
  ylabel('$F_x $','interpreter','latex','fontsize',16)
  xlabel('$t[s]$','interpreter','latex','fontsize',16)
  xlim([0 20]);
  %title([strrep(dname(21:end),'_','-'),' ',sprintf('%2.2f',Hrmsi(1)),' ',sprintf('%2.2f',stats.Tp)],'interpreter','latex','fontsize',12)
  hl = legend(hf,'Using Measured Data','Using Modeled Data','AutoUpdate','off');
  set(hl,'interpreter','latex','fontsize',14,'location','northeast') ; end
%  print('-dpng',[dname2(21:end-1),'_u_d_F.png']) EDITED CLARA
  
  %% Prepping Error Section
  modeltime = sav(j2).t+offset(j) ; % used to adjust times so they are the same
  errortimestart =  find(abs(t - modeltime(1)) == min(abs(t - modeltime(1)))) ; %381 ; 
  errortimeend = find(abs(t - modeltime(end)) == min(abs(t - modeltime(end))))-1 ; %1206
  errortime = t(errortimestart:errortimeend) ; 
  
  modeludum = mean(sav(j2).u) ;
  modelD = -sav(j2).Dx ; 
  modelF2 = sav(j2).Fx ; 
  uduminterpol = interp1(modeltime,modeludum,errortime,'linear');  % this interpolates modeled udum to find its predicted values in the times udum (errortime)
  Dinterpol = interp1(modeltime,modelD,errortime,'linear');  
  F2interpol = interp1(modeltime,modelF2,errortime,'linear');  
  %% Error Section, still in original loop
  errorHrms = modelHrms - Hrmsi; %wave height error
  erroreta = modeleta - eta0a ; %eta (mean water surface elevation) error
  errorudum = uduminterpol - udum(errortimestart:errortimeend)' ; %velocity error
  errorD = Dinterpol - -d(errortimestart:errortimeend)'; % total time-averaged force and dissipation error (opposite d is plotted)
  errorF2 = F2interpol - F2(errortimestart:errortimeend)'; % fluid force on vegetation error
  
%   errorData = [[errorHrms, NaN(1,abs(errortimestart-errortimeend)+1-numel(errorHrms))]; [erroreta, NaN(1,abs(errortimestart-errortimeend)+1-numel(erroreta))]; errorudum; errorD; errorF2] ; %need to make data the same length with NaN so they can be put in the same boxplot
%   errorlabel = {'Hrms', 'eta','U','D','F2'} ;
  errorData = [ [errorHrms, NaN(1,abs(errortimestart-errortimeend)+1-numel(errorHrms))];
                [erroreta, NaN(1,abs(errortimestart-errortimeend)+1-numel(erroreta))]; 
                [errorudum, NaN(1,abs(errortimestart-errortimeend)+1-numel(errorudum))]; 
                [errorD, NaN(1,abs(errortimestart-errortimeend)+1-numel(errorD))]; 
                [errorF2, NaN(1,abs(errortimestart-errortimeend)+1-numel(errorF2))] ; 
               ] ; %need to make data the same length with NaN so they can be put in the same boxplot
  errorlabel = {'Hrms', 'eta','U','D','F2'} ;
  unitslabel = {'Wave Height Attenuation [m]', 'Free Surface Elevation [mm?]', 'Velocity [m/s]', 'Dissipation [W?]', 'Force [N]'} ;
 % legendEntries = {'Data 1', 'Data 2', 'Data 3', 'Data 4', 'Data 5', 'Data 6'} ;
  legendEntries = arrayfun(@(x) ['Data ', num2str(x)], 1:length(dnames), 'UniformOutput', false); %x variable is local to this line, and can be used elsewhere
 % Calculate RMS and Bias errors for each dataset
  rmserrors = zeros(1, 5);
  biaserrors = zeros(1, 5);
  for ii = 1:5
      errorData2(ii, isnan(errorData(ii,:))) = 0 ;
      rmserrors(ii) = sqrt(nanmean( (errorData(ii,:) ).^2)); %
      biaserrors(ii) = nanmean(errorData(ii,:)); 
  end

% Prepare the results as strings
rms_text = sprintf('RMS Errors: %.2f, %.2f, %.2f, %.2f, %.2f', rmserrors);
bias_text = sprintf('Bias Errors: %.2f, %.2f, %.2f, %.2f, %.2f', biaserrors);

%% Phase Averaged Error Section
%     PAaxis1 = linspace(4,5, 360) ; %for error time, have to have range of times that overlap in all data sets
%     PAaxis2 = linspace(40,43.6,360) ; %for xi
    PAxitime = linspace(xi(1),xi(end),1000) ; 
    PAerrortime = linspace(errortime(1),errortime(end),1000) ;
    PAaxis = linspace(0, 360, 1000) ;
    errorHrmsPA = interp1(xi,errorHrms , PAxitime) ;
    erroretaPA = interp1(xi, erroreta, PAxitime) ; 
    errorudumPA = interp1(errortime, errorudum, PAerrortime) ; 
    errorDPA = interp1(errortime, errorD, PAerrortime) ; 
    errorF2PA = interp1(errortime, errorF2, PAerrortime) ; 
    PAerrorData(:,:,j) = [errorHrmsPA;erroretaPA;errorudumPA;errorDPA;errorF2PA] ; 
    HrmserrorData(:,j) = [errorHrms] ; 
    HrmsAxis(:,j) = [xi];

  %% Plotting Error Box Plots
  if plottingyn == 1 
      figure(cnt) ; %change to figure(cnt) if you want to see each boxplot
      hold off
      for num = 1:5
          subplot(1, 5, num);
          boxplot(errorData(num,:));  % Replace with your data (e.g., errorHrms)
          title(errorlabel(num));
          set(gca, 'Position', [(num-1)*0.19 + 0.05, 0.25, 0.15, 0.65]); %this is the text box positioning
          % Add the RMSE and Bias textboxes below each subplot
          annotation('textbox', [(num-1)*0.19 + 0.05, 0.1, 0.15, 0.1], 'String', 'RMSE: '+string(rmserrors(num)), 'EdgeColor', 'none', 'FontSize', 10, 'HorizontalAlignment', 'center');
          annotation('textbox', [(num-1)*0.19 + 0.05, 0.00, 0.15, 0.1], 'String', 'Bias: '+string(biaserrors(num)), 'EdgeColor', 'none', 'FontSize', 10, 'HorizontalAlignment', 'center');
      end
      sgtitle('Error Boxplots') ; end
  
%% Plotting Error Line Graphs NOT Phase Averaged
  if plottingyn == 1 %end is last line of plot
    figure(8);hold on;
    for num = 1:5
        subplot(5,1,num) ; 
        if num > 2 
            plot(errortime, errorData(num,:), 'LineWidth', 2) ; hold all
        else 
            plot(xi,errorData(num,1:length(xi)), 'LineWidth', 2) ; hold all %xi is what is used earlier against Hrms
        end
        title(errorlabel(num));
    end
    sgtitle('Error Plots') ;
    legend(legendEntries, 'Location', 'best'); end

 %% Storing Data so I can actually figure out what is going on
%     storagexistart(j) = xi(1);
%     storagexiend(j) = xi(end) ;
%     storageerrortime(j) = errortime(1);
%     storageerrortimeend(j) = errortime(end);
%     storage = [ "Hrmsi", "eta0a", "udum", "d", "F2"; 
%             string(Hrmsi(1)), string(eta0a(1)), string(udum(1)), string(d(1)), string(F2(1)) ];
end %all of this is done 6 times


%% Plotting Stored Data
if plottingpost ==1 
figure(cnt+10); hold on; % Plotting Phase Averaged Error Line Graphs
    cnt = cnt+1 ;  
    for num = 1:5
         subplot(3,2,num) ; 
            if num == 1 
                plot(HrmsAxis, HrmserrorData, 'LineWidth',2) ; hold all
                ylabel(unitslabel(num),'fontsize',10) ;
                xlabel('Space','fontsize',10) ;
            else
             for num2 = 1:length(dnames) %num 2 changed to length(dnames) from 6
             plot(PAaxis, PAerrorData(num,:, num2), 'LineWidth',2) ; hold all ; end
             ylabel(unitslabel(num),'fontsize',10) ;
             xlabel('Phase Angle','fontsize',10) ; 
         title(errorlabel(num)) ; 
            end
     end
     sgtitle('Phase Averaged Error Plots') ; 
     legend(legendEntries, 'Location', 'best') ; 
 figure(cnt+10) % scatterplot of Cd and Re
    cnt = cnt+1 ; 
    subplot(1,2,1) ; 
     scatter(Re, Cdtable, 'filled', 'g') ; hold on
     ReLof = polyfit(Re, Cdtable, 1) ;
     Cd_fit_Re = polyval(polyfit(Re, Cdtable, 1), Re); % line of best fit
     plot(Re, Cd_fit_Re, '--b')
     ylabel('Cd')
     xlabel('Re')
     title('Cd and Re Correlation')
     text(median(Re), median(Cdtable)*1.75, ['Slope = ', num2str(ReLof(1))], 'FontSize', 12, 'Color', 'red')
    subplot(1,2,2); % scatterplot of Cd and KC
     scatter(KC, Cdtable, 'filled', 'm') ; hold on
     KCLof = polyfit(KC, Cdtable, 1) ;
     Cd_fit_KC = polyval(polyfit(KC, Cdtable, 1), KC); % line of best fit
     plot(KC, Cd_fit_KC, '--r')
     ylabel('Cd')
     xlabel('KC')
     title('Cd and KC Correlation')
     text(median(KC), median(Cdtable)*1.75, ['Slope = ', num2str(KCLof(1))], 'FontSize', 12, 'Color', 'red')
end

