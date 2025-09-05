% data processing section of explore_prediction
% files created through ClaraSummarize
% created: 2025/02/24 ; last edit: 2025/07/21
% clear ; clc ; close all
plottingog = 0 ; plottingyn = 0 ; plottingpost = 0 ; %1 means plots, 0 means no plots
%og section  % plotting error    %plotting stored data
cnt = 0; %figure numbers
addpath('./ClaraFunctions') ; 
addpath('./data') ; 
addpath('./mfiles') ; 
load('keltyprofile.mat') ; 
kelty_At = transpose(kelty_At) ; 
kelty_z = kelty_z ; 
%% constants
cf = .05;
icorrect =1;
Dtrunk = 0.1143 ; Droot = 0.0286 ; Daverage = 0.038 ;%both in meters, diameter of trunk and roots of model from Kelty full paper, average is ave of effective diameters p52
fluiddensity = 1000 ; fluidviscosity = .001 ; %kg/m3, kg/ms
offset = [3.0 3.9 2.8 3.5 3.8 3.8 ]; %this line is why j can only go to 6
BWidth = 3.66;Length = 18; %width and length of flume
zu = [1.4040, 1.5500, 1.7200, 1.8580] ; % ADV 2,3,4,5 (order is correct and checked) ; from the bottom of the wave flume, not the false floor
flumesectorA = 3.66 * 3.7 ; %in meters from Kelty thesis, this is for the sectors (1-13)
%% accessing files
currentfolder = pwd ; 
basename = join([currentfolder, '/data/SummaryFiles/'], '') ; 
dnames = dir([basename, '*Summary.mat']) ; 

summaryfiles = organize_files(dnames, '^(HighDensity|LowDensity|Baseline)_[\w]+_[\w]+(_\w+)*(?=_Trial)') ; 
sumcategories = fieldnames(summaryfiles) ;
exceltable = readtable('Experimental_Info.xlsx') ;
exceltable.TrialName = strings(height(exceltable),1) ; 
for i = 1:height(exceltable)
    trialnum = exceltable.Trial(i) ; 
    if isnan(trialnum)
        exceltable.TrialName(i) = "" ;
    else 
        if exceltable.Trial(i) < 10
            exceltable.TrialName(i) = join([exceltable.ExperimentFileName{i},'_Trial0',string(trialnum),'_Summary.mat'], '') ; 
        else 
            exceltable.TrialName(i) = join([exceltable.ExperimentFileName{i},'_Trial',string(trialnum),'_Summary.mat'], '') ; 
        end
    end
end
        
aalldata = struct() ; 
skippedtrials = {} ; 
modeledTplist = {} ; 
%% data processing
cz_xptab = zeros(1,length(dnames)) ;
for j = 1:length(dnames)
      clear p u eta eta_p ubp wbp categoryname N 
%% % creating the aalldata structure aalldata.(categoryname).
  cnt = cnt+1;
  offset = [offset, zeros(1, length(dnames) - length(offset))] ; % making offset the same length as dnames


    dname = [basename,dnames(j).name,'/'] ;
    file = dnames(j).name ; %summary mat file name
    load([basename, dnames(j).name])
    categoryname = string(regexp(dnames(j).name,  '^(HighDensity|LowDensity|Baseline)_[\w]+_[\w]+(_\w+)*(?=_Trial)', 'match')) ; 
  
  hv = str2num(dname(strfind(dname,'hv')+2:strfind(dname,'hv')+4))/100+.03; %finding hv in flies, calculating it, changing to number from string
  p = [dat.press.press];%gets press field of data, makes into set list
  u = [dat.u.u];u = u(:,2:5);%same as line 22 but for velocity
  xu = [dat.u.x];xu = xu(:,2:5);%velocity in x dir
  zu = [dat.u.z];zu = zu(:,2:5);%velocity in z dir
  w = [dat.w.w];w = w(:,2:5); %only velocity within mangroves is being saved
  zw = [dat.w.z];zw = zw(:,2:5); % z coor of ADV
  xp = [dat.press.x];%access x field of press field of the data
  [~,index] = ismember(35.8930,xp) ; 
  cz_xptab(j) = index ; 

  eta = [dat.wg.eta]; % water level elevation (where the water line is)
  xwg = [dat.wg.x];
  etaus =[dat.uswg.eta];
  xuswg = [dat.uswg.x];
    
    if contains(categoryname, 'Baseline') || any(u(:,3)~=0) 
      activity = mean(std(eta),2); %average of standard deviation of each row
      startinds = find(mean(abs(eta(:,find(xwg>50))),2)>activity);startind = startinds(10);
      endinds = find(mean(abs(eta(:,find(xwg<50))),2)>activity);endind = endinds(end-10);
      p_init = mean(p(1:round(startind/2),:));
      eta_init = mean(eta(1:round(startind/2),:)); %average of eta rows 1 to startind/2 rounded to nearest integer and all columns
      etaus_init = mean(etaus(1:round(startind/2),:));%see line 39
      eta = eta(startind:endind,:);
      t = [0:size(eta,1)-1]./100;%element wise division of array of 1st col of eta
      stats = find_stats(t,eta(:,2),4);%function find_stats
      [k,nn,c] = dispersion (2*pi/stats.Tp,hv);
      etaus = etaus(startind:endind,:);
      p = p(startind:endind,:) ;
      u = u(startind:endind,:);
      w = w(startind:endind,:);
      for jj =1:size(u,2)
        ubp(:,jj) = bandpass(u(:,jj),1/100,.5*1/stats.Tp,2*1/stats.Tp);%filters signals so only those within the specified range can pass
        wbp(:,jj) = bandpass(w(:,jj),1/100,.5*1/stats.Tp,2*1/stats.Tp);
      end
    
      for i = 1:size(p,2)
        if dat.press(i).z-dat.press(i).swd<0
          eta_p(:,i) = p2eta([p(:,i)-mean(p(:,i))],1/100,hv,dat.press(i).z-dat.press(i).swd);
        else
          eta_p(:,i) = dat.press(i).z+p(:,i)/9810;
        end
      end
      dx = .1;
      xi = dat.press(1).x:dx:dat.press(6).x; 
      Hrmsi = interp1(xwg,sqrt(8)*std(eta),xi);
      %Hrmsi = interp1(xp,sqrt(8)*std(eta_p),xi);
      [Sxy, Sxx, Syy] = radstress (Hrmsi, 0*ones(size(xi)), nn*ones(size(xi)), c*ones(size(xi)));
      ddxSxx = cdiff(dx,Sxx);
      udum = u(:,3); %for j = 44, 75, ... the 3rd col is just 0
      %udum = udum-mean(udum);
      taub = 1000*cf*udum.*abs(udum);
      dissb = mean(-taub.*udum);

      Cd = 1;
      
      if contains(dnames(j).name,'Base')
          N = 0;
      elseif contains(dnames(j).name,'High') 
          N= 10 / flumesectorA; %number of trees/m2
      elseif contains(dnames(j).name, 'Low')
          N = 5/flumesectorA; %number of trees/m2
      end
      if contains(dnames(j).name, '185') || contains(dnames(j).name, '182')
          zz = zu(4) ; 
          Atm = mean(kelty_At(1:zz)) ;
      elseif contains(dnames(j).name, '148')
          zz = zu(3) ; 
          Atm = mean(kelty_At(1:zz)) ;
      elseif contains(dnames(j).name, '103') || contains(dnames(j).name, '100')
          zz = zu(2) ; 
          Atm = mean(kelty_At(1:zz)) ;
      elseif contains(dnames(j).name, '073') || contains(dnames(j).name, '070')
          zz = zu(1) ; 
          Atm = mean(kelty_At(1:zz)) ;
      else 
          disp("Area term undefined based on given title/hv parameters") ;
      end

      [j1, j2] = min(abs(xi-43));
      
    % Cd Calculations (and some things for plots)
      % get exact value of Cd
      F2overCd = Atm*N*1000/2*(hv+eta_p(:,3)).*udum.*abs(udum);% uses p(3) for eta
      dissvegoverCd2 = mean(-F2overCd.*udum);
      Cdexact2 =((9810*nn*c/8)*(Hrmsi(end)^2-Hrmsi(1)^2) - dissb*Length)/(dissvegoverCd2*Length); 
      if N==0;Cdexact = 1;Cdexact2 = 1;end
      %Cd = Cdexact;
      F2 = Atm*N*1000*Cdexact2/2*(hv+eta_p(:,3)).*udum.*abs(udum); %fluid force on vegetation
      d = -(abs(F2).*abs(udum)); %total time averaged force and dissipation
      dissveg2 = -mean(abs(F2).*abs(udum));
      Lwavelength = 2*pi / k ; 
% smoothed start
    uu = u ; % num is supposed to be trials
    tt = t ; 
    smr = 10 ; 
    uu0 = smoothNoise(uu, std(uu)/3, smr) ; %std(uu)/4
    uu0 = transpose(uu0) ; 
    tt0 = tt(1:end-smr) ; %make sure you subtract the time in smoothNoise from tt in this line
    uu_new = spline(tt0, uu0,  1:.01:tt(end-smr)) ; %same comment as above
    uu_new = transpose(uu_new) ; 
    uu_n_xp = 1:.01:tt(end-smr) ; %for plot
    dissbsmoothed = mean(-taub(1:end-100-smr).*uu_new);
    F2overCdsmoothed = Atm*N*1000/2*(hv+eta_p(1:end-100-smr,3)).*uu_new.*abs(uu_new);% uses p(3) for eta
    dissvegoverCd2smoothed = mean(-F2overCdsmoothed.*uu_new);
    Cdexact2smoothed =((9810*nn*c/8)*(Hrmsi(end)^2-Hrmsi(1)^2) - dissbsmoothed*Length)/(dissvegoverCd2smoothed*Length); 
    F2smoothed = Atm*N*1000*Cdexact2smoothed/2*(hv+eta_p(1:end-100-smr,3)).*uu_new.*abs(uu_new); %fluid force on vegetation
    dsmoothed = -(abs(F2smoothed).*abs(uu_new)); %total time averaged force and dissipation
    dissveg2smoothed = -mean(abs(F2smoothed).*abs(uu_new));

    Resmoothed = fluiddensity*mean(abs(uu_new))*Daverage/fluidviscosity ; 
% smoothed end
       
      
      modelEf(1) = (9810*nn*c/8)*Hrmsi(1)^2;
      modelHrms(1) = sqrt(8*modelEf(1)/(9810*nn*c)); 
      datEf(1) = (9810*nn*c/8)*Hrmsi(1)^2;
      datHrms(1) = sqrt(8*modelEf(1)/(9810*nn*c));
      in.cd = 1*Cdexact2;    in.d = hv;
      in.T = stats.Tp;  in.N = N;  in.alpha = 0;  in.vegheight = Inf; in.D = Atm;
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
      [Sxy, Sxx(1), Syy] = radstress(modelHrms(1),0, nn, c);
      for ii = 2:length(xi)
        datEf(ii) = datEf(ii-1)+dx*(dissb+dissveg2);
        modelEf(ii) = modelEf(ii-1)+dx*(dissb+modeldissveg(ii-1)); %Clara @k=2, modeldissveg is NaN
        modelHrms(ii) = sqrt(8*modelEf(ii)/(9810*nn*c));
        datHrms(ii) = sqrt(8*datEf(ii)/(9810*nn*c));
        in.Hrms = modelHrms(ii);
        [out] = veg_stress_dissipation(in);
        sav(ii) = out;
        modeldissveg(ii) = out.meanDx;
        modelF(ii) = out.meanFx;
        [Sxy, Sxx(ii), Syy] = radstress(modelHrms(ii),0, nn, c);
        modeleta(ii)= modeleta(ii-1)-1*(Sxx(ii)-Sxx(ii-1))/(9810*hv)  - (dx/(9810*hv))*(mean(taub)+modelF(ii)); 
        dateta(ii)= dateta(ii-1)  - (dx/(9810*hv))*(ddxSxx(ii)+mean(taub)+mean(F2)); 
        eta0a(ii)= eta0a(ii-1)-1*(Sxx(ii)-Sxx(ii-1))/(9810*hv)  - (dx/(9810*hv))*(mean(taub)+0*modelF(ii)); 
        eta0b(ii)= eta0b(ii-1)  - (dx/(9810*hv))*(ddxSxx(ii)+mean(taub)+0*mean(F2)); 
      end
% Kelty Cd Calculation (written by Clara)
    rownum = find(strcmp(exceltable.TrialName, file),1) ; 
    wavetype = exceltable.WaveType{rownum} ;
    Tp = exceltable.Tp(rownum) ; 
    if isnan(Tp)
        Tp = exceltable.TargetT_Tp_FullStrokeDuration(rownum) ; 
        if isempty(modeledTplist)
            modeledTplist{1} = file ; 
        else
        modeledTplist{end+1} = file ; 
        end
    end

    alpha = exceltable.waveHeightDecayCoefficient(rownum) ; 
    Atm = Daverage ; %average projected area per height per tree (so i am using the weighted average diameter
    vegweth = hv ; %mean wetted height of vegetation (d)
    avewatd = hv ; %mean water depth at vegetation (h)
    Hi = 1.4*Hrmsi(1) ; %calculating incident wave height
    % if contains(file, 'HighDensity') || contains(file, 'LowDensity')
    if strcmp(wavetype, 'Regular')
        denominator = 9*pi*(sinh(k*avewatd)*(sinh(2*k*avewatd)+2*k*avewatd)) ; 
        numerator = 4*Atm*N*Hi*k*((sinh(k*vegweth)^3)+3*sinh(k*vegweth)) ;
        CdKelty = alpha * denominator / numerator ; 
        Hsjohnson = sqrt(2)*std(eta_p(:,1)); 
    elseif strcmp(wavetype, 'Random')
        denominator = 3*sqrt(pi)*(sinh(k*avewatd)*(sinh(2*k*avewatd)+2*k*avewatd)) ; 
        numerator = 4*Atm*N*Hrmsi(1)*k*((sinh(k*vegweth)^3)+3*sinh(k*vegweth)) ;
        CdKelty = alpha * denominator / numerator ;
        Hsjohnson = 4*std(eta_p(:,1)) ; 
    else 
        disp('Kelty Cd not calculated. Wave type not specified') ; end
    if contains(file, 'Baseline')
        CdKelty = 0 ; 
        alpha = 0 ; 
    end
        

% Re Calculations
    Re = fluiddensity*mean(abs(udum))*Daverage/fluidviscosity ; % Reynold's number
% KC Calculations
   waveperiod = range(sav(j2).t+offset(j)) ; %range(xi) 7.45 from Kelty paper
   KC = mean(abs(udum)) * waveperiod / Daverage ; 

% saving the data into a structure array
      aalldata = structure_variables(aalldata, categoryname, 'alpha', alpha) ;
      aalldata = structure_variables(aalldata, categoryname, 'Cdexact2', Cdexact2) ;
      aalldata = structure_variables(aalldata, categoryname, 'Cdsmoothed', Cdexact2smoothed) ;
      aalldata = structure_variables(aalldata, categoryname, 'd', d) ;
      aalldata = structure_variables(aalldata, categoryname, 'CdKelty', CdKelty) ;
      aalldata = structure_variables(aalldata, categoryname, 'datEf', datEf) ;
      aalldata = structure_variables(aalldata, categoryname, 'dateta', dateta) ;
      aalldata = structure_variables(aalldata, categoryname, 'datHrms', datHrms) ; 
      aalldata = structure_variables(aalldata, categoryname, 'eta', eta) ;
      aalldata = structure_variables(aalldata, categoryname, 'eta0a', eta0a) ; 
      aalldata = structure_variables(aalldata, categoryname, 'eta0b', eta0b) ;
      aalldata = structure_variables(aalldata, categoryname, 'eta_init', eta_init) ;
      aalldata = structure_variables(aalldata, categoryname, 'eta_p', eta_p) ;
      aalldata = structure_variables(aalldata, categoryname, 'eta_psmoothed', eta_p(1:end-smr,3)) ;
      aalldata = structure_variables(aalldata, categoryname, 'F2', F2) ;
      aalldata = structure_variables(aalldata, categoryname, 'F2smoothed', F2smoothed) ;
      % aalldata = structure_variables(aalldata, categoryname, 'F2overCd', F2overCd) ;
      aalldata = structure_variables(aalldata, categoryname, 'Hrmsi', Hrmsi) ; 
      aalldata = structure_variables(aalldata, categoryname, 'Hsjohnson', Hsjohnson) ;
      aalldata = structure_variables(aalldata, categoryname, 'hv', hv) ;
      aalldata = structure_variables(aalldata, categoryname, 'k', k) ;
      aalldata = structure_variables(aalldata, categoryname, 'KC', KC) ;
      aalldata = structure_variables(aalldata, categoryname, 'modelHrms', modelHrms) ;
      aalldata = structure_variables(aalldata, categoryname, 'modeleta', modeleta) ;
      aalldata = structure_variables(aalldata, categoryname, 'modeledTplist', modeledTplist) ;
      aalldata = structure_variables(aalldata, categoryname, 'p', p) ;
      aalldata = structure_variables(aalldata, categoryname, 'p_init', p_init) ; 
      aalldata = structure_variables(aalldata, categoryname, 'Re', Re) ;
      aalldata = structure_variables(aalldata, categoryname, 'Resmoothed', Resmoothed) ;
      aalldata = structure_variables(aalldata, categoryname, 'sav', sav) ;
      aalldata = structure_variables(aalldata, categoryname, 'stats', stats) ;
      aalldata = structure_variables(aalldata, categoryname, 't', t) ;
      aalldata = structure_variables(aalldata, categoryname, 'ttsmoothed', tt0) ;
      aalldata = structure_variables(aalldata, categoryname, 'Tp', Tp) ; %period from excel table
      aalldata = structure_variables(aalldata, categoryname, 'u', u) ; % ADV horizontal velocities (2:5), parallel to flow
      aalldata = structure_variables(aalldata, categoryname, 'uu0', uu0) ;
      aalldata = structure_variables(aalldata, categoryname, 'uu_new', uu_new) ;
      aalldata = structure_variables(aalldata, categoryname, 'uu_n_xp', uu_n_xp) ;
      % aalldata = structure_variables(aalldata, categoryname, 'udum', udum) ; % ADV 3/4 only
      % aalldata = structure_variables(aalldata, categoryname, 'w', w) ; %ADV vertical velocities
      aalldata = structure_variables(aalldata, categoryname, 'wavetype', wavetype) ;
      aalldata = structure_variables(aalldata, categoryname, 'waveperiod', waveperiod) ;
      aalldata = structure_variables(aalldata, categoryname, 'xi', xi) ;
      aalldata = structure_variables(aalldata, categoryname, 'xp', xp) ;
      aalldata = structure_variables(aalldata, categoryname, 'xwg', xwg) ; 
      aalldata = structure_variables(aalldata, categoryname, 'zu', zu) ; %ADV z coordinates

    else 
        disp("Trial " +dnames(j).name+ " has been skipped")
        if ~isempty(skippedtrials)
            skippedtrials{end+1} = dnames(j).name ; 
        else 
            skippedtrials{1} = dnames(j).name ; 
        end
    end 
end %CLARA

%%
  % save(['/home/elizabeth/Desktop/cshorex-main/osu_mangrove/data/', 'aalldata_20250418.mat'], 'aalldata') ; 
    % save(['/MATLAB Drive/ClaraZwolanek/data/', 'aalldatasmoothed_20250822'], 'aalldata') ; %MATLAB ONLINE