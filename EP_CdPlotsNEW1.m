clearvars -except aalldata corrdata; clc ; close all ; 
addpath('./ClaraFunctions') ; 
addpath('./data') ; 
addpath('./mfiles') ; 
     if ~exist('aalldata', "var") ;  load('aalldatanewCdcalc_20250721_2.mat') ; end
CdtableRand = readtable('Excel_SummaryTable_Random.xlsx') ; 
CdtableReg = readtable('Excel_SummaryTable_Regular.xlsx') ;

set(0,'defaultTextInterpreter','latex')
set(groot,'defaultAxesTickLabelInterpreter','latex') 
set(groot, 'defaultLegendInterpreter','latex')


currentfolder = pwd ; %MATLABONLINE
savfolderpath = join([currentfolder, '/ClaraFigures/RootGraphs/'], '') ; 

CdtableRand.sd = [NaN;NaN;0.29; 0.41;0.35;0.6;0.4;0.7;0.43;0.39;0.39;0.74;0.73;0.47;0.49;0.54;0.57;1.29;1.15;0.53;0.6] ;
totstruct.Regular = CdtableReg ; 
totstruct.Random = CdtableRand ; 
%% actual new code
clear fieldnames ; fieldnames = fieldnames(aalldata) ;
graphcolors = hsv(length(fieldnames)-4) ; 
legendnames = strrep(strrep(fieldnames, '_',' '), 'NoWall', '') ; 
legendnames{end+1} = 'Perfect Correlation' ; 
perfectcoor = linspace(0,4,10) ; 

liztable = [] ; 
loopmin = 5; %(should be the min value in the for loop for totalnum, all categories)
for totalnum = 5:length(fieldnames)

categoryname =  fieldnames{totalnum} ; 
set_category_variables
NumofTrials = length(t) ; 
for num = 1:NumofTrials
    [CdvalueK, sdvalueK, strJ, strK, Khrmsi, Jhrmsi, Ktp, Jtp] = match2KeltyTable(aalldata, totstruct, categoryname, num) ; 
    if isnan(CdvalueK)
        disp(categoryname + "Trial " + num + ". Cd value: NaN. Errorbars value: NaN") ;
    else
        disp(categoryname + "Trial " + num + ". Cd value: " + CdvalueK +". Errorbars value: "+ sdvalueK) ;
    end
    
    ltnames = {"Johnson", "Kelty", "Hrmsi_J","Hrmsi_K","Tp_J","Tp_K"} ; 
    var2 = [CdvalueK, sdvalueK, strJ, strK, Jhrmsi, Khrmsi, Jtp, Ktp] ; 
    if isempty(liztable)
        for num1 = 1:length(ltnames)
            jkjk = ltnames(num1) ; 
            val2 = var2(num1+2) ; 
            liztable.(jkjk{1})(1) = val2 ; 
        end
    end
    
    liztable.Johnson(end+1, 1) = strJ ; 
    liztable.Kelty(end+1,1) = strK ; 
    liztable.Hrmsi_J(end+1,1) = Jhrmsi ; 
    liztable.Hrmsi_K(end+1,1) = Khrmsi ; 
    liztable.Tp_J(end+1,1) = Jtp ; 
    liztable.Tp_K(end+1,1) = Ktp ; 

    uu = u{num} ; 
    tt = t{num} ; 
    uu0 = smoothNoise(uu, std(uu)/3, 10) ; %std(uu)/4
    uu0 = transpose(uu0) ; 
    tt0 = tt(1:end-10) ; %make sure you subtract the time in smoothNoise from tt in this line
    uu_new = spline(tt0, uu0,  1:.01:tt(end-10)) ;
    uu_n_xp = 1:.01:tt(end-10) ; %for plot
    eta_pp = eta_p{num} ; 
    eta_pp = eta_pp(:,3) ; 
%%%%%%% HERE
    F2overCd = Atm*N*1000/2*(hv+eta_p(:,3)).*udum.*abs(udum);% uses p(3) for eta
      dissvegoverCd2 = mean(-F2overCd.*udum);
      Cdexact2smoothed =((9810*nn*c/8)*(Hrmsi(end)^2-Hrmsi(1)^2) - dissb*Length)/(dissvegoverCd2*Length); 

    

%% plotting in for loop
   figure(1)
   graph(totalnum-(loopmin-1)) = scatter(CdvalueK, Cdexact2{num}, 36, graphcolors(totalnum-(loopmin-1),:), "filled") ; hold on
   xlabel('Kelty $C_D$')
   ylabel('Johnson $C_D$')
   title("$C_D$ Comparison: Johnson and Kelty values") 

   figure(2)
   graph2(totalnum-(loopmin-1)) = scatter(Jhrmsi, Cdexact2{num}, 36, graphcolors(totalnum-(loopmin-1),:), "filled") ; hold on
   graph2(totalnum-(loopmin-1)) = scatter(Khrmsi, CdvalueK, 36, graphcolors(totalnum-(loopmin-1),:), "filled", "^", 'MarkerEdgeColor', 'black') ; hold on
   xlabel("$H_s$")
   ylabel("$C_D$")
   title("$C_D$ vs. $H_s$")

   figure(3)
   if ~isnan(CdvalueK)
   graph3(totalnum-(loopmin-1)) = scatter(Re{num}, Cdexact2{num}, 36, graphcolors(totalnum-(loopmin-1),:), "filled") ; hold on
   end
   graph3(totalnum-(loopmin-1)) = scatter(Re{num} , CdvalueK, 36, graphcolors(totalnum-(loopmin-1),:), "filled", "^", 'MarkerEdgeColor', 'black') ; hold on
   xlabel("$Re$")
   ylabel("$C_D$")
   title("$C_D$ vs. $Re$ Matched Trials")

end
end
figure(1); legend(graph, legendnames(loopmin:end))
graph(end+1) = plot(perfectcoor, perfectcoor, 'LineWidth', 2, 'Color', '#c0c0c0') ; hold on
liztable2 = struct2table(liztable) ; 

figure(2); legend(graph2, legendnames(loopmin:end))
graph4(1) = scatter(-1,-1, 36, 'black', 'filled');
graph4(2) = scatter(-1,-1, 36, 'black', 'filled', '^');
xlim([0 1])
ylim([0 60])
legend(graph4, {"Johnson","Kelty"})

figure(3); legend(graph3, legendnames(loopmin:end))
graph4(1) = scatter(-1,-1, 36, 'black', 'filled');
graph4(2) = scatter(-1,-1, 36, 'black', 'filled', '^') ;
xlim([0 2e4])
ylim([0 60])
legend(graph4, {"Johnson","Kelty"})


% save('KeltyandJohnsonMatched2.mat', 'liztable2')
% writetable(liztable2)