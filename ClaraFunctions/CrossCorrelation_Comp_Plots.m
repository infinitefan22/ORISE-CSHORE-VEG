% errbars = 0 ;
% indvlayoutplot = 0 ; %plotting indv layout graphs. 1 for yes, anything else for no (0)
%% Plotting prep
    startingval = Tp{num} ; % should change to what you want the color bar to match to
    kkhh = round(startingval/kfactor) ; 
   if contains(categoryname, "High")
       shapetype = 'o' ;
       shorttit = 'High Density ' ;
   elseif contains(categoryname, "Low")
       shapetype = 'square' ;
       shorttit = 'Low Density ' ;
   else
       shapetype = '^' ;
       shorttit = '' ;
   end
   if Qcosh ==1 ; etatext = '$eta*cosh$' ; else ; etatext = '$eta$' ; end
  
%% All Combined
   figure(200) ;
   for num = 1:NumofTrials
       if clrbr == 1 
           scatter(xc.etavu{num}, zu(1:totADV), 36, khcolors(kkhh,:), shapetype, 'LineWidth', 1.5) ; hold on 
       else 
           scatter(xc.etavu{num}, zu(1:totADV), 36, khcolors(kkhh,:), shapetype, 'LineWidth', 1.5) ; hold on ; 
       end
   end
   % gg1 = scatter([xc.mean{1:totADV}], zu(1:totADV), 100, graphcolors(gcnum,:),shapetype, 'LineWidth', 2) ; hold on %Mean
   % if ~exist('graphes1','var') ; graphes1(1) = gg1 ; else ; graphes1(end+1) = gg1 ; end
   if ~exist('legendentry','var') ; legendentry{1} = strrep(strrep(categoryname, '_NoWall', ''), '_', '-') ; else;legendentry{end+1} = strrep(strrep(categoryname, '_NoWall', ''), '_', '-') ; end
   if errbars ==1 ; errorbar([xc.mean{1:totADV}], zu(1:totADV), xc.sdmean(1:totADV),'horizontal','Color', graphcolors(gcnum,:), 'LineStyle', 'none') ; end
   title("Correlation "+etatext+" vs. $u$")
   xlim([.5 1.01])
   xlabel('Cross Correlation')
   ylim([1.2 2])
   yticks(zu) ; yticklabels(string(zu))
   ylabel('$z$ [m]')
   if clrbr ==1 
       colormap(gca,jet)
       colorbar = colorbar ; 
       colorbar.Label.String = cblab ; colorbar.Label.FontSize = 12 ; colorbar.Label.Interpreter = 'latex' ;
       clim([kmapmin kmapmax]) ; 
   end
    
    %see EP_Cosh for reference on legend for multiple trials but one shape being graphed
    hold on ; plot([1 1 1 1], [1 1.5 2 2.5], 'k')
    gg2(1) = scatter(0, 0, 1, 'black', 'o') ; gg2(2) = scatter(0,0,1,'black', 'square') ; 
    legend(gg2,{'High Density', 'Low Density'} )