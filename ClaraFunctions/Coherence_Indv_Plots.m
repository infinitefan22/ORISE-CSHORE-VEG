errbars=1;
%% Plotting prep
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

 %% Plotting Indv Layouts (like EP_VelocityProfiles)
    % if indvlayoutplot == 1  
        figure(cnt) ; cnt = cnt+1 ;
        trialcolors = hsv(NumofTrials) ; 
        clear trialnumbers
        for num = 1:NumofTrials
            trialnumbers(num) = join(['Trial ',string(num)], '') ; trialnumbers(end+1) = '$mean$ all trials' ; 
            if errbars == 0; scatter(ms.etavu{num}, zu(1:totADV), 36, trialcolors(num,:), 'o', 'filled', 'MarkerEdgeColor', 'black') ; hold on ; end
        end
        % scatter(vrmscategories.(categoryname), zu, 100, 'black', 'filled') %RMS
        scatter([ms.mean{1:totADV}], zu(1:totADV), 100, graphcolors(gcnum,:), 'filled', 'MarkerEdgeColor', "black") ; hold on %Mean
        if errbars ==1 ; errorbar([ms.mean{1:totADV}], zu(1:totADV), ms.sdmean(1:totADV),'horizontal','Color', graphcolors(gcnum,:), 'LineStyle', 'none') ; end
        if errbars ==1; legend({'Mean Coherence', 'Standard Deviation'}); else; legend(trialnumbers);end
        % xlim([0 1.01])
        xlabel('Coherence')
        ylim([1.2 2.2])
        yticks(zu) ; yticklabels(string(zu)) 
        ylabel('$z$ [m]')
        title(strrep(strrep(categoryname, '_NoWall', ''), '_', '-') + " Correlation "+etatext+" vs. $u$")
        set(gcf, 'Position', [10, 10, 600, 400]);
        savfigname = categoryname + "_ErrorBars_EtavU.png" ;
        hold on ; explotcohere1 = plot([1 1 1 1], [1 1.5 2 2.5], 'black')
        legend(explotcohere1, 'Perfect Coherence')
        if savfigs==1; pause(3); exportgraphics(gcf, fullfile(savfolderpath, savfigname), 'Resolution', 300) ; close all ;end
    % end