function scatterNicePlot(namePrefix,figSubtitle,dataName,dataType,label,ymin,ymax,var1,var2,colorToUse)

% namePrefix: namePrefix
% figSubtitle: baseline or odorID
% dataName: expData(1).name
% dataType: dataType
% label: label
% ymin: -ylimit.(dataType)
% ymax: ylimit.(dataType)
% var1: baselineAvgsToCompare(:,1)
% var2: baselineAvgsToCompare(:,2)

figure('Name',strcat(namePrefix, "_", dataType, "_", figSubtitle, "_niceplot"))

hold on;
scatter(var1,var2,50,'o','MarkerEdgeColor','none','MarkerFaceColor', colorToUse, 'MarkerFaceAlpha', 0.5)
% scatter(var1,var2,50,'o','MarkerEdgeColor','k','MarkerEdgeAlpha',0.5,'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.5)
plot([ymin,ymax],[ymin,ymax], '--k')
axis([ymin,ymax,ymin,ymax])
title(dataName,figSubtitle,'Interpreter','none')
axis square
set(gcf, 'Position', [50 50 300 400])    % x y width height       
yticks([ymin,0 ymax]);
yticks([ymin,0 ymax]);
set(gca,'TickDir','in');
set(gca,'FontName','Arial');
set(gca,'TickLength', [0.025, 0.025]);
set(findall(gcf,'-property','FontSize'),'FontSize',20)
ylabel(strcat(label.(dataType), " Post-injection"),'FontSize',25)
xlabel(strcat(label.(dataType), " Pre-injection"),'FontSize',25)
hold off;

end
