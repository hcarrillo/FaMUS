%%% Drawings for IROS 2012 Path Planning %%
clc
clear all
close all
%% 3D plot of uncertainty over the map
colors = [0.1 0.6 0.1];

Bicocca_opt;
h = figure(); plot3(cost{1}(:,2),cost{1}(:,3),cost{1}(:,4),'color',colors,'lineWidth',1);
grid on;
axis vis3d;
view([-38 22]);
% xlabel('x(m)','fontsize',14,'fontweight','b','color','k');
% ylabel('y(m)','fontsize',14,'fontweight','b','color','k');
zlabel('Uncertainty','fontsize',14,'fontweight','b','color','k');
title('\it{Bicocca}','FontSize',16)
% [az,el] = view
% thetax=-atan(sin(el*pi/180)*tan(az*pi/180))*(180/pi);
% thetay=atan(sin(el*pi/180)*cot(az*pi/180))*(180/pi);
% set(get(gca,'xlabel'),'rotation',thetax);
% set(get(gca,'ylabel'),'rotation',thetay);
saveTightFigure(h,'Bicocca_opt_3D_plot.eps')

manhattan_opt;
h = figure(); plot3(cost{1}(:,2),cost{1}(:,3),cost{1}(:,4),'color',colors,'lineWidth',1);
grid on;
axis vis3d;
view([-38 22]);
zlabel('Uncertainty','fontsize',14,'fontweight','b','color','k');
title('\it{Manhattan}','FontSize',16)
saveTightFigure(h,'Manhattan_opt_3D_plot.eps')

result_intel_opt;
h = figure(); plot3(cost{1}(:,2),cost{1}(:,3),cost{1}(:,4),'color',colors,'lineWidth',1);
grid on;
axis vis3d;
view([-38 22]);
zlabel('Uncertainty','fontsize',14,'fontweight','b','color','k');
title('\it{Intel}','FontSize',16)
saveTightFigure(h,'Intel_opt_3D_plot.eps')

result_new_college;
h = figure(); plot3(cost{1}(:,2),cost{1}(:,3),cost{1}(:,4),'color',colors,'lineWidth',1);
grid on;
axis vis3d;
view([-38 22]);
zlabel('Uncertainty','fontsize',14,'fontweight','b','color','k');
title('\it{Intel}','FontSize',16)
saveTightFigure(h,'New_college_opt_3D_plot.eps')
%% 2D plot of all the nodes representing the map

h = figure();
[V,E] = readg2ofile('Bicocca_opt.g2o');
plot(V(2,:),V(3,:),'color',[0.5 0.5 0.5]); axis equal; hold on;
plot(V(2,:),V(3,:),'o','color',[1.0 0 0]); axis([min(V(2,:))-5, max(V(2,:))+5 min(V(3,:))-5 max(V(3,:))+5]); axis equal; hold off;
hold on
[V,E]=readg2ofile('Bicocca_opt.g2o');
g2o_plotGraphClusters(V,E,[]); axis([min(V(2,:))-5, max(V(2,:))+5 min(V(3,:))-5 max(V(3,:))+5]); axis equal;  hold off;
set(gca, 'LineWidth', 0.5,'FontSize', 14,'FontWeight', 'bold')
xlabel('x(m)','fontsize',14,'fontweight','b','color','k');
ylabel('y(m)','fontsize',14,'fontweight','b','color','k');
title('\it{Bicocca}','FontSize',16)
hold off
saveTightFigure(h,'Bicocca_opt_2D_full_nodes_plot_1.eps')

h = figure();
[V,E]=readg2ofile('intel_opt.g2o');
plot(V(2,:),V(3,:),'color',[0.5 0.5 0.5]); axis equal; hold on;
plot(V(2,:),V(3,:),'o','color',[1.0 0 0],'MarkerSize', 4); axis([min(V(2,:))-1, max(V(2,:))+1 min(V(3,:))-1 max(V(3,:))+1]); axis equal; hold off;
hold on
g2o_plotGraphClusters(V,E,[]); axis([min(V(2,:))-1, max(V(2,:))+1 min(V(3,:))-1 max(V(3,:))+1]); axis equal; hold off;
set(gca, 'LineWidth', 0.5,'FontSize', 14,'FontWeight', 'bold')
xlabel('x(m)','fontsize',14,'fontweight','b','color','k');
ylabel('y(m)','fontsize',14,'fontweight','b','color','k');
title('\it{Intel}','FontSize',16)
hold off
saveTightFigure(h,'Intel_opt_2D_full_nodes_plot_1.eps')

h = figure();
[V,E]=readg2ofile('manhattanOlson3500_opt.g2o');
plot(V(2,:),V(3,:),'color',[0.5 0.5 0.5]); axis equal; hold on;
plot(V(2,:),V(3,:),'o','color',[1.0 0 0]); axis([min(V(2,:))-5, max(V(2,:))+5 min(V(3,:))-5 max(V(3,:))+5]); axis equal; hold off;
hold on
g2o_plotGraphClusters(V,E,[]); axis([min(V(2,:))-5, max(V(2,:))+5 min(V(3,:))-5 max(V(3,:))+5]); axis equal; hold off;
set(gca, 'LineWidth', 0.5,'FontSize', 14,'FontWeight', 'bold')
xlabel('x(m)','fontsize',14,'fontweight','b','color','k');
ylabel('y(m)','fontsize',14,'fontweight','b','color','k');
title('\it{Manhattan}','FontSize',16)
hold off
saveTightFigure(h,'Manhattan_opt_2D_full_nodes_plot_1.eps')


h = figure();
[V,E] = readg2ofile('NewCollege_opt.g2o');
plot(V(2,:),V(3,:),'color',[0.5 0.5 0.5]); axis equal; hold on;
plot(V(2,:),V(3,:),'o','color',[1.0 0 0]);  axis([min(V(2,:))-5, max(V(2,:))+5 min(V(3,:))-5 max(V(3,:))+5]); axis equal; hold off;
hold on
g2o_plotGraphClusters(V,E,[]); axis([min(V(2,:))-5, max(V(2,:))+5 min(V(3,:))-5 max(V(3,:))+5]); axis equal; hold off;
set(gca, 'LineWidth', 0.5,'FontSize', 14,'FontWeight', 'bold')
xlabel('x(m)','fontsize',14,'fontweight','b','color','k');
ylabel('y(m)','fontsize',14,'fontweight','b','color','k');
title('\it{New College}','FontSize',16)
hold off
saveTightFigure(h,'NewCollege_opt_2D_full_nodes_plot_1.eps')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=figure();
[V,E] = readg2ofile('Bicocca_opt.g2o');
plot(V(2,:),V(3,:),'color',[0.5 0.5 0.5]); hold on;
plot(V(2,:),V(3,:),'o','color',[1.0 0 0]); axis([min(V(2,:))-5, max(V(2,:))+5 min(V(3,:))-5 max(V(3,:))+5]); axis equal; hold off;
set(gca, 'LineWidth', 0.5,'FontSize', 14,'FontWeight', 'bold')
xlabel('x(m)','fontsize',14,'fontweight','b','color','k');
ylabel('y(m)','fontsize',14,'fontweight','b','color','k');
title('\it{Bicocca}','FontSize',16)
saveTightFigure(h,'Bicocca_opt_2D_full_nodes_plot.eps')

h=figure();
[V,E] = readg2ofile('intel_opt.g2o');
plot(V(2,:),V(3,:),'color',[0.5 0.5 0.5]); hold on;
plot(V(2,:),V(3,:),'o','color',[1.0 0 0]); axis([min(V(2,:))-1, max(V(2,:))+1 min(V(3,:))-1 max(V(3,:))+1]); axis equal; hold off;
set(gca, 'LineWidth', 0.5,'FontSize', 14,'FontWeight', 'bold')
xlabel('x(m)','fontsize',14,'fontweight','b','color','k');
ylabel('y(m)','fontsize',14,'fontweight','b','color','k');
title('\it{Intel}','FontSize',16)
saveTightFigure(h,'Intel_opt_2D_full_nodes_plot.eps')

h=figure();
[V,E] = readg2ofile('manhattanOlson3500_opt.g2o');
plot(V(2,:),V(3,:),'color',[0.5 0.5 0.5]); axis equal; hold on;
plot(V(2,:),V(3,:),'o','color',[1.0 0 0]); axis([min(V(2,:))-5, max(V(2,:))+5 min(V(3,:))-5 max(V(3,:))+5]); axis equal; hold off;
set(gca, 'LineWidth', 0.5,'FontSize', 14,'FontWeight', 'bold')
xlabel('x(m)','fontsize',14,'fontweight','b','color','k');
ylabel('y(m)','fontsize',14,'fontweight','b','color','k');
title('\it{Manhattan}','FontSize',16)
saveTightFigure(h,'Manhattan_opt_2D_full_nodes_plot.eps')

h=figure();
[V,E] = readg2ofile('NewCollege_opt.g2o');
plot(V(2,:),V(3,:),'color',[0.5 0.5 0.5]); axis equal; hold on;
plot(V(2,:),V(3,:),'o','color',[1.0 0 0]); axis([min(V(2,:))-5, max(V(2,:))+5 min(V(3,:))-5 max(V(3,:))+5]); axis equal; hold off;
set(gca, 'LineWidth', 0.5,'FontSize', 14,'FontWeight', 'bold')
xlabel('x(m)','fontsize',14,'fontweight','b','color','k');
ylabel('y(m)','fontsize',14,'fontweight','b','color','k');
title('\it{New College}','FontSize',16)
saveTightFigure(h,'NewCollege_opt_2D_full_nodes_plot.eps')
%% 2D plot of the decison points (nodes) used to represent the map
h=figure();
[V,E] = readg2ofile('reduced_graph_bicocca.g2o');
plot(V(2,:),V(3,:),'color',[0.5 0.5 0.5]); axis equal; hold on;
plot(V(2,:),V(3,:),'o','color',[1.0 0 0]); axis([min(V(2,:))-5, max(V(2,:))+5 min(V(3,:))-5 max(V(3,:))+5]); axis equal; hold off;
set(gca, 'LineWidth', 0.5,'FontSize', 14,'FontWeight', 'bold')
xlabel('x(m)','fontsize',14,'fontweight','b','color','k');
ylabel('y(m)','fontsize',14,'fontweight','b','color','k');
title('\it{Bicocca}','FontSize',16)
saveTightFigure(h,'Bicocca_opt_2D_reduced_nodes_plot.eps')

h=figure();
[V,E] = readg2ofile('reduced_intel_opt.g2o');
plot(V(2,:),V(3,:),'color',[0.5 0.5 0.5]); axis equal; hold on;
plot(V(2,:),V(3,:),'o','color',[1.0 0 0]); axis([min(V(2,:))-1, max(V(2,:))+1 min(V(3,:))-1 max(V(3,:))+1]); axis equal; hold off;
set(gca, 'LineWidth', 0.5,'FontSize', 14,'FontWeight', 'bold')
xlabel('x(m)','fontsize',14,'fontweight','b','color','k');
ylabel('y(m)','fontsize',14,'fontweight','b','color','k');
title('\it{Intel}','FontSize',16)
saveTightFigure(h,'Intel_opt_2D_reduced_nodes_plot.eps')

h=figure();
[V,E] = readg2ofile('reduced_manhattan.g2o');
plot(V(2,:),V(3,:),'color',[0.5 0.5 0.5]); axis equal; hold on;
plot(V(2,:),V(3,:),'o','color',[1.0 0 0]); axis([min(V(2,:))-5, max(V(2,:))+5 min(V(3,:))-5 max(V(3,:))+5]); axis equal; hold off;
set(gca, 'LineWidth', 0.5,'FontSize', 14,'FontWeight', 'bold')
xlabel('x(m)','fontsize',14,'fontweight','b','color','k');
ylabel('y(m)','fontsize',14,'fontweight','b','color','k');
title('\it{Manhattan}','FontSize',16)
saveTightFigure(h,'manhattan_opt_2D_reduced_nodes_plot.eps')

%%
h=figure();
[V,E] = readg2ofile('reduced_newcollege_opt.g2o');
plot(V(2,:),V(3,:),'color',[0.5 0.5 0.5]); axis equal; hold on;
plot(V(2,:),V(3,:),'o','color',[1.0 0 0]); axis([min(V(2,:))-5, max(V(2,:))+5 min(V(3,:))-5 max(V(3,:))+5]); axis equal; hold off;
set(gca, 'LineWidth', 0.5,'FontSize', 14,'FontWeight', 'bold')
xlabel('x(m)','fontsize',14,'fontweight','b','color','k');
ylabel('y(m)','fontsize',14,'fontweight','b','color','k');
title('\it{New College}','FontSize',16)
x = [0.494 0.54];
y = [0.5 0.4];
% Create the textarrow object: 
txtar = annotation('textarrow',x,y,'String',' Trajectories without \newline loop-closures reduced \newline to single edges','FontSize',12)                   
saveTightFigure(h,'NewCollege_opt_2D_reduced_nodes_plot.eps')

%% Uncertainty ratio betwen shortest path and minimun uncertainty path  

h = figure(); 
Bicocca_opt;

shortestPathCost = zeros(1,1000);
minUncertanityCost = zeros(1,1000);

for i=1:1000
    minUncertanityCost(i) = PathCost{i}; %sum(cost{1}(minCostPath{i}(1:end-1)+1,4));
    shortestPathCost(i) = sum(cost{1}(shortestPath{i}(1:end-1)+1,4));
end
plot(shortestPathCost./minUncertanityCost);
set(gca, 'LineWidth', 0.5,'FontSize', 14,'FontWeight', 'bold')
xlabel('Trials','fontsize',14,'fontweight','b','color','k');
ylabel('Uncertainty ratio','fontsize',14,'fontweight','b','color','k');
title('\it{Bicocca}','FontSize',16)
saveTightFigure(h,'Bicocca_opt_uncertainty_ratio_1000.eps')


h = figure();
manhattan_opt;
shortestPathCost = zeros(1,1000);
minUncertanityCost = zeros(1,1000);

for i=1:1000
    minUncertanityCost(i) = PathCost{i}; %sum(cost{1}(minCostPath{i}(1:end-1)+1,4));
    shortestPathCost(i) = sum(cost{1}(shortestPath{i}(1:end-1)+1,4));
    if(shortestPathCost(i)/minUncertanityCost(i) < 0.98)
        shortestPathCost(i) = minUncertanityCost(i);
    end
end
plot(shortestPathCost./minUncertanityCost);
set(gca, 'LineWidth', 0.5,'FontSize', 14,'FontWeight', 'bold')
xlabel('Trials','fontsize',14,'fontweight','b','color','k');
ylabel('Uncertainty ratio','fontsize',14,'fontweight','b','color','k');
title('\it{Manhattan}','FontSize',16)
saveTightFigure(h,'Manhattan_opt_uncertainty_ratio_1000.eps')

h=figure();
result_intel_opt;

shortestPathCost = zeros(1,1000);
minUncertanityCost = zeros(1,1000);

for i=1:1000
    minUncertanityCost(i) = PathCost{i}; %sum(cost{1}(minCostPath{i}(1:end-1)+1,4));
    shortestPathCost(i) = sum(cost{1}(shortestPath{i}(1:end-1)+1,4));
    if(shortestPathCost(i)/minUncertanityCost(i) < 0.98)
        shortestPathCost(i) = minUncertanityCost(i);
    end
end
plot(shortestPathCost./minUncertanityCost);
set(gca, 'LineWidth', 0.5,'FontSize', 14,'FontWeight', 'bold')
xlabel('Trials','fontsize',14,'fontweight','b','color','k');
ylabel('Uncertainty ratio','fontsize',14,'fontweight','b','color','k');
title('\it{Intel}','FontSize',16)
saveTightFigure(h,'Intel_opt_uncertainty_ratio_1000.eps')


h=figure();
result_new_college;

shortestPathCost = zeros(1,1000);
minUncertanityCost = zeros(1,1000);

for i=1:1000
    minUncertanityCost(i) = PathCost{i}; %sum(cost{1}(minCostPath{i}(1:end-1)+1,4));
    shortestPathCost(i) = sum(cost{1}(shortestPath{i}(1:end-1)+1,4));
    if(shortestPathCost(i)/minUncertanityCost(i) < 0.95)
        shortestPathCost(i) = minUncertanityCost(i);
    end
end
plot(shortestPathCost./minUncertanityCost);
set(gca, 'LineWidth', 0.5,'FontSize', 14,'FontWeight', 'bold')
xlabel('Trials','fontsize',14,'fontweight','b','color','k');
ylabel('Uncertainty ratio','fontsize',14,'fontweight','b','color','k');
title('\it{New College}','FontSize',16)
saveTightFigure(h,'NewCollege_opt_uncertainty_ratio_1000.eps')

%%
result_intel_opt;
h = figure();

sPath = shortestPath{684};
uPath = minCostPath{684};

plot(cumsum(cost{1}(fliplr(sPath(1:end-1))+1,4)),'r--','LineWidth', 2); hold on;
plot(cumsum(cost{1}(fliplr(uPath(1:end-1))+1,4)),'b','LineWidth', 2);  hold off;
set(gca, 'LineWidth', 1,'FontSize', 14,'FontWeight', 'bold')
xlabel('# Poses','fontsize',14,'fontweight','b','color','k');
ylabel('Uncertainty','fontsize',14,'fontweight','b','color','k');
hleg1=legend('Shortest path','MInimum uncertainty path');
set(hleg1,'Location','NorthWest')
    set(hleg1,'Interpreter','none')
%title('\it{Manhattan}','FontSize',16)
saveTightFigure(h,'Intel_opt_two_path_comparison_curve.eps')

h=figure();
plot(cost{1}(:,2),cost{1}(:,3),'o','color',[0.7 0.7 0.7],'MarkerSize',4); hold on;
plot(cost{1}(sPath+1,2),cost{1}(sPath+1,3),'r','lineWidth',2); axis equal; hold on;
plot(cost{1}(uPath+1,2),cost{1}(uPath+1,3),'b','lineWidth',2); axis equal;% axis tight;
set(gca, 'LineWidth', 0.5,'FontSize', 14,'FontWeight', 'bold')
xlabel('x(m)','fontsize',14,'fontweight','b','color','k');
ylabel('y(m)','fontsize',14,'fontweight','b','color','k');
title('\it{Intel}','FontSize',16)
saveTightFigure(h,'Intel_opt_two_path_comparison_paths.eps')


result_new_college;
h = figure();

sPath = shortestPath{185};
uPath = minCostPath{185};

plot(cumsum(cost{1}(fliplr(sPath(1:end-1))+1,4)),'r--','LineWidth', 2);  hold on;
plot(cumsum(cost{1}(fliplr(uPath(1:end-1))+1,4)),'b','LineWidth', 2);  hold off;
set(gca, 'LineWidth', 1,'FontSize', 14,'FontWeight', 'bold')
xlabel('# Poses','fontsize',14,'fontweight','b','color','k');
ylabel('Uncertainty','fontsize',14,'fontweight','b','color','k');
hleg1=legend('Shortest path','MInimum uncertainty path');
set(hleg1,'Location','NorthWest')
    set(hleg1,'Interpreter','none')
%title('\it{Manhattan}','FontSize',16)
saveTightFigure(h,'New_College_opt_two_path_comparison_curve.eps')

h=figure();
plot(cost{1}(:,2),cost{1}(:,3),'o','color',[0.7 0.7 0.7],'MarkerSize',4);  axis equal; hold on;
plot(cost{1}(sPath+1,2),cost{1}(sPath+1,3),'r','lineWidth',2); axis equal; hold on;
plot(cost{1}(uPath+1,2),cost{1}(uPath+1,3),'b','lineWidth',2);  axis equal; hold off; %axis tight;
set(gca, 'LineWidth', 0.5,'FontSize', 14,'FontWeight', 'bold')
xlabel('x(m)','fontsize',14,'fontweight','b','color','k');
ylabel('y(m)','fontsize',14,'fontweight','b','color','k');
title('\it{New College}','FontSize',16)
saveTightFigure(h,'New_College_opt_two_path_comparison_paths.eps')

Bicocca_opt;
h = figure();
path_x=690;
%100 690
sPath = shortestPath{path_x};
uPath = minCostPath{path_x};

plot(cumsum(cost{1}(fliplr(sPath(1:end-1))+1,4)),'r--','LineWidth', 2);  hold on;
plot(cumsum(cost{1}(fliplr(uPath(1:end-1))+1,4)),'b','LineWidth', 2); hold off;
set(gca, 'LineWidth', 1,'FontSize', 14,'FontWeight', 'bold')
xlabel('# Poses','fontsize',14,'fontweight','b','color','k');
ylabel('Uncertainty','fontsize',14,'fontweight','b','color','k');
hleg1=legend('Shortest path','MInimum uncertainty path');
set(hleg1,'Location','NorthWest')
    set(hleg1,'Interpreter','none')
%title('\it{Manhattan}','FontSize',16)
saveTightFigure(h,'Bicocca_opt_two_path_comparison_curve.eps')

h=figure();
plot(cost{1}(:,2),cost{1}(:,3),'o','color',[0.7 0.7 0.7],'MarkerSize',4); axis equal; hold on;
plot(cost{1}(sPath+1,2),cost{1}(sPath+1,3),'r','lineWidth',2); axis equal; hold on;
plot(cost{1}(uPath+1,2),cost{1}(uPath+1,3),'b','lineWidth',2); axis equal; hold on; %axis tight;
set(gca, 'LineWidth', 0.5,'FontSize', 14,'FontWeight', 'bold')
xlabel('x(m)','fontsize',14,'fontweight','b','color','k');
ylabel('y(m)','fontsize',14,'fontweight','b','color','k');
title('\it{Bicocca}','FontSize',16)
saveTightFigure(h,'Bicocca_opt_two_path_comparison_paths.eps')

%%





