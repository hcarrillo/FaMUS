 %
 % plottingFunctions.m
 %
 %  Created on: Feb 6, 2012
 %  Author: Henry
 %
 
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
%% Comparison of two paths
Bicocca_opt;
h = figure();
path_x=1;
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
