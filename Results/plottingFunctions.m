%%% Drawings for IROS 2012 Path Planning %%

colors = [0.1 0.6 0.1];

Bicocca_opt;

h = figure(1); plot3(cost{1}(:,2),cost{1}(:,3),cost{1}(:,4),'color',colors,'lineWidth',1.5);
grid on;
axis vis3d;
view([-38 22]);

%print(h,'Bicocca_25b_cost_3d.eps','-depsc2');

h = figure(2);
[V,E]=readg2ofile('Bicocca_opt.g2o');
g2o_plotGraphClusters(V,E,[]); axis equal; axis tight;


h=figure(3);
[V,E] = readg2ofile('reduced_graph.g2o');
plot(V(2,:),V(3,:),'color',[0.5 0.5 0.5]); hold on;
plot(V(2,:),V(3,:),'o','color',[1.0 0 0]); axis equal; axis tight;


h = figure(4); 
Bicocca_opt;

shortestPathCost = zeros(1,1000);
minUncertanityCost = zeros(1,1000);

for i=1:1000
    minUncertanityCost(i) = PathCost{i}; %sum(cost{1}(minCostPath{i}(1:end-1)+1,4));
    shortestPathCost(i) = sum(cost{1}(shortestPath{i}(1:end-1)+1,4));
end
plot(shortestPathCost./minUncertanityCost);

h=figure(5);
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


h = figure(6);

sPath = shortestPath{684};
uPath = minCostPath{684};

plot(cumsum(cost{1}(fliplr(sPath(1:end-1))+1,4)),'r-'); hold on;
plot(cumsum(cost{1}(fliplr(uPath(1:end-1))+1,4)),'g'); 


h=figure(7);
plot(cost{1}(:,2),cost{1}(:,3),'o','color',[0.7 0.7 0.7],'MarkerSize',4); hold on;
plot(cost{1}(sPath+1,2),cost{1}(sPath+1,3),'r','lineWidth',2); hold on;
plot(cost{1}(uPath+1,2),cost{1}(uPath+1,3),'g','lineWidth',2);
axis equal ; axis tight;


h = figure(10);
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




