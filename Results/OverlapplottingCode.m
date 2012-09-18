result_new_college;
for i=1:length(minCostPath)
    figure(1); plot(cost(:,2),cost(:,3),'.'); hold on;  plot(cost(minCostPath{i}+1,2),cost(minCostPath{i}+1,3),'r'); hold off;
    drawnow;
    pause(0.33);
end

i=28;

plot(cost{1}(:,2),cost{1}(:,3),'x'); hold on;  plot(cost{1}(shortestPath{i}+1,2),cost{1}(shortestPath{i}+1,3),'r');
hold on; 
plot(cost{1}(minCostPath{i-1}+1,2),cost{1}(minCostPath{i-1}+1,3),'g')


%%% Overlap 
commonVertices = zeros(1,1000);
for i=1:1000
    a = minCostPath{i};
    b = shortestPath{i};
    if(length(b)==2)
        commonVertices(i) = 1;
        continue;
    end
    common = intersect(a,b);
    commonVertices(i) = (length(common)-2)./(length(b)-2);
end
percentOverlap = mean(commonVertices)  
allsame = length(find(commonVertices==1))
allDifferent = length(find(commonVertices==0))


dopt = cost{1}(:,4);
trace = cost{1}(:,5);
dopt/max(dopt(:))
dopt_n = dopt/max(dopt(:));
trace_n = trace/max(trace(:));
plot3(cost{1}(:,2),cost{1}(:,3),dopt_n);
hold on; plot3(cost{1}(:,2),cost{1}(:,3),trace_n,'r');

