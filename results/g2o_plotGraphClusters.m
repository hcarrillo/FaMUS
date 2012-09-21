function g2o_plotGraph(VertexList, EdgeList,membership, odometryColor, height)

if(~exist('odometryColor'))
    odometryColor = [0.7 0.7 0.7];
end

if(~exist('height'))
    height = zeros(1,length(VertexList));
end

line_width_size =1;
% to take care of zero-based indexing in the loops
%VertexList = [VertexList(:,1) VertexList];

num_loopClosures = length(membership);
membership = membership+1;
colors = rand(num_loopClosures,3);
colors(1,:) = [0.5,0.5,0.5];
colors(2,:) = [0.5,0.5,0.5];

[R C] = size(VertexList);

%plot3(VertexList(2,:),VertexList(3,:),'gr.');
beginVertex = VertexList(:,EdgeList(1,:)+1);
endVertex = VertexList(:,EdgeList(2,:)+1);


%loop closures
start = length(EdgeList)-num_loopClosures+1;
for i=start:length(EdgeList)
    if(membership(i-start+1)==1)
       %line([beginVertex(2,i) endVertex(2,i)],[beginVertex(3,i) endVertex(3,i)],[height(beginVertex(1,i)+1) height(endVertex(1,i)+1)],'Color',[1.0 0.2 0.2],'LineWidth',1.0)
    end
end
for i=start:length(EdgeList)
    if(membership(i-start+1)==2)
       line([beginVertex(2,i) endVertex(2,i)] ,[beginVertex(3,i) endVertex(3,i)],[height(beginVertex(1,i)+1) height(endVertex(1,i)+1)],'Color',[0 0.6 0.1],'LineWidth',line_width_size)
    end
end

for i=1:1:length(EdgeList)-num_loopClosures
    if(height(beginVertex(1,i)+1)/50 == 0)
        odometryColor = 'blue';
    elseif(height(beginVertex(1,i)+1)/50 == 1)
        odometryColor = 'green';
    elseif(height(beginVertex(1,i)+1)/50 == 2)
        odometryColor = 'red';
    elseif(height(beginVertex(1,i)+1)/50 == 3)
        odometryColor = 'cyan';
    elseif(height(beginVertex(1,i)+1)/50 == 4)
        odometryColor = 'magenta';
    end
line([beginVertex(2,i) endVertex(2,i) ],...
        [beginVertex(3,i) endVertex(3,i) ],...
        [height(beginVertex(1,i)+1) height(endVertex(1,i)+1)],'LineWidth',line_width_size, 'Color',odometryColor)
end
end

% for i=start:length(EdgeList)
%     if(membership(i-start+1)==2)
%        line([ beginVertex(2,i) endVertex(2,i)],[ beginVertex(3,i) endVertex(3,i)],'Color','blue','MarkerSize',4)
%     end
% end





