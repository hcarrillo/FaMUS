function [vertices, edges] = readg2ofile(filename)

fid = fopen(filename,'r');
x=0;
tag = 'sds';

vertices = zeros(4,1000);
edges = zeros(2,1000);

vertex_count = 1;
edges_count = 0;

while(~isempty(tag))
    tag = fscanf(fid,'%s',1);
    data = fscanf(fid,'%d %f %f %f');
    if(strcmp(tag,'VERTEX_SE2'))
        vertices(:,vertex_count)=data;
        vertex_count = vertex_count+ 1;
    elseif(strcmp(tag,'EDGE_SE2'))
        edges(:,edges_count+1)=data(1:2,1);
        edges_count = edges_count+1;
    end
end
vertices = vertices(:,1:vertex_count-1);
edges = edges(:,1:edges_count-1);