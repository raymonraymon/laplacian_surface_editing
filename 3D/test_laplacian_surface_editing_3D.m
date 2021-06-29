%test_laplacian_surface_editing_3D
clc;
close all;
clear ;
dbstop if error ;
warning('off');

[vertex,faces]=readOBJ('box.obj');
BI = [1 4 5];
BC = vertex(BI,:)+rand(length(BI),3);

U = laplacian_surface_editing_3D(vertex,faces,BI,BC);

drawMesh(vertex,faces);
hold on 
drawMesh(U,faces);