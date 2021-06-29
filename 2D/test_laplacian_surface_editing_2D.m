clc;
close all;
clear ;
dbstop if error ;
warning('off');
% assume xy is n x 2

% uncomment the curve you would like to try
%xy = load('unicorn.txt');
xy = load('T.txt');
anchors = [99:101 240:242 171 172];
anch_pos = xy(anchors,:)+rand((length(anchors)),2);

U = laplacian_surface_editing_2D(xy,anchors,anch_pos);

figure(1);
plot(xy(:,1),xy(:,2),'g-',...
     U(:,1),U(:,2),'-b',...
     U(anchors,1),U(anchors,2),'*k',...
     anch_pos(:,1),anch_pos(:,2), 'mo',...
     U(anchors(length(anchors)),1),U(anchors(length(anchors)),2),'*r');
axis equal;
title('Editing result');