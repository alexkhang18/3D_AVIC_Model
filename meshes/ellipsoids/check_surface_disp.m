clc
clear all
close all

%% 
v1 = readmatrix('vertices_ellipsoidal_cell.txt');
% v2 = readmatrix('vertices_ellipsoidal_cell2.txt');
centroid = mean(v1);

f = readmatrix('faces.txt');

l = 0.98;
F = [1/sqrt(l), 0, 0; 0, 1/sqrt(l), 0; 0, 0, l];

v2 = (F*(v1-centroid)')'+centroid;

d = v2-v1;

figure
trisurf(f,v1(:,1),v1(:,2),v1(:,3),'FaceAlpha',0.3); hold on;
trisurf(f,v2(:,1),v2(:,2),v2(:,3),'FaceAlpha',0.3); hold on;

quiver3(v1(:,1),v1(:,2),v1(:,3),d(:,1),d(:,2),d(:,3),'r')