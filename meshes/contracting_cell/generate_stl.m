clc
clear all
close all

%%
vertices = readmatrix('predicted_endo1_vertices.txt');
faces = readmatrix('faces.txt');

TR = triangulation(faces,vertices(:,1),vertices(:,2),vertices(:,3));

stlwrite(TR,'predicted_endo1.stl')