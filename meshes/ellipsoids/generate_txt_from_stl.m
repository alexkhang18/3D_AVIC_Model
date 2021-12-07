clc
close all
clear all

%% creates .stl
TR = stlread("ellipsoidal_cell2.stl"); 

vertices = TR.Points;
faces = TR.ConnectivityList;

writematrix(vertices,'vertices_ellipsoidal_cell2.txt')
writematrix(faces,'faces_ellipsoidal_cell2.txt')
