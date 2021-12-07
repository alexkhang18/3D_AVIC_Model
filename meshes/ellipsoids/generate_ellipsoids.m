clc
close all
clear all

%%
a = 70/4;
b = 70/4;
c = 70/4;

% creates ellipsoid centers at the origin
[x,y,z] = ellipsoid(149.95/2,149.95/2,140/2,c,b,a,25);

% creates triangular mesh
fvc = surf2patch(x,y,z,'triangles');


vertices_n = fvc.vertices;
faces_n = fvc.faces;

writematrix(vertices_n,'vertices_ellipsoidal_cell2.txt')
writematrix(faces_n,'faces.txt')

% writes an STL
TR = triangulation(faces_n, vertices_n);
stlwrite(TR,"ellipsoidal_cell2.stl")

a = 1.3;
F = [sqrt(1/a), 0, 0; 0, sqrt(1/a), 0; 0, 0, a];

centroid = mean(vertices_n);
vertices_n2 = (F*(vertices_n-centroid)')' + centroid;
writematrix(vertices_n2,'vertices_ellipsoidal_cell.txt')

% writes an STL
TR = triangulation(faces_n, vertices_n2);
stlwrite(TR,"ellipsoidal_cell.stl")

vertices_n3 = vertices_n2 - centroid;
[azimuth,elevation,r] = cart2sph(vertices_n3(:,1),vertices_n3(:,2),vertices_n3(:,3));
[x,y,z] = sph2cart(azimuth,elevation,r*0.5);
vertices_n3 = [x,y,z] + centroid;
TR = triangulation(faces_n, vertices_n3);
stlwrite(TR,"ellipsoidal_nucleus.stl")

% figure
% trimesh(faces_n,vertices_n(:,1),vertices_n(:,2),vertices_n(:,3),'FaceColor','b'); hold on;
% trimesh(faces_n,vertices_n2(:,1),vertices_n2(:,2),vertices_n2(:,3),'FaceColor','r'); hold off;

% xlim([0 150])
% ylim([0 150])
% zlim([0 140])