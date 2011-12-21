%% Filename :movingplate2d.m
% Written by : amneet bhalla on taylor@mech.northwestern.edu
clear all;
clc;
%% Mesh parameters.
Lx = 32.0; Ly = 22.0;
Nx = 64*4*4; Ny = 64*4*4;
dx = Lx/Nx; dy = Ly/Ny;

%%plate parameters.
Length_plate = 1.0;
num_pts_y = ceil(Length_plate/dy);

%% write out the points.
idx = 0;
    
    for j = 1:num_pts_y
        y = (j-0.5)*dy ;
       
        idx = idx+1;
        X_array(idx) = 0.0; Y_array(idx) = y;
    end

XCOM = sum(X_array)/length(X_array);
YCOM = sum(Y_array)/length(Y_array);

X_array = X_array - XCOM;
%Y_array = Y_array - YCOM;
%% plot the plate
plot(X_array, Y_array, '.');
%axis 'equal'

%% write the coordinates in the file
fid = fopen('../movingplate2d.vertex','wt');
fprintf(fid,'%d\n', length(X_array));

for i = 1:length(X_array)
    fprintf(fid,'%f\t%f\n',X_array(i),Y_array(i));
end

