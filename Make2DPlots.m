% close all;clear;clc;

filename = 'result.txt';
col = 7;

fid = fopen(filename);
fgets(fid);
data = fscanf(fid,'%g\t',[2 1]);
Nx = data(1);
Ny = data(2);

fgets(fid);
data = fscanf(fid,'%g\t',[col inf]);
fclose(fid);

x = reshape(data(1,:), Nx, Ny);
y = reshape(data(2,:), Nx, Ny);
p = reshape(data(3,:), Nx, Ny);
u = reshape(data(4,:), Nx, Ny);
v = reshape(data(5,:), Nx, Ny);
a = reshape(data(6,:), Nx, Ny);
rho = reshape(data(7,:), Nx, Ny);

figure(1);pcolor(x,y,p);title('pressure');colorbar;shading interp;
figure(2);pcolor(x,y,sqrt(u.^2+v.^2));title('velocity');colorbar;shading interp;
figure(3);pcolor(x,y,a);title('speed of sound');colorbar;shading interp;
figure(4);pcolor(x,y,rho);title('density');colorbar;shading interp;