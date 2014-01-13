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

figure(1);pcolor(x,y,p);title('pressure');colorbar;%shading interp;
figure(2);pcolor(x,y,sqrt(u.^2+v.^2));title('velocity');colorbar;%shading interp;
figure(3);pcolor(x,y,a);title('speed of sound');colorbar;%shading interp;
figure(4);pcolor(x,y,rho);title('density');colorbar;%shading interp;

% Plot grid

% x = reshape(x, 1, Nx*Ny);
% y = reshape(y, 1, Nx*Ny);
% figure(5);plot(x,y,'.');title('grid');axis equal;

% filename = 'interfaces.txt';
% col = 5;
% fid = fopen(filename);
% data = fscanf(fid,'%g\t',[col inf]);
% fclose(fid);
% 
% x1 = data(1,:);
% x2 = data(2,:);
% y1 = data(3,:);
% y2 = data(4,:);
% boundary = data(5,:);
% figure(5);hold on;
% for counter = 1:length(x1)
%     color = 'b';
%     if boundary(counter) == -1
%         color = 'r';
%     end
%     if boundary(counter) == -2
%         color = 'k';
%     end
%     if boundary(counter) == -3;
%         color = 'm';
%     end
%     plot([x1(counter),x2(counter)],[y1(counter),y2(counter)],color, 'LineWidth', 2.0);
% end
% y = -0.5:0.005:0.5;
% plot(-sqrt(0.5^2-y.^2),y,'g','LineWidth',2);