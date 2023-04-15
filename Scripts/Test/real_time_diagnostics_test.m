clear
close all

% r = 50;
% n = 1000;
% X = repmat(linspace(-r,r,n), n, 1);
% Y = rot90(X);
% C = X.^2 + Y.^2 <= r^2;
% imagesc(C)
% set(gcf,'Renderer','painters');
% movegui(gcf,'center');
% axis image
% 
% 
% nc = 10;
% MC = repmat(C, nc, nc); % multiple circles
% imagesc(MC)
% set(gcf,'Renderer','painters');
% movegui(gcf,'center');
% axis image

delta = 10;

x_begin = -delta/2;
y_begin = -delta/2;

Nx = 23;
Ny = 18;

x_end = (Nx-1)*delta;
y_end = (Ny-1)*delta;

x_pos = x_begin:delta:x_end;
y_pos = y_begin:delta:y_end;

colors = {'red', 'yellow', 'green'};  % Cell array of strings
A = ones(Ny, Nx);  % Sample data
cM = colors(A);   % Ny-by-Nx cell array of strings

set(gcf,'Renderer','painters');
for nx = 1:Nx
    for ny = 1:Ny
        rectangle('Position',[x_pos(nx),y_pos(ny),delta,delta],'FaceColor',cM{ny,nx});
    end
end
axis equal

for nx = 1:Nx
    for ny = 1:Ny
        pause(0.1)
        A(ny,nx) = 2;
        cM = colors(A);
        rectangle('Position',[x_pos(nx),y_pos(ny),delta,delta],'FaceColor',cM{ny,nx});
        axis equal
        drawnow
        pause(0.1)
        A(ny,nx) = 3;
        cM = colors(A);
        rectangle('Position',[x_pos(nx),y_pos(ny),delta,delta],'FaceColor',cM{ny,nx});
        axis equal
        drawnow
    end
    
end



