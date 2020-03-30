clear
close all
%space variables (ENTER)
Nx=100;           %number of columns
Ny=50;           %number of rows 
%discretization variables (ENTER)
dx = 1;          % x-grid size
dy = 1;          % y-grid size
%fluid density & viscosity (ENTER)
density = 1000;
viscosity = 1;
ratio = density/viscosity;       %ratio of density and dynamic viscosity 
%size and dimension of pressure and velocity (AUTO)
p = zeros(Ny,Nx);      %pressure
u = zeros(Ny+1,Nx+1);  %x-velocity  
v = zeros(Ny+1,Nx+1);  %y-velocity
residual = zeros(Ny*Nx,1); %residuals from continuity
dp = zeros(Ny*Nx,1);  %changes in pressures
%initial conditions (ENTER)
u = zeros(Ny+1,Nx+1)+1;
%constant value of 0.1 (initializes velocity field)
%temporary variables (AUTO)
u1=u;  
v1=v;
dp1=zeros(Ny,Nx);
residual1=zeros(Ny,Nx);
%apply boundary conditions

% u1(Ny/2-Ny/5+2:Ny/2+Ny/5,Nx/2-Nx/5+2:Nx/2+Nx/5) = 0.0; %for rectangle/square                  
% v1(Ny/2-Ny/5+2:Ny/2+Ny/5,Nx/2-Nx/5+2:Nx/2+Nx/5) = 0.0;

% radius = 0.25*(min(Nx,Ny)); %for circles
% center = [Nx/2,Ny/2];
% for i=1:Ny         
%     for j=1:Nx
%         distance = sqrt((i-center(2))^2 + (j-center(1))^2);
%         if (distance <= radius)
%         u1(i,j) = 0;
%         v1(i,j) = 0;
%         end
%     end
% end
%spy(u1)

% A = imread('image_star.bmp');
% A = imread('image_aero.bmp');
A = imread('image_wavey.bmp');
%spy(A)

for i=1:Ny         
    for j=1:Nx
        if (A(Ny-i+1,j) == 0)
        u1(i,j) = 0;
        v1(i,j) = 0;
        end
    end
end
%spy(u1)

        
%timestep value, relaxation factor, number of iterations (ENTER)
dt=1;             
relaxation_factor=1;            
total_iterations=1;            
residual_max = zeros(total_iterations,1);
%check CFL criteria (CHECK!)
CFL_x = max(max(u))*dt/dx;
CFL_y = max(max(u))*dt/dy;
%calculate sparse matrix (AUTO)
J_a = 2*(1/dx^2+1/dy^2);
J_b = -1/dy^2;
J_c = -1/dx^2;
J=spalloc(Nx*Ny,Nx*Ny,(Nx-2)*(Ny-2)*4+Nx*Ny);
for i=1:Nx*Ny-1
    J(i,i+1)=J_b;
end
for i=2:Nx*Ny
    J(i,i-1)=J_b;
end
for i=1:Nx*Ny-Nx
    J(i,i+Nx)=J_c;
end
for i=1:Nx*Ny-Nx
    J(i+Nx,i)=J_c;
end
for i=1:Ny-1
J(i*Nx+1,i*Nx)=0;
J(i*Nx,i*Nx+1)=0;
end
for i=1:Nx*Ny
    J(i,i)=-sum(J(i,:));
end
%spy(J)   %for checking sparse matrix
for j=1:Ny
    for i=1:Nx
        
        residual1(j,i)=-((u1(j,i+1)-u1(j,i))/dx+(v1(j+1,i)-v1(j,i))/dy)/dt;    %calculate residuals from continuity
   
    end
end
for j=1:Ny
    for i=1:Nx
        
        residual(Nx*(j-1)+i,1)=residual1(j,i);                          %converting residual from a matrix to a vector
    
    end
end
dp=J\residual;                                              %changes in pressure field
for j=1:Ny
    for i=1:Nx
        dp1(j,i)=dp(Nx*(j-1)+i,1);                             %converitng changes in pressure field from a vector to a matrix
    end
end
for j=2:Ny
    for i=2:Nx
        u1(j,i)=u1(j,i)+relaxation_factor*(dp1(j,i-1)-dp1(j,i))*dt/dx;                 %u velocity correction
    end
end
for j=2:Ny
    for i=2:Nx+1
        v1(j,i)=v1(j,i)+relaxation_factor*(dp1(j-1,i-1)-dp1(j,i-1))*dt/dy;             %v velocity correction
    end
end
p = p + relaxation_factor*dp1;                                      %pressure field correction
u = u1;
v = v1;
% 
figure()
%contourf (p)  %plot pressure field
hold on
UU = u(2:Ny-1,3:Nx);  %select u velocity field (adjust for staggered grid)
VV = v(2:Ny-1,3:Nx);  %select v velocity field (adjust for staggered grid)
[X,Y]=meshgrid(2:1:Nx-1,2:1:Ny-1);   %vector plot
q=quiver(X,Y,UU,VV,1);
q.Color = 'black';
axis equal;
%draw a square
% v_ = [Nx/2-Nx/5+1 Ny/2-Ny/5+1.5;  Nx/2-Nx/5+1 Ny/2+Ny/5+0.5; Nx/2+Nx/5 Ny/2-Ny/5+1.5; Nx/2+Nx/5 Ny/2+Ny/5+.5];
% f_ = [2 1 3 4];
% p_=patch('Faces',f_,'Vertices',v_,'FaceColor','red');
% p_.EdgeColor='none';
% p_.FaceColor='white';
%draw a circle (requires MATLAB image processing toolbox)
% center1 = center;
% center1(1)=center1(1) - 0.5;
% radius1 = radius + 0.25;
% h = drawcircle('Center',center1,'Radius',radius1,'Color','white','FaceAlpha',1);
xlabel ('x-dimension')
ylabel ('y-dimension')
title ('Pressure and velocity field around the circle')
colorbar