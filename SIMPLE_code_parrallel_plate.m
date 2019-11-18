clear
close all

%space variables (ENTER)
Nx=300+20+3;           %number of columns (mm)
Ny=60+3;           %number of rows  (mm)

%discretization variables (ENTER)
dx = 1;          % x-grid size
dy = 1;          % y-grid size


%fluid density & viscosity (ENTER)
density = 1000;         %kg/m^3
viscosity = 1e-3;       %about pa.s
ratio = density/viscosity;       %ratio of density and dynamic viscosity 

%size and dimension of pressure and velocity (AUTO)
p = zeros(Ny,Nx);      %pressure
u = zeros(Ny+1,Nx+1);  %x-velocity  
v = zeros(Ny+1,Nx+1);  %y-velocity
residual = zeros(Ny*Nx,1); %residuals from continuity
dp = zeros(Ny*Nx,1);  %changes in pressures

%initial conditions (ENTER)
u = zeros(Ny+1,Nx+1)+7000;      %mm/s
%constant value of 0.1 (initializes velocity field)

%temporary variables (AUTO)
u1=u;  
v1=v;
dp1=zeros(Ny,Nx);
residual1=zeros(Ny,Nx);

%apply boundary conditions 
u1(2,2:Nx) = 0.0;
u1(Ny,2:Nx) = 0.0;
v1(2,2:Nx) = 0.0;
v1(Ny,2:Nx) = 0.0;

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

%removing numerical artifacts
u1(2,2:Nx-1) = 0.0;
u1(Ny,2:Nx-1) = 0.0;
v1(2,2:Nx-1) = 0.0;
v1(Ny,2:Nx-1) = 0.0;

%resultant velocity
vel = sqrt(u.^2+v.^2);




figure(1)
subplot(3,1,1);
contourf (p(2:Ny,13:Nx-13)*1e-3)  %plot pressure field
axis equal
xlabel ('x-dimension (mm)')
ylabel ('y-dimension (mm)')
title ('pressure field in (pa) in the parrallel plate channel')
colorbar


subplot(3,1,2);
contourf (vel(2:Ny,13:Nx-13))    %plot velocity magnitude field
axis equal
xlabel ('x-dimension (mm)')
ylabel ('y-dimension (mm)')
title ('velocity magnitude in (mm/s) with zero slip in the parrallel plate channel (boundary kept)')
colorbar

subplot(3,1,3);
contourf (vel(3:Ny-1,13:Nx-13))    %plot velocity magnitude field
axis equal
xlabel ('x-dimension (mm)')
ylabel ('y-dimension (mm)')
title ('velocity magnitude in (mm/s) with zero slip in the parrallel plate channel (boundary removed)')
colorbar

