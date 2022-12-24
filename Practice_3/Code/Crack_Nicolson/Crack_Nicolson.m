%% ---PDEs Heate Equations with Crack-Nicolson---

clear 
close all
format long
clc

%% --------- Initial information ----------------
ax = 0;
bx = 1;
T = 1;
dx = 0.2;  %h
dt = 0.01; %k
Nt = 100;
Nx = 5;
kap = 1; % Example 1
% kap = 1/2; % Example 2

number_mesh=6;
number_mesh_point=zeros(number_mesh,1);
error_H0=zeros(number_mesh,1);
error_L2=zeros(number_mesh,1);

%% ----------------Create vectors---------------

for jj = 1:number_mesh
    number_mesh_point(jj)=Nx;
    x= zeros(Nx+1,1);
    for i = 1:Nx+1
        x(i) = ax + (bx - ax)*(i-1)/Nx;
    end
    
    t = zeros(Nt+1,1);
    for i = 1:Nt
        t(i) = dt*(i-1);
    end
    
    u0 = zeros(Nx-1,1);
    for i = 2:Nx
        u0(i-1) = fu0(x(i));
    end
    
    %% -------------Create matric A-------------
    A = sparse(Nx-1,Nx-1);
    I = eye(Nx-1,Nx-1);
    for i = 1:Nx-1
        if i == 1
            A(i,i) = -2;
            A(i,i+1) = 1;
        else
            if i == Nx-1
                A(i,i) = -2;
                A(i,i-1) = 1;
            else
                A(i,i)= -2;
                A(i,i+1) = 1;
                A(i,i-1) = 1;
            end
        end
    end
    
    A = A*(kap*dt/dx^2);
    %% -------------Crack_Nicolson-------------
    G1 = I - 0.5*A;
    G2 = I + 0.5*A;
    G = (G1^(-1))*(G2);

    %% ------Solve and calculate errors--------
    F = zeros(Nx-1,1);
    u_dis = zeros(Nx+1,1);
    u_ex = zeros(Nx+1,1);
    error_H0(jj) = 0;
    error_L2(jj) = 0;
    for it = 2:Nt    
        for i = 2:Nx
            F(i-1) = 0.5*f0(x(i),t(it))*dt+ 0.5*f0(x(i),t(it-1))*dt;
        end
        
        U = G*u0 + (G1^(-1))*F;
        u0 = U;
    
        figure(jj)
        for i = 2:Nx
            u_dis(i) = U(i-1);
        end
        for i = 2:Nx
            u_ex(i) = u_exact(x(i),t(it));
        end
        plot(x,u_dis,'r',x,u_ex,'b');

        for i = 2:Nx
            error_H0 = error_H0 + (((u_dis(i+1)-u_ex(i+1))-(u_dis(i)-u_ex(i))))^2*dx/dt;
        end
        
        for i = 2:Nx
            error_L2 = error_L2 + (u_dis(i)-u_ex(i))^2*dx*dt;        
        end
    end
    %% ------------------Error------------------

    error_H0(jj) = (error_H0(jj))^(1/2);
    error_L2(jj) = (error_L2(jj))^(1/2);

    dx = dx/2;
    Nx = Nx*2;
    dt = dt/4;
    Nt = Nt*4;
    
end
%%
error_H0
error_L2

%%
figure(number_mesh+1)
plot(log(number_mesh_point), -log(error_H0),'blue',...
    log(number_mesh_point), -log(error_L2), 'red',...
    log(number_mesh_point), 2*log(number_mesh_point),'black',...
    log(number_mesh_point), log(number_mesh_point),'green',...
    log(number_mesh_point), 3/2*log(number_mesh_point),'yellow');
xlabel('Log(MeshPoint)');ylabel('-Log(Error)');
title('Errors');
legend('norm_{H_0}','norm_{L_2}','2x','x','3/2x','Location','NorthEastOutside');  













