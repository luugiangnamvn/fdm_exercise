%% ---PDEs Heate Equations with Crack-Nicolson---

clear 
close all
format long
clc

%% --------- Initial information ----------------
ax = 0;
bx = 1;
T = 1;
dx = 0.1;  %h
dt = 0.0025; %k
Nt = 400;
Nx = 10;
theta = 1/2;
number_mesh=4;
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
        
    
    

    %% ------Solve and calculate errors--------
    A = sparse(Nx-1,Nx-1);
    I = eye(Nx-1,Nx-1);
    J = sparse(Nx-1,Nx-1);
    F = zeros(Nx-1,1);
    u_dis = zeros(Nx+1,1);
    u_ex = zeros(Nx+1,1);
    error_H0(jj) = 0;
    error_L2(jj) = 0;
    for it = 2:Nt    
        for i = 1:Nx-1
            r = dt*(fb1(x(i+1),t(it))+c1(x(i+1),t(it)));
            r1 = dt*a1(x(i+1),t(it))/dx^2-dt/(2*dx)*(fa1(x(i+1),t(it))+b1(x(i+1),t(it)));
            r2 = -2*dt/dx^2*a1(x(i+1),t(it));
            r3 = dt*a1(x(i+1),t(it))/dx^2+dt/(2*dx)*(fa1(x(i+1),t(it))+b1(x(i+1),t(it)));
            if i == 1
                A(i,i) = r2;
                A(i,i+1) = r3;
            else
                if i == Nx-1
                    A(i,i) = r2;
                    A(i,i-1) = r1;
                else
                    A(i,i)= r2;
                    A(i,i+1) = r3;
                    A(i,i-1) = r1;
                end
            end
            J(i,i) = r;
        end
        G1 = I - theta*A;
        G2 = I + J + (1-theta)*A;
        G = (G1^(-1))*(G2);
        for i = 2:Nx
            r1 = dt*a1(x(i),t(it))/dx^2-dt/(2*dx)*(fa1(x(i),t(it))+b1(x(i),t(it)));
            r2 = -2*dt/dx^2*a1(x(i),t(it));
            r3 = dt*a1(x(i),t(it))/dx^2+dt/(2*dx)*(fa1(x(i),t(it))+b1(x(i),t(it)));
            if(i==2)
                F(i-1) = 0.5*f0(x(i),t(it))*dt+ 0.5*f0(x(i),t(it-1))*dt + r1*theta*g0(t(it)) + (1-theta)*r1*g0(t(it-1));
            elseif (i==Nx)
                F(i-1) = 0.5*f0(x(i),t(it))*dt+ 0.5*f0(x(i),t(it-1))*dt + r3*theta*g1(t(it)) + (1-theta)*r3*g1(t(it-1));
            else
                F(i-1) = 0.5*f0(x(i),t(it))*dt+ 0.5*f0(x(i),t(it-1))*dt;

            end
        end
        
        U = G*u0 + (G1^(-1))*F;
        u0 = U;
        
        u_dis(1) = g0(t(it));
        u_dis(Nx+1) = g1(t(it));
        
        for i = 2:Nx
                u_dis(i) = U(i-1);
        end
        for i = 1:Nx+1
            u_ex(i) = u_exact(x(i),t(it));
        end
        
        figure(jj)
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
plot(log(number_mesh_point), -log(error_H0)+2,'blue',...
    log(number_mesh_point), -log(error_L2), 'red',...
    log(number_mesh_point), 2*log(number_mesh_point)+4,'black',...
    log(number_mesh_point), log(number_mesh_point)+4,'green',...
    log(number_mesh_point), 3/2*log(number_mesh_point)+2,'magenta',...
    log(number_mesh_point), 1/2*log(number_mesh_point)+4,'yellow');
xlabel('Log(MeshPoint)');ylabel('-Log(Error)');
title('Errors');
legend('norm_{H_0}','norm_{L_2}','2x','x','3/2x','1/2x','Location','NorthEastOutside');  













