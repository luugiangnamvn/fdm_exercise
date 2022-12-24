%% Solve equation -u''(x)=f(x) with the Dirichlet boundary condition 
clear all
close all
clc
%% Initial informations
ax=0.0;
bx=1.0;
cases=2;
al = exact_solution(0,cases);
be = exact_solution(1,cases);
N=100;% number of mesh points of first mesh
number_mesh=6;
number_mesh_point=zeros(number_mesh,1);
norm_max=zeros(number_mesh,1);
norm_l2=zeros(number_mesh,1);
norm_maxh1=zeros(number_mesh,1);
norm_h1=zeros(number_mesh,1);
%% Solve discrite solution and refine mesh

for inumber_mesh=1:number_mesh
    
    number_mesh_point(inumber_mesh)=N;
    del = zeros(N,1);
%% Create mesh point    
    x=zeros(N+1,1);
    for i=1:N+1
        x(i)=1-cos(pi*(i-1)/(2*N));
    end
        
 %% Create matrix A   
    A=sparse(N-1,N-1);
%A = zeros(N-1,N-1);
    for i=1:N-1
        a = (4*(x(i+1)-x(i)))/((x(i+2)-x(i))*((x(i+1)-x(i))^2+(x(i+2)-x(i+1))^2));
        b = -4/((x(i+1)-x(i))^2+(x(i+2)-x(i+1))^2);
        c = (4*(x(i+2)-x(i+1)))/((x(i+2)-x(i))*((x(i+1)-x(i))^2+(x(i+2)-x(i+1))^2));
        if (i==1)
            A(i,i)=-b;
            A(i,i+1)=-a;
        elseif(i==N-1)
            A(i,i)=-b;
            A(i,i-1)=-c;
        else
            A(i,i+1)=-a;
            A(i,i)=-b;
            A(i,i-1)=-c;
        end
    end
    
%% Create vector b    
    b=zeros(N-1,1);
    for i=1:N-1
        a = (4*(x(i+1)-x(i)))/((x(i+2)-x(i))*((x(i+1)-x(i))^2+(x(i+2)-x(i+1))^2));
        c = (4*(x(i+2)-x(i+1)))/((x(i+2)-x(i))*((x(i+1)-x(i))^2+(x(i+2)-x(i+1))^2));
        if(i==1)
            b(i)=functionf(x(i,1),cases) + c*al;
        elseif(i==N-1)
            b(i)=functionf(x(i,1),cases) + a*be;
        else
            b(i)=functionf(x(i,1),cases);
        end
    end
%% Solve discrete solution
    u=A\b;
%% Get exact solution    
    u_ex=zeros(N+1,1);
    for i=1:N+1
        u_ex(i)=exact_solution(x(i,1),cases);
    end
%% Create discrete solution with boundary 
    u_dis=zeros(N+1,1);
    for i=1:N+1
        if (i==1)
            u_dis(i)=al;
        elseif(i==N+1)
            u_dis(i)=be;
        else
            u_dis(i)=u(i-1,1);
        end
    end
    
    
%% Calculate the error on L^infinity
    norm_max(inumber_mesh)=0.0;
    for i=1:N+1
        if (abs(u_dis(i)-u_ex(i)) > norm_max(inumber_mesh))
            norm_max(inumber_mesh)=abs(u_dis(i)-u_ex(i));
        end
    end
    
    norm_max(inumber_mesh)
%%  Calculate the error on L^2 

    norm_l2(inumber_mesh)=0;
    for i=1:N+1
        if (i==1)
            norm_l2(inumber_mesh)=norm_l2(inumber_mesh)+...
                (u_dis(i)-u_ex(i))^2*(x(i+1)-x(i));
        elseif(i==N+1)
            norm_l2(inumber_mesh)=norm_l2(inumber_mesh)+...
                (u_dis(i)-u_ex(i))^2*(x(i)-x(i-1));
        else
            norm_l2(inumber_mesh)=norm_l2(inumber_mesh)+...
                (u_dis(i)-u_ex(i))^2*((x(i-1)+x(i+1)-2*x(i))/2);
        end
    end
    
    norm_l2(inumber_mesh)=(norm_l2(inumber_mesh))^(1/2);
    norm_l2(inumber_mesh)
%% Calculate the error on maxH1    

    norm_maxh1(inumber_mesh)=0;
    for i=1:N
        if (abs(((u_dis(i+1)-u_ex(i+1))-(u_dis(i)-u_ex(i)))/(x(i+1)-x(i))) ...
                > norm_maxh1(inumber_mesh))
            norm_maxh1(inumber_mesh)=abs(((u_dis(i+1)-u_ex(i+1))-...
                (u_dis(i)-u_ex(i)))/(x(i+1)-x(i)));
        end
    end
    norm_maxh1(inumber_mesh)

%% Calculate the error on H1

    norm_h1(inumber_mesh)=0;
    for i=1:N
        norm_h1(inumber_mesh)=norm_h1(inumber_mesh)+...
            (((u_dis(i+1)-u_ex(i+1))-(u_dis(i)-u_ex(i)))/(x(i+1)-x(i)))^2*...
            (x(i+1)-x(i));
    end
    norm_h1(inumber_mesh)=(norm_h1(inumber_mesh))^(1/2);
%% Figure exact and dicrete solutions    
    figure
    plot(x,u_ex,'blue', x,u_dis,'red');
    xlabel('x');ylabel('value');
    title('Comparison between exact and discrete solutions with N=');
%     ylim([-0.02 0.12]);
    legend('Exact solution','Discrete solution');
    
%% Refine mesh (increse mesh point)    
    N=2*N;
end



%% Figure for errors respect to number of mesh point
figure
plot(log(number_mesh_point), -log(norm_max),'blue', log(number_mesh_point), -log(norm_l2), 'red',...
    log(number_mesh_point), -log(norm_maxh1), 'cyan', log(number_mesh_point), -log(norm_h1), 'magenta',...
    log(number_mesh_point), 2*log(number_mesh_point),'black',...
    log(number_mesh_point), 3/2*log(number_mesh_point),'green',...
    log(number_mesh_point), log(number_mesh_point),'yellow');
xlabel('Log(MeshPoint)');ylabel('-Log(Error)');
title('Errors');
legend('norm_max','norm_l2','norm_maxh1','norm_h1','2x','3/2x','x','Location','NorthEastOutside');  
