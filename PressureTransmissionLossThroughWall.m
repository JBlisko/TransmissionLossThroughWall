%Author: Jonathan Blisko
%Date: 11/29/20

clear all
close all
clc

%% Testing
road={'Road'};
room={'Room'};
wall={'Wall'};

figure
line([0.5 0.5],[0 0.5],'Color','k','LineStyle','-.')
line([0 0.5],[0.5 0.5],'Color','k','LineStyle','-.')
line([0.45 0.45],[0 0.45],'Color','k','LineStyle','-.')
line([0 0.45],[0.45 0.45],'Color','k','LineStyle','-.')
line([1 1],[0 1],'Color','k','LineStyle','--')
line([0 1],[1 1],'Color','k','LineStyle','--')
xlim([0 1])
ylim([0 1])
pbaspect([1 1 1])
text(0.45,0.95,road,'Interpreter','Latex')
text(0.9,0.5,road,'Interpreter','Latex')
text(0.2,0.25,room,'Interpreter','Latex')
text(0.52,0.25,wall,'Interpreter','Latex')
title('Reference Diagram','Interpreter','Latex','Color','white')
xlabel('X - Direction, \textit{x/L}','Interpreter','Latex','Color','w')
ylabel('Y - Direction, \textit{y/L}','Interpreter','Latex','Color','w')
set(gca,'xcolor','w')
set(gca,'ycolor','w')

%% Initialize Variables
% Constants
c=343;
SteadyConvergence=1e-6;
k=1;
k_2=1;
max_it=50000;
diff=1;
C=0.69;
p_0=0.1;
P_ref=20e-6;
rho=1.225; %kg/m^3
omega=2*pi*1000;

% X-Direction
L=1;
deltax=0.01;
nx=L/deltax+1;
x=0:deltax:L;

% Y-Direction
H=L;
deltay=deltax;
ny=H/deltay+1;
y=0:deltay:H;

% Time
deltat=deltax/(sqrt(2)*c);

% Wall Properties
h=0.05;
rho_w=1000; %kg/m^3

%Transmission Loss Through Wall (Constant Thickness)
TL=10*log10(1+(rho_w*h*omega/(2*rho*c))^2);

%% Initial Conditions
p=zeros(nx,ny)+0.01;
spl=zeros(nx,ny);

% Boundary Constant
% a=0 is pure sine, a=1 is constant
a=0;
% Number of Peaks on Each Boundary Face
b=6;

for i=1:nx
    p(i,ny)=(1-a)*abs(sin(b*pi*x(i)))+a;
    p(nx,i)=(1-a)*abs(sin(b*pi*y(i)))+a;
end

figure
[X,Y]=meshgrid(x,y);
contourf(Y,X,p)
pbaspect([1 1 1])

%% Governing Equations
while (diff>=SteadyConvergence && k<max_it)
    p_old=p;
    
    for i=nx-1:-1:51
        % Outside to Right
        for j=ny-1:-1:2
            p(i,j)=2*p(i,j)-p_old(i,j)+C^2*(p(i+1,j)+p(i-1,j)+p(i,j+1)+p(i,j-1)-4*p(i,j));
            spl(i,j)=10*log10((p_0*p(i,j))^2/(2*P_ref^2));
        end
        p(i,1)=p(i,2);
        spl(i,1)=spl(i,2);
    end
    
    for i=50:-1:46
        % Outside Above Side Wall
        for j=ny-1:-1:51
            p(i,j)=2*p(i,j)-p_old(i,j)+C^2*(p(i+1,j)+p(i-1,j)+p(i,j+1)+p(i,j-1)-4*p(i,j));
            spl(i,j)=10*log10((p_0*p(i,j))^2/(2*P_ref^2));
        end
    end
    
    for i=45:-1:2
        % Outside Above Room
        for j=ny-1:-1:51
            p(i,j)=2*p(i,j)-p_old(i,j)+C^2*(p(i+1,j)+p(i-1,j)+p(i,j+1)+p(i,j-1)-4*p(i,j));
            spl(i,j)=10*log10((p_0*p(i,j))^2/(2*P_ref^2));
            p(1,j)=p(2,j);
            spl(1,j)=spl(2,j);
        end
        
        % Room
        for j=45:-1:2
            p(i,j)=2*p(i,j)-p_old(i,j)+C^2*(p(i+1,j)+p(i-1,j)+p(i,j+1)+p(i,j-1)-4*p(i,j));
            spl(i,j)=10*log10((p_0*p(i,j))^2/(2*P_ref^2));
            p(1,j)=p(2,j);
            spl(1,j)=spl(2,j);
        end
        
        p(i,1)=p(i,2);
        spl(i,1)=spl(i,2);
    end
    
    % Upper Wall
    spl(1:45,45)=spl(1:45,51)-TL;
    p(1:45,45)=sqrt(10.^(spl(1:45,45)/10)*P_ref^2)/p_0;
    
    % Side Wall
    spl(45,1:45)=spl(51,1:45)-TL;
    p(45,1:45)=sqrt(10.^(spl(45,1:45)/10)*P_ref^2)/p_0;
    
    % Advance While Conditions
    diff=max(max(abs(p-p_old)));
    k=k+1;
end

%% Plot
figure
[X,Y]=meshgrid(x,y);
contourf(Y,X,p)
title('Normalized Pressure, \textit{${p/p_0}$}','Interpreter','Latex','Color','white')
xlabel('X - Direction, \textit{x/L}','Interpreter','Latex','Color','w')
ylabel('Y - Direction, \textit{y/L}','Interpreter','Latex','Color','w')
colorbar
colormap('parula')
set(gcf,'RendererMode','manual')
xlim([0 1])
ylim([0 1])
pbaspect([1 1 1])

figure
[X,Y]=meshgrid(x,y);
contourf(Y,X,spl,35)
title('Sound Pressure Level (dB)','Interpreter','Latex','Color','white')
xlabel('X - Direction, \textit{x/L}','Interpreter','Latex','Color','w')
ylabel('Y - Direction, \textit{y/L}','Interpreter','Latex','Color','w')
colorbar
colormap('parula')
set(gcf,'RendererMode','manual')
xlim([0 1])
ylim([0 1])
pbaspect([1 1 1])
hold on
line([0.5 0.5],[0 0.5],'Color','k','LineStyle','-.')
line([0 0.5],[0.5 0.5],'Color','k','LineStyle','-.')
line([0.45 0.45],[0 0.45],'Color','k','LineStyle','-.')
line([0 0.45],[0.45 0.45],'Color','k','LineStyle','-.')
line([1 1],[0 1],'Color','k','LineStyle','--')
line([0 1],[1 1],'Color','k','LineStyle','--')
text(0.45,0.95,road,'Interpreter','Latex')
text(0.88,0.5,road,'Interpreter','Latex')
text(0.175,0.225,room,'Interpreter','Latex')
text(0.34,0.075,wall,'Interpreter','Latex')
hold off
