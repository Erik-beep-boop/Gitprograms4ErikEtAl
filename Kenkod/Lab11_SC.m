%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                               %%%
%%% Lab 1, part 1                 %%%
%%%                               %%%
%%% Scientific Computing for      %%%
%%% PDE                           %%%
%%%                               %%%
%%%                               %%%
%%% Author: Ken Mattsson          %%%
%%% Date:   2022-06-30            %%%
%%%                               %%%
%%% Solve u_t+u_x=0               %%%
%%%                               %%%
%%% Use periodic operators        %%%
%%% Both implicit (spectral)      %%%
%%% and explicit of order         %%%
%%% 2,4,6,8,10 and 12             %%%
%%%                               %%%
%%%                               %%%
%%% Use RK4 to time-integrate     %%%
%%% Test also  convergence        %%%
%%%                               %%%
%%% Here use a smoothy            %%%
%%% initial data                  %%%
%%%                               %%%
%%% The code generates an avi-    %%%
%%% file for numerical solution   %%%
%%%                               %%%
%%% We compare dispersion error   %%%
%%% for explicit and implicit     %%%
%%% finite difference stencils    %%%
%%%                               %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all;
close all;
clc;
% s = settings;
% s.matlab.appearance.figure.GraphicsTheme.TemporaryValue = "light";
disp(' ');   
disp('---------------------------------------------')
disp('Code that solve the first order wave equation in 1-D with periodic boundary condition.')
disp('Here using explicit stencils of orders 2,4,6,8,10,12')
disp('We here test for smooth initial data to compare the dispersion errors (difference between the exact and numerical solution).')
disp(' '); 
disp('The red line is the analytic solution and the blue the numerical solution. An animation of the numerical solution is generated and stored as a mp4-file.')
%disp('The time-step is set to dt=1/80 and the end time of the simulation is t_1=4.')


 order=0;
 while (order ~=2 &&  order ~=4 && order ~=6 && order ~=8 && order ~=10 && order ~=12)
 order = input('Order of accuracy explicit (2) (4) (6) (8) (10) (12): ');
 end
 
 m=0;
 while ( m <=15 )
 m = input('How many grid-intervals (>15) : ');
 end
 
 

scrsz = get(0,'ScreenSize');
%figure('Position',[scrsz(3)/2 scrsz(4) scrsz(3)/2 scrsz(4)])
vidObj = VideoWriter('Wave1.mp4','MPEG-4');
open(vidObj);

x_l=0;x_r=1;    % The boundaries of the domain
domain=x_r-x_l;

%m=401;           % Number of grid-intervals, first grid
t_1=25;           % End time
t=0;

h=domain/m;      % Grid-step
k=h/10;

max_itter=floor(t_1/k);
%max_itter=2000;

% How often to update the movie
n_step=100;

theAxes=[x_l x_r -0.3 1.3]; % Regarding the figure
%theAxes=[x_l x_r -1 1]; % Regarding the figure



[H,Q,M] = Periodic_E_SC ( domain,m, order );
D=Q/h;

temp=zeros(m,1);    % Temporary vector in RK4

w1=zeros(m,1);      % Step 1 vector in RK4
w2=zeros(m,1);      % Step 2 vector in RK4
w3=zeros(m,1);      % Step 3 vector in RK4
w4=zeros(m,1);      % Step 4 vector in RK4

x=linspace(x_l,x_r-h,m);	  % Discrete x-values   

felet=zeros(max_itter+1,1);	      % Error vector in time

exact=Initial_P1(t,x,0);    %Exact solution smooth
  
V=exact;                          % Numerical solution
felet(1)=sqrt(h)*norm(V-exact);

tt=rem(t,domain);

for nr_itter=1:max_itter

%%%--------Start R-K utan dissipation----------------------------------------
  w1=RHS_P(D,V);
  
  w2=RHS_P(D,V+k/2*w1);
   
  w3=RHS_P(D,V+k/2*w2);
 
  w4=RHS_P(D,V+k*w3);
    
  V=V+k/6*(w1+2*w2+2*w3+w4);
  
  t=t+k; 
  tt=rem(t,domain);
  
  
  exact=Initial_P1(tt,x,domain)+Initial_P1(tt,x,0)+Initial_P1(tt,x,-domain);
  
  felet(nr_itter+1)=sqrt(h)*norm(V-exact);
  
  if mod(nr_itter,n_step)==0
      Plotta_1D_SC; % Plott the solution and generate a movie
  end
  
end

  %felet(nr_itter+1);

  disp(['The l_2-error is: ', num2str(felet(nr_itter+1))])
 
 tiden=linspace(0,t_1,max_itter+1);
 close(vidObj);
 
 
    figure(2);
    plot(x,exact,'r',x,V,'b--','LineWidth',1);
    xlabel('x');
    ylabel('u');
    title(['Numerical and exact smooth solution at t = ',num2str(t), ' Order = ', num2str(order)])
    legend('Exact solution','Numerical solution')
    %title(['Wave propagation',ordningstyp])
    axis(theAxes);
    ax = gca; % current axes
    ax.FontSize = 16;
    
    
