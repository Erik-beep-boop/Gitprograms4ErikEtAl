clear, clc, close all
s = settings;
s.matlab.appearance.figure.GraphicsTheme.TemporaryValue = "light";

ordning=6;
m = 201;
rho = 1;
c = 1;
beta = 0;
BC = 1;
% How often to update the movie
n_step=10;

scrsz = get(0,'ScreenSize');
clc;

vidObj = VideoWriter('Wave.mp4','MPEG-4');
open(vidObj);

C=[1/rho*(c^2) 0
   0 rho];

A=[0 1
   1 0];

D_mat =[beta 0
   0 0];


e1=[1 0];e2=[0 1]; % Pick out variable
I_2=eye(2);

CFL=0.05; %CFL=k/h
rr=0.1; %Width of Gaussian

t_1=1.8;
x_l=-1;x_r=1;bredd=x_r-x_l;
y_d=-2.1;y_u= 2.1;


h=bredd/(m-1);
n=2*m;

Val_operator_SC_PDE;  % Change here to make use of upwind SBP

new=2;		  	        % tid (n+1)
old=1;	        	    % tid (n)
temp_tid=0;		        

t=0;

dt=CFL*h;               % CFL

max_itter=floor(t_1/dt);

e_1=zeros(m,1);e_1(1)=1;
e_m=zeros(m,1);e_m(m)=1;

I_m=eye(m);
HI2=kron(I_2,HI);
I_m2=kron(I_2,I_m);


if BC == 1
    L=[kron(e1,e_1');
       kron(e1,e_m')];
elseif BC == 2
    L=[kron(e2,e_1');
       kron(e2,e_m')];
elseif BC == 3
    L=[kron(e1-beta*e2,e_1');
       kron(e1+beta*e2,e_m')];
elseif BC == 4
    L=[kron(e1,e_1');
       kron(e1,e_m')
       kron(e2,e_1');
       kron(e2,e_m')];        
end


P=I_m2-HI2*L'*inv(L*HI2*L')*L; % Projection operator

if (ordning == 5 | ordning == 7 | ordning == 9)

    Dpm = zeros(2*m,2*m);
    Dpm(1:m, m+1:2*m) = Dp;
    Dpm(m+1:2*m, 1:m) = Dm;
    T=P*Dpm*P;
    %size(Dpm) % byt till minus T
else
    T=P*kron(A,D1)*P; 
end


D=sparse(kron(inv(C),I_m)*T);
   
temp=zeros(n,1);  % temporary solutionvector in RK4 w1,..,w4

w1=zeros(n,1);      % First step in RK
w2=zeros(n,1);      % Second step in RK
w3=zeros(n,1);      % Third step in RK
w4=zeros(n,1);      % Fourt step in RK

x=linspace(x_l,x_r,m);	  % x points
    
felet=zeros(max_itter+1,1);	% Error vector

%if BC==1
uc1=exp(-((x+t)/rr).^2)';
uc2=-exp(-((x-t)/rr).^2)';


%figure(1);
%plot(x,uc1,'r',x,uc2,'b');

exakt=zeros(n,1);	        % Analytiska lösningen
exakt(1:m)=uc1-uc2;
exakt(m+1:n)=uc2+uc1;

V=zeros(n,2);	                % Lösningsvektorn
V(1:m,old)=uc1-uc2; % pressure
V(m+1:n,old)=uc2+uc1; % velocity
felet(1)=sqrt(h)*norm(V(:,old)-exakt);

     
Plotta_Maxwell_1D; % Plot the solution and generate a movie



for nr_itter=1:max_itter

%%%--------Start R-K utan dissipation----------------------------------------
  w1=D*V(:,old);
  temp=V(:,old)+dt/2*w1;
      
  w2=D*temp;
  temp=V(:,old)+dt/2*w2;
    
  w3=D*temp;
  temp=V(:,old)+dt*w3;
  
  w4=D*temp;
    
  V(:,new)=V(:,old)+dt/6*(w1+2*w2+2*w3+w4);
  
  t=t+dt; 

  temp_tid=old;old=new;new=temp_tid;
  
  if mod(nr_itter,n_step)==0
      Plotta_Maxwell_1D; % Plot the solution and generate a movie
  end
  
exakt=zeros(n,1);	        % Analytiska l?sningen
exakt(1:m)=uc2-uc1;
exakt(m+1:n)=uc2+uc1;

% Analytic solution after t=1.8 with BC type 1
tt=2.0-t;
uc1=exp(-((x-tt)/rr).^2)';
uc2=-exp(-((x+tt)/rr).^2)';
exakt(1:m)=uc1-uc2;
exakt(m+1:n)=uc1+uc2;
felet(nr_itter+1)=sqrt(h)*norm(-V(:,old)-exakt);
  
end

  %differens(i)=felet(max_itter);
  
 close(vidObj);
 tiden=linspace(0,t_1,max_itter+1);
 
 if BC==1
     disp(['The l_2 error is given by : ', num2str(felet(nr_itter+1))])
     % decimal = num2str(felet(nr_itter+1), '%.20f');
     % disp(['The l_2 error is given by : ', decimal ]);   
 end
 
figure(2);
plot(x,V(1:m,old),'r',x,-V(m+1:n,old),'b--','LineWidth',1);
title(['Numerical solution at t = ',num2str(t_1)]);
axis([x_l x_r y_d y_u]);
grid;xlabel('x');
legend('Pressure','Velocity')
ax = gca;          % current axes
ax.FontSize = 16;


% plot starting values
t = 0;
uc1=exp(-((x-t)/rr).^2)';
uc2=-exp(-((x+t)/rr).^2)';
V(1:m,old)=uc1-uc2;
V(m+1:n,old)=uc1+uc2;

figure(3);
plot(x,V(1:m,old),'r',x,-V(m+1:n,old),'b--','LineWidth',1);
title(['Numerical solution at t = ',num2str(t)]);
axis([x_l x_r y_d y_u]);
grid;xlabel('x');
legend('Pressure','Velocity')
ax = gca;          % current axes
ax.FontSize = 16;




