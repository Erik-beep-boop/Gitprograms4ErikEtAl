clear all;
close all;
clc;
s = settings;
s.matlab.appearance.figure.GraphicsTheme.TemporaryValue = "light";

ordning=7;
m = 51;
rho = 1;
c = 1;
beta = 0;


BC = 4;

C=[1/rho*(c^2) 0
   0 rho];
Cinv = [1 0
        0 1];

A=[0 1
   1 0];

D_mat =[beta 0
   0 0];
alpha=-1; % alpha<=0 for well-posed


e1=[1 0];e2=[0 1]; % Pick out variable
I_2=eye(2);

CFL=0.1; %CFL=k/h
rr=0.1; %Width of Gaussian

t_1=1.8;
x_l=-1;x_r=1;
bredd=x_r-x_l;
y_d=-2.1;y_u= 2.1;


h=bredd/(m-1);
n=2*m;

Val_operator_SC_PDE;  % Change here to make use of upwind SBP

new=2;		  	        % tid (n+1)
old=1;	        	    % tid (n)
temp_tid=0;		        % tempor?r vid tidskiftningen

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
    L=[kron(e1+beta*e2,e_1');
       kron(e1-beta*e2,e_m')];
elseif BC == 4
    L=[kron(e1,e_1');
       kron(e1,e_m')
       kron(e2,e_1');
       kron(e2,e_m')];        
end


P=I_m2-HI2*L'*inv(L*HI2*L')*L; % Projection operator
Dx = kron(A, eye(m)) .* [zeros(m), Dp; Dm zeros(m)];
%disp(Dx);

M = - P * eye(102) * Dx * P;
Meig = eig(M);
disp(Meig);