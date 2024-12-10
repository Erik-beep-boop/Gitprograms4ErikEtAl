clear all;
close all;
clc;
gridsize = 79;
m = 2*gridsize;
h = 1;
k = 0.5;
T = 0;

I = eye(m,m);

e_1=zeros(m,1);e_1(1)=1;
e_m=zeros(m,1);e_m(m)=1;




H=diag(ones(gridsize,1),0);
H(1:6,1:6)=[0.19087e5 / 0.60480e5 0 0 0 0 0; 0 0.84199e5 / 0.60480e5 0 0 0 0; 0 0 0.18869e5 / 0.30240e5 0 0 0; 0 0 0 0.37621e5 / 0.30240e5 0 0; 0 0 0 0 0.55031e5 / 0.60480e5 0; 0 0 0 0 0 0.61343e5 / 0.60480e5;];
H(gridsize-5:gridsize,gridsize-5:gridsize)=fliplr(flipud(H(1:6,1:6)));
H=H*h;
HI=inv(H);

Qp=(-1/105*diag(ones(gridsize-3,1),-3)+1/10*diag(ones(gridsize-2,1),-2)-3/5*diag(ones(gridsize-1,1),-1)-1/4*diag(ones(gridsize,1),0)+1*diag(ones(gridsize-1,1),+1)-3/10*diag(ones(gridsize-2,1),+2)+1/15*diag(ones(gridsize-3,1),+3)-1/140*diag(ones(gridsize-4,1),+4));

Q_U =[-0.265e3 / 0.300272e6 0.1587945773e10 / 0.2432203200e10 -0.1926361e7 / 0.25737600e8 -0.84398989e8 / 0.810734400e9 0.48781961e8 / 0.4864406400e10 0.3429119e7 / 0.202683600e9; -0.1570125773e10 / 0.2432203200e10 -0.26517e5 / 0.1501360e7 0.240029831e9 / 0.486440640e9 0.202934303e9 / 0.972881280e9 0.118207e6 / 0.13512240e8 -0.231357719e9 / 0.4864406400e10; 0.1626361e7 / 0.25737600e8 -0.206937767e9 / 0.486440640e9 -0.61067e5 / 0.750680e6 0.49602727e8 / 0.81073440e8 -0.43783933e8 / 0.194576256e9 0.51815011e8 / 0.810734400e9; 0.91418989e8 / 0.810734400e9 -0.53314099e8 / 0.194576256e9 -0.33094279e8 / 0.81073440e8 -0.18269e5 / 0.107240e6 0.440626231e9 / 0.486440640e9 -0.365711063e9 / 0.1621468800e10; -0.62551961e8 / 0.4864406400e10 0.799e3 / 0.35280e5 0.82588241e8 / 0.972881280e9 -0.279245719e9 / 0.486440640e9 -0.346583e6 / 0.1501360e7 0.2312302333e10 / 0.2432203200e10; -0.3375119e7 / 0.202683600e9 0.202087559e9 / 0.4864406400e10 -0.11297731e8 / 0.810734400e9 0.61008503e8 / 0.1621468800e10 -0.1360092253e10 / 0.2432203200e10 -0.10677e5 / 0.42896e5;];
Qp(1:6,1:6)=Q_U;
Qp(gridsize-5:gridsize,gridsize-5:gridsize)=flipud( fliplr(Q_U(1:6,1:6) ) )'; %%% This is different fromstandard SBP

Qm=-Qp';

e_1_gridsize=zeros(gridsize,1);e_1_gridsize(1)=1;
e_gridsize=zeros(gridsize,1);e_gridsize(gridsize)=1;

Dp=HI*(Qp-1/2*e_1_gridsize*e_1_gridsize'+1/2*e_gridsize*e_gridsize') ;
Dm=HI*(Qm-1/2*e_1_gridsize*e_1_gridsize'+1/2*e_gridsize*e_gridsize') ;



D_X = zeros(m,m);
D_X(gridsize+1:2*gridsize, 1:gridsize) = Dp;
D_X(1:gridsize, gridsize+1:2*gridsize) = Dm;










%ORDER 7 BIG H
H=diag(ones(m,1),0);
H(1:6,1:6)=[0.19087e5 / 0.60480e5 0 0 0 0 0; 0 0.84199e5 / 0.60480e5 0 0 0 0; 0 0 0.18869e5 / 0.30240e5 0 0 0; 0 0 0 0.37621e5 / 0.30240e5 0 0; 0 0 0 0 0.55031e5 / 0.60480e5 0; 0 0 0 0 0 0.61343e5 / 0.60480e5;];
H(m-5:m,m-5:m)=fliplr(flipud(H(1:6,1:6)));
H=H*h;
HI=inv(H);

L = zeros(2, m);
%L(1,1) = 1;
L(1,gridsize+1) = 1;
%L(2,gridsize) = 1;
L(2,2*gridsize) = 1; %-1
LT = transpose(L);

C = eye(m, m); %BECAUSE WE ASSUME p = 1, c = 1
%C(1:floor(gridsize/2), 1:floor(gridsize/2)) = 5*eye(floor(gridsize/2),floor(gridsize/2));

%C(1:20, 1:20) = 2*eye(20,20);
%C(21, 21) = 1.5;
%C(22, 22) = 1.1;
%C(23, 23) = 0.8;
%C(24, 24) = 0.55;
%C(25, 25) = 0.35;


%C(5,5) = 0.5;
%C(6,6) = 1;
%C(7,7) = 1.5;
%C(8,8) = 2.5;
%C(9,9) = 3.5;
%C(10,10) = 4;
%C(11,11) = 3.5;
%C(12,12) = 2.5;
%C(13,13) = 1.5;
%C(14,14) = 1;
%C(15,15) = 0.5;


PART = L*HI*LT;
SANDWIDSH = LT * inv(PART) * L;
P = eye(m, m) - HI * SANDWIDSH;
% image((P+1)*255/2)




u = zeros(m, 1);
for i = 0:20
    u(i+30) = (1-cos(2*i/6.366));
end
for i = 0:gridsize-1
    u(i+1) = (cos(2*pi*i/(gridsize-1)));
end
%u(30:50) = 1;

MATRIX = -P*(D_X/C)*P;
% image(MATRIX*100 + 255/2)
pause(1.5)

while T < 100
    %RUNGE KUTTA FYRA
    k1 = k*(MATRIX*u);
    k2 = k*(MATRIX*(u+k1/2));
    k3 = k*(MATRIX*(u+k2/2));
    k4 = k*(MATRIX*(u+k3));
    u_t = (1/6)*(k1 + 2*k2 + 2*k3 + k4);
    
    uTHEORETICAL = zeros(m, 1);

    u = u + k*u_t;
    T = T + k;
    plot(u(gridsize+1:m))
    %plot(1.1*uTHEORETICAL(1:gridsize))
    hold on
    plot(u(1:gridsize))
    axis([0 gridsize -2 2])
    title(['Waves at t = ',num2str(T)])
    legend('Pressure','Velocity')
    hold off
    pause(0.005)
      
end


