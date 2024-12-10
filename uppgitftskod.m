clear all;
close all;
clc;
s = settings;
s.matlab.appearance.figure.GraphicsTheme.TemporaryValue = "light";

order=0;
m = 201;


scrsz = get(0,'ScreenSize');
vidObj = VideoWriter('Wave1.mp4','MPEG-4');
open(vidObj);

I = eye(m,m);

e_1=zeros(m,1);e_1(1)=1;
e_m=zeros(m,1);e_m(m)=1;

t_1=25;           % End time
t=0;

h=domain/m;      % Grid-step
k=h/10;

max_itter=floor(t_1/k);

% How often to update the movie
n_step=100;