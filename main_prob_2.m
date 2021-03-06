clc
clear all
% close all
%%
%=========================Problem 2======================================
l2=12;l3=18;
om2=.5;om3=.75;
delt=.24;
t=0:delt:24;
tit20=deg2rad(30);tit30=deg2rad(15);
tit2=om2*t+tit20;
tit3=om3*t+tit30;

x=l2*cos(tit2)+l3*cos(tit3);
y=l2*sin(tit2)+l3*sin(tit3);
vx=-l2*om2*sin(tit2)-l3*om3*sin(tit3);
vy=l2*om2*cos(tit2)+l3*om3*cos(tit3);
ax=-l2*om2^2*cos(tit2)-l3*om3^2*cos(tit3);
ay=-l2*om2^2*sin(tit2)-l3*om3^2*sin(tit3);

plot(x,y,'m','LineWidth',2)
% figure(1)
% subplot(2,2,3),plot(x,y,'m','LineWidth',2)
% subplot(2,2,1),plot(t,x,'c-.','LineWidth',2)
% subplot(2,2,2),plot(t,y,'g:','LineWidth',2)
% 
% 
% figure(2)
% subplot(2,2,3),plot(vx,vy,'m','LineWidth',2)
% subplot(2,2,1),plot(t,vx,'c-.','LineWidth',2)
% subplot(2,2,2),plot(t,vy,'g:','LineWidth',2)
% 
% 
% figure(3)
% subplot(2,2,3),plot(ax,ay,'m','LineWidth',2)
% subplot(2,2,1),plot(t,ax,'c-.','LineWidth',2)
% subplot(2,2,2),plot(t,ay,'g:','LineWidth',2)


% figure(1)
% subplot(2,2,1)
% title('X-cordinate motion vs time of End Effector')
% xlabel('Time(s)')
% ylabel('X-cordinate motion(cm)')
% subplot(2,2,2)
% title('Y-cordinate motion vs time of End Effector')
% xlabel('Time(s)')
% ylabel('Y-cordinate motion(cm)')
% subplot(2,2,3)
% title('Trajectory of End effector')
% xlabel('X-cordinate motion(cm)')
% ylabel('Y-cordinate motion(cm)')
%  hi1 = suptitle('Position Analysis of End Effector');
% set(hi1,'FontSize',16,'FontWeight','bold')
% figure(2)
% subplot(2,2,1)
% title(' Velocity of X-cordinate vs time of End Effector')
% xlabel('Time(s)')
% ylabel('V_x (cm/s)')
% subplot(2,2,2)
% title('Velocity of Y-cordinate vs time of End Effector')
% xlabel('Time(s)')
% ylabel('V_y (cm/s)')
% subplot(2,2,3)
% title('Trajectory of Velocity components of End effector')
% xlabel('V_x (cm/s)')
% ylabel('V_y (cm/s)')
%  hi2 = suptitle('Velocity Analysis of End Effector');
% set(hi2,'FontSize',16,'FontWeight','bold')
% figure(3)
% subplot(2,2,1)
% title('Accelration of X-cordinate vs time of End Effector')
% xlabel('Time(s)')
% ylabel('A_x (cm^2/s)')
% subplot(2,2,2)
% title(' Acceleration of Y-cordinate vs time of End Effector')
% xlabel('Time(s)')
% ylabel('A_y (cm^2/s)')
% subplot(2,2,3)
% title('Trajectory of Acceleration components of End effector')
% xlabel('A_x (cm^2/s)')
% ylabel('A_y (cm^2/s)')
%  hi3 = suptitle('Accelration Analysis of End Effector');
% set(hi3,'FontSize',16,'FontWeight','bold')
