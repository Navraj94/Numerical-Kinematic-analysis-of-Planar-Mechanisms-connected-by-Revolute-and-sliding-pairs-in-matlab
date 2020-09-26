clc
clear all
close all

%%
%========================First Question-slider Crank=======================

l2=3;l3=9;h=1;
tol=1e-6;maxit=5000;

%% ====================== NEWTON RAPHSON-==================================
% initialize all values above and just copy paste this section where ever
% required have written with generality
% Dont forget
% 1.system equations function (named as per requirement just change tat lazy)
% 2.jacobain(x) function (named as per requirement just change tat lazy)
% func=@(param,x)slicrank(param,x);
% jacob=@(param,x)sli_jaco(param,x);
tita=0:.1:2*pi+.1;
param=zeros(4,1);
param=[l2;0;l3;h];
for i=1:length(tita)
    param(2)=tita(i);
    x=[];
    x(:,1)=[1;10];
    f=slicrank(param,x);
    jac=sli_jaco(param,x);
    conv=2;
    j=2;
    while conv > tol && j < maxit
        delt=jac\f;
        x(:,j)=x(:,j-1)-delt;
        f=slicrank(param,x(:,j));
        jac=sli_jaco(param,x(:,j));
        conv=norm(delt);
        j=j+1;
    end
    xc(:,i)=x(:,j-1);
    vc(:,i)=inv(jac)*150*[l2*sin(tita(i));-l2*cos(tita(i))];
    rh=-150^2*[l2*cos(tita(i));l2*sin(tita(i))]+vc(1,i)^2*[l3*cos(xc(1,i));l3*sin(xc(1,i))];
    ac(:,i)=inv(jac)*rh;
    
    
end

    figure(1)
    subplot(3,2,1),plot(tita,xc(1,:),'--','LineWidth',2)
    xticks([0 pi/2 pi 3*pi/2 2*pi])
    xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
    title('\theta3 vs \theta2(varied from 0 to 2*\pi)')
    hold on
    grid on
    subplot(3,2,2),plot(tita,xc(2,:),'-.','LineWidth',2)
    xticks([0 pi/2 pi 3*pi/2 2*pi])
    xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
    title('x4 vs \theta2(varied from 0 to 2*\pi)')
    hold on
    grid on
    subplot(3,2,3),plot(tita,vc(1,:),'--','LineWidth',2)
    xticks([0 pi/2 pi 3*pi/2 2*pi])
    xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
    title(' Velocity of \theta3 vs \theta2(varied from 0 to 2*\pi)')
    hold on
    grid on
    subplot(3,2,4),plot(tita,vc(2,:),'-.','LineWidth',2)
    xticks([0 pi/2 pi 3*pi/2 2*pi])
    xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
    title(' Velocity of x4 vs \theta2(varied from 0 to 2*\pi)')
    hold on
    grid on
    subplot(3,2,5),plot(tita,ac(1,:),'--','LineWidth',2)
    xticks([0 pi/2 pi 3*pi/2 2*pi])
    xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
    title(' Acceleration of \theta3 vs \theta2(varied from 0 to 2*\pi)')
    hold on
    grid on
    subplot(3,2,6),plot(tita,ac(2,:),'-.','LineWidth',2)
    xticks([0 pi/2 pi 3*pi/2 2*pi])
    xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
    title(' Acceleration of x4 vs \theta2(varied from 0 to 2*\pi)')
    hold on
    grid on

 hi = suptitle('Solution by Newton-Raphson Iteration Method')
set(hi,'FontSize',16,'FontWeight','bold')


% xc=x(:,j-1);
% jac=sli_jaco(param,xc);
% f=slicrank(param,xc);
% vel=jac*xc*150;
% acc=jac*vel;

 %%
% %=========================Problem 2======================================
% l2=12;l3=18;
% om2=.5;om3=.75;
% delt=.24;
% t=0:delt:24;
% tit20=deg2rad(30);tit30=deg2rad(15);
% tit1=om2*t+tit20;
% tit2=om3*t+tit30;
% 
% x=l2*cos(tit1)+l3*cos(tit2);
% y=l2*sin(tit1)+l3*sin(tit2);
% figure(40)
% plot(x,y,'m','LineWidth',2)


