close all
clear
load('B1654run35.mat');
%%
%Data Process
index=find(FZ<-110&FZ>-1150);
Fz0=FZ(index)';
Fx0=FX(index)';
SL0=SL(index)';
%Inclination Angle
gamma=(IA(index)*pi/180)';
%scatter3(SL0,Fz0,Fx0);
beta=0.7*ones(36,1);
%%
%Nonlinear Regression MDL
x(1,:)=Fz0';
x(2,:)=SL0';
x(3,:)=gamma';
opts=statset('Display','iter','TolFun',1e-40);
mdl=fitnlm(x',Fx0,@LongitudinalForceFx,beta,'Options',opts);
%%
%Examine
betaFitted=mdl.Coefficients.Estimate;
y=LongitudinalForceFx(betaFitted,x');
%%
scatter3(SL0,Fz0,Fx0,'r');
hold on;
scatter3(x(2,:),x(1,:),y,'b');
xlabel('纵向滑移率');
ylabel('垂直载荷(N)');
zlabel('纵向力(N)');
