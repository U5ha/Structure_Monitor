clc
clear
dt=0.005;
M=[10 0;0 10];
k1=10000;k2=10000;
c1=31.76;c2=31.76;
K=[k1+k2 -k2 ;-k2 k2];
C=[c1+c2 -c2;-c2 c2];
z=zeros(2);
I1=eye(2);
Ac=[z I1;
    -M\K -M\C];
Bc=[z;inv(M)];
%Displacement measurements
Cc=[I1 z];
Dc=z;
%Discrete system for displacement measurements
sysc=ss(Ac,Bc,Cc,Dc);
sysd=c2d(sysc,dt,"zoh");
A=sysd.A;B=sysd.B;C=sysd.C;D=sysd.D;
load("HW1_Q2_El_Centro.txt")
load("HW1_Q2_Measurements.txt")
t_el_centro=HW1_Q2_El_Centro(:,1);
ug_el_centro=HW1_Q2_El_Centro(:,2);
t=HW1_Q2_Measurements(:,1);
y1=HW1_Q2_Measurements(:,2);
y2=HW1_Q2_Measurements(:,3);
yd1=HW1_Q2_Measurements(:,4);
yd2=HW1_Q2_Measurements(:,5);
%Interpolating El Centro ground Motion
ug=interp1(t_el_centro,ug_el_centro,t);
Ig=[1;1];
uk=-(M*Ig)*ug';
%Initializing
z_hat_00=[0 0 0 0]';
z_hat(:,1)=z_hat_00;
R=[2.642e-4 0;0 4.251e-4];
y=[y1 y2]';
P_00=10^2*eye(4);
I=eye(4);
P{1}=P_00;
Q=eye(4)*0.0;
%%
for i=1:length(t)-1
    %Prediction zhat k+1,k
    z_hat(:,i+1)=A*z_hat(:,i)+B*uk(:,i);
    %Pk+1,k
    P{i+1}=A*P{i}*A'+Q;
    %Kalman Gain
    kalman_gain=P{i+1}*C'*inv(C*P{i+1}*C'+R);
    Kk{i}=kalman_gain;
    %Updation zhat k+1,k+1
    z_hat(:,i+1)=z_hat(:,i+1)+kalman_gain*(y(:,i+1)-C*z_hat(:,i+1)-D*uk(:,i+1));
    %Pk+1,k+1
    P{i+1}=(I-kalman_gain*C)*P{i+1}*(I-kalman_gain*C)'+kalman_gain*R*kalman_gain';
end
%%
figure(1)
plot(t,z_hat(1,:),"LineWidth",2)
hold on
plot(t,y1,'--')
xlabel('time (s)','FontSize',16)
ylabel('y_1 m','FontSize',16)
box on
grid on
title('y_1 vs time','FontSize',16)
ldg=legend('Using Kalman Filter','Measured');
ldg.FontSize=13;
figure(2)
plot(t,z_hat(2,:),"LineWidth",2)
hold on
plot(t,y2,'--')
xlabel('time (s)','FontSize',16)
ylabel('y_2 m','FontSize',16)
box on
grid on
title('y_2 vs time','FontSize',16)
ldg=legend('Using Kalman Filter','Measured');
ldg.FontSize=13;
figure(3)
plot(t,z_hat(3,:),"LineWidth",1)
xlabel('time (s)','FontSize',16)
ylabel('$\dot y_1$ m/s','Interpreter','latex','FontSize',16)
box on
grid on
title('$\dot y_1$ vs time','Interpreter','latex','FontSize',16)
ldg=legend('Using Kalman Filter');
ldg.FontSize=13;
figure(4)
plot(t,z_hat(3,:),"LineWidth",1)
xlabel('time (s)','FontSize',16)
ylabel('$\dot y_2$ m/s','Interpreter','latex','FontSize',16)
box on
grid on
title('$\dot y_2$ vs time','Interpreter','latex','FontSize',16)
ldg=legend('Using Kalman Filter');
ldg.FontSize=13;