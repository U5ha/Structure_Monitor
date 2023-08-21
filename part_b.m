clc
clear
dt=0.01;
A=[1 dt 0 0;
    0 1 0 0;
    0 0 1 dt;
    0 0 0 1];
B=[0;0;dt^2/2;dt];
C=[0 1 0 0;
    0 0 0 1];
D=[0;0];
sysd=ss(A,B,C,D,dt);
load("HW1_Q1_Velocity_Measurements.txt")
t=HW1_Q1_Velocity_Measurements(:,1);
vx=HW1_Q1_Velocity_Measurements(:,2);
vy=HW1_Q1_Velocity_Measurements(:,3);
y=[vx vy]';
R=[4 0;
    0 16];
g=-9.81;
u=g*ones(1,length(t));
%Initializing
P_00=10^4*eye(4);
z_hat_00=[300 100 200 60]';
%z_hat=zeros(4,length(t));
z_hat(:,1)=z_hat_00;
P{1}=P_00;
Q=zeros(4);
I=eye(4);
A=sysd.A;B=sysd.B;C=sysd.C;D=sysd.D;

for i=1:length(t)-1
    %Prediction zhat k+1,k
    z_hat(:,i+1)=A*z_hat(:,i)+B*u(i);
    %Pk+1,k
    P{i+1}=A*P{i}*A'+Q;
    %Kalman Gain
    kalman_gain=P{i+1}*C'*inv(C*P{i+1}*C'+R);
    Kk{i}=kalman_gain;
    %Updation zhat k+1,k+1
    z_hat(:,i+1)=z_hat(:,i+1)+kalman_gain*(y(:,i+1)-C*z_hat(:,i+1)-D*u(:,i+1));
    %Pk+1,k+1
    P{i+1}=(I-kalman_gain*C)*P{i+1}*(I-kalman_gain*C)'+kalman_gain*R*kalman_gain';
end
load("HW1_Q1_Displacement_Measurements.txt")
ux=HW1_Q1_Displacement_Measurements(:,2);
uy=HW1_Q1_Displacement_Measurements(:,3);
%%
figure(1)
plot(t,z_hat(1,:),'b','LineWidth',2)
hold on
plot(t,ux,'--')
xlabel('time (s)','FontSize',16)
ylabel('u_x m','FontSize',16)
box on
grid on
title('Horizontal Displacement vs time','FontSize',16)
ldg=legend('Using Kalman Filter','Measured');
ldg.FontSize=13;
figure(2)
plot(t,z_hat(2,:),'b','LineWidth',2)
hold on
plot(t,vx,'--')
xlabel('time (s)','FontSize',16)
ylabel('v_x m/s','FontSize',16)
box on
grid on
title('Horizontal Velocity vs time','FontSize',16)
ldg=legend('Using Kalman Filter','Measured');
ldg.FontSize=13;
figure(3)
plot(t,z_hat(3,:),'b','LineWidth',2)
hold on
plot(t,uy,'--')
xlabel('time (s)','FontSize',16)
ylabel('u_y m','FontSize',16)
box on
grid on
title('Vertical Displacement vs time','FontSize',16)
ldg=legend('Using Kalman Filter','Measured');
ldg.FontSize=13;
figure(4)
plot(t,z_hat(4,:),'b','LineWidth',2)
hold on
plot(t,vy,'--')
xlabel('time (s)','FontSize',16)
ylabel('v_y m/s','FontSize',16)
box on
grid on
title('Vertical Velocity vs time','FontSize',16)
ldg=legend('Using Kalman Filter','Measured');
ldg.FontSize=13;