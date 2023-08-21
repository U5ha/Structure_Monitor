clc
clear
default_plot
%Measeurement File
load("HW2_measurements.txt")
t=HW2_measurements(:,1);
dt=t(2)-t(1);
ydd_1=HW2_measurements(:,2);
ydd_2=HW2_measurements(:,3);
y=[ydd_1 ydd_2]';
m1=10;m2=10;
M=[10 0;0 10];
%Initializing P
P_choice=4;
Q_choice=2;
switch P_choice
    case 1
        P_00=zeros(8);
        P_00(1,1)=10^-2;P_00(2,2)=10^-2;P_00(3,3)=10^-2;P_00(4,4)=10^-2;    
        P_00(5,5)=10^10;P_00(6,6)=10^10;P_00(7,7)=10^10;P_00(8,8)=10^10;
    case 2
        P_00=zeros(8);
        P_00(1,1)=10^-6;P_00(2,2)=10^-6;P_00(3,3)=10^-6;P_00(4,4)=10^-6;    
        P_00(5,5)=10^6;P_00(6,6)=10^6;P_00(7,7)=10^6;P_00(8,8)=10^6;
    case 3
        P_00=zeros(8);
        P_00(1,1)=10^-6;P_00(2,2)=10^-6;P_00(3,3)=10^-6;P_00(4,4)=10^-6;    
        P_00(5,5)=10^4;P_00(6,6)=10^4;P_00(7,7)=10^4;P_00(8,8)=10^4;
    case 4
        P_00=zeros(8);
        P_00(1,1)=10^-10;P_00(2,2)=10^-10;P_00(3,3)=10^-8;P_00(4,4)=10^-8;    
        P_00(5,5)=10^8;P_00(6,6)=10^8;P_00(7,7)=10^4;P_00(8,8)=10^4;
end
%Initializing Q
switch Q_choice
    case 1
        Q=zeros(8);
    case 2
        Q=zeros(8);
        Q(5,5)=50;Q(6,6)=50;Q(7,7)=5e-6;Q(8,8)=5e-6;
    case 3
        Q=zeros(8);
        Q(5,5)=5000;Q(6,6)=5000;Q(7,7)=5e-4;Q(8,8)=5e-4;
end

%El centro file
load("HW1_Q2_El_Centro.txt")
t_el_centro=HW1_Q2_El_Centro(:,1);
ug_el_centro=HW1_Q2_El_Centro(:,2);
%Interpolating El Centro ground Motion
ug=interp1(t_el_centro,ug_el_centro,t);
Ig=[1;1];
uk=(Ig)*ug';

%Modelling and calculating Jacobian
syms y_1 y_2 yd_1 yd_2 c1 c2 k1 k2 ugdd1 ugdd2
f=[yd_1;
   yd_2;
   -ugdd1-(1/m1)*(c1+c2)*yd_1+(c2/m1)*yd_2-(1/m1)*(k1+k2)*y_1+(k2/m1)*y_2;
   -ugdd2+(c2/m2)*yd_1-(c2/m2)*yd_2+(k2/m2)*y_1-(k2/m2)*y_2;
   0;0;0;0];
state=[y_1; y_2; yd_1; yd_2;k1; k2; c1; c2];
H=[-(1/m1)*(c1+c2)*yd_1+(c2/m1)*yd_2-(1/m1)*(k1+k2)*y_1+(k2/m1)*y_2;
    (c2/m2)*yd_1-(c2/m2)*yd_2+(k2/m2)*y_1-(k2/m2)*y_2];
A1=dt*f;
A2=dt*subs(f,state,state+A1*0.75);
F=state+(1/3)*A1+(2/3)*A2;

J=jacobian(F,[y_1 y_2 yd_1 yd_2 k1 k2 c1 c2]);
G=jacobian(H,[y_1 y_2 yd_1 yd_2 k1 k2 c1 c2]);

%Initializing State vector
x_hat_00=[0;0;0;0;5000;5000;50;50];
x_hat(:,1)=x_hat_00;
P{1}=P_00;
I=eye(8);
R=[0.0438 0; 0 0.0967];
for i=1:length(t)-1
    %Prediction zhat k+1,k
    x_hat(:,i+1)=Substitution_F(x_hat(:,i),uk(:,i));
    J1=Substitution_J(x_hat(:,i+1),uk(:,i));
    %Pk+1,k
    P{i+1}=J1*P{i}*J1'+Q;
    %Kalman Gain
    G1=Substitution_G(x_hat(:,i+1),uk(:,i));
    kalman_gain=P{i+1}*G1'*inv(G1*P{i+1}*G1'+R);
    Kk{i}=kalman_gain;
    %Updation zhat k+1,k+1
    H1=Substitution_H(x_hat(:,i+1),uk(:,i+1));
    x_hat(:,i+1)=x_hat(:,i+1)+kalman_gain*(y(:,i+1)-H1);
    %Pk+1,k+1
    P{i+1}=(I-kalman_gain*G1)*P{i+1}*(I-kalman_gain*G1)'+kalman_gain*R*kalman_gain';
end
head=["Displacement 1st Floor","Displacement 2nd Floor","Velocity 1st Floor","Velocity 2nd Floor","Stiffness of 1st Floor K_1","Stiffness of 2nd Floor K_2","Damping 1st Floor C_1","Damping 2nd Floor C_2"];
y_def=["$y_1 (m)$","$y_2 (m)$","$\dot y_1 (m/s)$","$\dot y_2 (m/s)$","$K_1 (N/m)$","$K_2 (N/m)$","$C_1 (Ns/m)$","$C_2 (Ns/m)$"];

for i=1:8
    H1=gca;
    H1.LineWidth=2;
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperSize', [16 9]);
    % End figure setup
    set(gca, 'Fontsize', 25)
    figure(i)
    hold on
    box on
    plot(t,x_hat(i,:))
    xlabel('$time (s)$','Interpreter','latex')
    ylabel(y_def(i),'Interpreter','latex')
    title(head(i),'FontSize',16)
end
%%
for i=1:length(t)
    if t(i)<=10
        K1_original(i)=10000;
    else
        K1_original(i)=7500;
    end
    K2_original(i)=10000;
    C1_original(i)=31.76;
    C2_original(i)=31.76;
end
figure(5)
hold on
plot(t,K1_original)
legend('EKF Estimated Parameter','Original')
figure(6)
hold on
ylim([0 12000])
plot(t,K2_original)
legend('EKF Estimated Parameter','Original')
figure(7)
hold on
plot(t,C1_original)
legend('EKF Estimated Parameter','Original')
figure(8)
H1=gca;
H1.LineWidth=2;
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [16 9]);
% End figure setup
set(gca, 'Fontsize', 25)
hold on
plot(t,C2_original)
legend('EKF Estimated Parameter','Original')