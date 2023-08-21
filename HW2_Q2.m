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
   0;
   0;
   0;
   0];
state=[y_1; y_2; yd_1; yd_2;k1; k2; c1; c2];
H=[-(1/m1)*(c1+c2)*yd_1+(c2/m1)*yd_2-(1/m1)*(k1+k2)*y_1+(k2/m1)*y_2;
    (c2/m2)*yd_1-(c2/m2)*yd_2+(k2/m2)*y_1-(k2/m2)*y_2];
A1=f;
A2=subs(f,state,state+A1*0.75*dt);
F=state+(1/3)*dt*A1+(2/3)*dt*A2;
%Initializing State vector
x_hat_00=[0;0;0;0;5000;5000;50;50];
x_hat(:,1)=x_hat_00;
P{1}=P_00;
I=eye(8);
R=[0.0438 0;
   0 0.0967];
%UKF parameters
N=8;lambda=-7.97;alpha=0.1;beta=2;
Wm0=lambda/(N+lambda);
Wc0=lambda/(N+lambda)+(1-alpha^2+beta);
Wmi=1/(2*(N+lambda));
Wci=1/(2*(N+lambda));
Wm=ones(1,2*N+1);Wc=ones(1,2*N+1);
Wm(1)=Wm0;Wc(1)=Wc0;
Wm(2:end)=Wmi;Wc(2:end)=Wci;

for i=1:length(t)-1
    %Selecting sigma points
    temp1=(N+lambda).*P{i};
    [root1]=chol(temp1)';
    sigma1=x_hat(:,i).*ones(N,2*N+1);
    sigma1(:,2:N+1)=sigma1(:,2:N+1)+root1;
    sigma1(:,N+2:2*N+1)=sigma1(:,N+2:2*N+1)-root1;
    %Passing sigma points through the function F
    for e1=1:2*N+1
        [Xi(:,e1)]=Substitution_F(sigma1(:,e1),uk(:,i));
    end
    %Predicting Z_hat
    sum1=zeros(8,1);
    for j=1:2*N+1
        sum1=sum1+Wm(j).*Xi(:,j);
    end
    x_hat(:,i+1)=sum1;
    %Predicting P
    sum2=0;diff=0;
    for m=1:2*N+1
        diff=Xi(:,m)-x_hat(:,i+1);
        sum2=sum2+Wc(m).*(diff*diff');
    end
    P{i+1}=sum2+Q;
    %Updation
    %Redrawing sigma points from predicted P matrix and state vector to pass through function H 
    temp2=(N+lambda).*P{i+1};
    [root2]=chol(temp2)';
    sigma2=x_hat(:,i+1).*ones(N,2*N+1);
    sigma2(:,2:N+1)=sigma2(:,2:N+1)+root2;
    sigma2(:,N+2:2*N+1)=sigma2(:,N+2:2*N+1)-root2;
    %Passing sigma points through the function F
    for e3=1:2*N+1
        [Xi_new(:,e3)]=Substitution_F(sigma2(:,e3),uk(:,i));
    end
    
    for e2=1:2*N+1
        [yi(:,e2)]=Substitution_H(Xi_new(:,e2),uk(:,i+1));
    end
    sum3=0;
    for n1=1:2*N+1
        sum3=sum3+Wm(n1).*yi(:,n1);
    end
    y_hat=sum3;
    %Calculating the 2 covarinace Pyy and Pxy
    sum4=0;sum5=0;diff1=zeros(2,1);diff2=zeros(8,1);
    for a1=1:2*N+1
        diff1=yi(:,a1)-y_hat;
        sum4=sum4+Wc(a1).*(diff1*diff1.');
       % a1
        diff2=Xi(:,a1)-x_hat(:,i+1);
        sum5=sum5+Wc(a1).*(diff2*diff1');
    end
    sum4=sum4+R;
    Pyy=sum4;
    Pxy=sum5;
    Kk=Pxy/(Pyy);
    %Updating state
    x_hat(:,i+1)=x_hat(:,i+1)+Kk*(y(:,i+1)-y_hat);
    %Updating P matrix
    P{i+1}=P{i+1}-Kk*Pyy*Kk';
end
%%
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
legend('UKF Estimated Parameter','Original')
figure(6)
hold on
plot(t,K2_original)
legend('UKF Estimated Parameter','Original')
figure(7)
hold on
plot(t,C1_original)
legend('UKF Estimated Parameter','Original')
figure(8)
H1=gca;
H1.LineWidth=2;
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [16 9]);
% End figure setup
set(gca, 'Fontsize', 25)
hold on
plot(t,C2_original)
legend('UKF Estimated Parameter','Original')