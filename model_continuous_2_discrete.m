clc
clear
dt=0.005;
M=[10 0;0 10];
k1=10000;k2=10000;
c1=31.76;c2=31.76;
K=[k1+k2 -k2 ;-k2 k2];
C=[c1+c2 -c2;-c2 c2];
z=zeros(2);
I=eye(2);
Ac=[z I;
    -M\K -M\C];
Bc=[z;inv(M)];
%%
%Acceleration measurements
Cc_a=[-M\K -M\C];
%Dc_a=inv(M);
Dc_a=z;
%%
%Displacement measurements
Cc_d=[I z];
Dc_d=z;
%%
%Discrete system for accleration measurements
sysc_a=ss(Ac,Bc,Cc_a,Dc_a);
sysd_a=c2d(sysc_a,dt,"zoh");
%%
%Discrete system for displacement measurements
sysc_d=ss(Ac,Bc,Cc_d,Dc_d);
sysd_d=c2d(sysc_d,dt,"zoh");
%%
%manual conversion of A and B matrix to discrete system as C and D remains
%the same
[psi,lambda]=eig(Ac);
for i=1:4
    s(i)=lambda(i,i);
end
exponential_matrix=zeros(4);
for i=1:4
    exponential_matrix(i,i)=exp(s(i)*dt);
end
Ad=real(psi*exponential_matrix*inv(psi));
I1=eye(4);
Bd=inv(Ac)*(Ad-I1)*Bc;