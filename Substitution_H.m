function [H_value] = Substitution_H(z_hat,force)
y_1=z_hat(1,1);
y_2=z_hat(2,1);
yd_1=z_hat(3,1);
yd_2=z_hat(4,1);
k1=z_hat(5,1);
k2=z_hat(6,1);
c1=z_hat(7,1);
c2=z_hat(8,1);
ugdd1=force(1,1);
ugdd2=force(2,1);
H_value(1,1)=(c2*yd_2)/10 + (k2*y_2)/10 - yd_1*(c1/10 + c2/10) - y_1*(k1/10 + k2/10);
H_value(2,1)=(c2*yd_1)/10 - (c2*yd_2)/10 + (k2*y_1)/10 - (k2*y_2)/10;
end