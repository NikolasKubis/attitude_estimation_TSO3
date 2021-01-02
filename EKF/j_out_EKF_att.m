function H = j_out_EKF_att(x,tracking_direction1,tracking_direction2)

q1=x(1);
q2=x(2);
q3=x(3);
q4=x(4);

t11=tracking_direction1(1);
t21=tracking_direction1(2);
t31=tracking_direction1(3);



t12=tracking_direction2(1);
t22=tracking_direction2(2);
t32=tracking_direction2(3);



H=[ 2*q1*t11 - 2*q4*t21 + 2*q3*t31, 2*q2*t11 + 2*q3*t21 + 2*q4*t31, 2*q2*t21 - 2*q3*t11 + 2*q1*t31, 2*q2*t31 - 2*q1*t21 - 2*q4*t11;
 2*q4*t11 + 2*q1*t21 - 2*q2*t31, 2*q3*t11 - 2*q2*t21 - 2*q1*t31, 2*q2*t11 + 2*q3*t21 + 2*q4*t31, 2*q1*t11 - 2*q4*t21 + 2*q3*t31;
 2*q2*t21 - 2*q3*t11 + 2*q1*t31, 2*q4*t11 + 2*q1*t21 - 2*q2*t31, 2*q4*t21 - 2*q1*t11 - 2*q3*t31, 2*q2*t11 + 2*q3*t21 + 2*q4*t31;
 2*q1*t12 - 2*q4*t22 + 2*q3*t32, 2*q2*t12 + 2*q3*t22 + 2*q4*t32, 2*q2*t22 - 2*q3*t12 + 2*q1*t32, 2*q2*t32 - 2*q1*t22 - 2*q4*t12;
 2*q4*t12 + 2*q1*t22 - 2*q2*t32, 2*q3*t12 - 2*q2*t22 - 2*q1*t32, 2*q2*t12 + 2*q3*t22 + 2*q4*t32, 2*q1*t12 - 2*q4*t22 + 2*q3*t32;
 2*q2*t22 - 2*q3*t12 + 2*q1*t32, 2*q4*t12 + 2*q1*t22 - 2*q2*t32, 2*q4*t22 - 2*q1*t12 - 2*q3*t32, 2*q2*t12 + 2*q3*t22 + 2*q4*t32];
 
H=[H zeros(6,3)];

end