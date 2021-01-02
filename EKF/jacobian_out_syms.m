syms   q1 q2 q3 q4 tracking_direction1 t11 t21 t31 tracking_direction2 t12 t22 t32 real;

q=[q1;q2;q3;q4];

tracking_direction1=[t11 t21 t31]';

tracking_direction2=[t12 t22 t32]';

D=[q1^2+q2^2-q3^2-q4^2  2*(q2*q3+q1*q4)  2*(q2*q4-q1*q3);...
         2*(q2*q3-q1*q4) q1^2-q2^2+q3^2-q4^2  2*(q4*q3+q1*q2);...
         2*(q2*q4+q1*q3) 2*(q4*q3-q1*q2) q1^2-q2^2-q3^2+q4^2];

h1=D'*tracking_direction1;

h2=D'*tracking_direction2;

h=[h1;h2];

v=[q1;q2;q3;q4];
J=jacobian(h,v)