function F = j_sys_EKF_att(x,inertia)

q1=x(1);
q2=x(2);
q3=x(3);
q4=x(4);
w1=x(5);
w2=x(6);
w3=x(7);
I1=inertia(1,1);
I2=inertia(2,2);
I3=inertia(3,3);

F=[  0 -w1/2 -w2/2 -w3/2 -q2/2 -q3/2 -q4/2;
    w1/2 0 w3/2 -w2/2 q1/2 -q4/2 q3/2;
    w2/2 -w3/2 0 w1/2 q4/2 q1/2 -q2/2;
    w3/2 w2/2 -w1/2 0 -q3/2 q2/2 q1/2;
    0 0 0 0  0 (I2*w3 - I3*w3)/I1  (I2*w2 - I3*w2)/I1;
    0 0 0 0 -(I1*w3 - I3*w3)/I2 0 -(I1*w1 - I3*w1)/I2;
    0 0 0 0 (I1*w2 - I2*w2)/I3 (I1*w1 - I2*w1)/I3 0];

end