syms  w w1 w2 w3 OmegaM skew skew2 q q1 q2 q3 q4 I1 I2 I3 Inertia T t1 t2 t3 f1 f2 f real;

q=[q1;q2;q3;q4];

T =[t1;t2;t3];

w=[w1;w2;w3];

skew =[0 -w3 w2;w3 0 -w1;-w2 w1 0];

skew2 =[0 -I3*w3 I2*w2;I3*w3 0 -I1*w1;-I2*w2 I1*w1 0];

OmegaM=0.5*[0 -w'; w -skew];

Inertia=diag([I1, I2, I3]);

f1=OmegaM*q;

f2=inv(Inertia)*((skew2*w)+T);

f=[f1;f2];

v=[q1 q2 q3 q4 w1 w2 w3];

J=jacobian(f,v);





