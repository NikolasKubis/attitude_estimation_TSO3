function mq= q_mult(q1,q2)

s1=q1(1);
s2=q2(1);
v1=q1(2:4);
v2=q2(2:4);

mq=[s1*s2-(v1'*v2);(s1*v2+s2*v1)+(skew(v1)*v2)];


end