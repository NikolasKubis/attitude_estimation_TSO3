function q_inv= q_ant(q)

s=q(1);

v=q(2:4);

q_inv=[s  -v']';


end