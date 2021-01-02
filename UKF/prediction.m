function  predicted_sigma_point=prediction(point,nominal_input,inertia,dt)


qk=point(1:4);

wk=point(5:7);

w_imp2=wk+((dt/2)*inv(inertia)*((skew(inertia*wk)*wk)+nominal_input));

omega2=0.5*[0 -w_imp2'; w_imp2 -skew(w_imp2)];

qk1=expm(dt*omega2)*qk;

norma=100;

w_n2=w_imp2;

while norma>(10^(-4))
    
    fw_n2=-w_n2+w_imp2+(dt/2)*inv(inertia)*(skew(inertia*w_n2)*w_n2+nominal_input);
    
    jw_n2= -eye(3)+((dt/2)*inv(inertia)*(skew(inertia*w_n2)-skew(w_n2)*inertia));
    
    w_n_12=w_n2-inv(jw_n2)*fw_n2;
    
    w_n2=w_n_12;
    
    norma=norm(fw_n2,1);
end

wk1=w_n_12;

predicted_sigma_point(1:4)=qk1;

predicted_sigma_point(5:7)=wk1;



end