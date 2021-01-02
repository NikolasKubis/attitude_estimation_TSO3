function prediction_sigma=prediction(sp,nominal_input)

prediction_sigma=zeros(7,size(sp,2));

for l=1:size(sp,2)
    
    qk=sp(1:4,l);
    
    wk=sp(5:7,l);
    
    model_errork=sp(8:10,l);
    
    w_imp=wk+((dt/2)*inv(inertia)*((skew(inertia*wk)*wk)+nominal_input+ model_errork));
    
    omega=0.5*[0 -w_imp'; w_imp -skew(w_imp)];
    
    qk1=expm(dt*omega)*qk;
    
    norma=100;
    
    w_n=w_imp;
    
    while norma>(10^(-10))
        
        fw_n= -w_n+w_imp+(dt/2)*inv(inertia)*(skew(inertia*w_n)*w_n+nominal_input+model_errork);
        
        jw_n= -eye(3)+((dt/2)*inv(inertia)*(skew(inertia*w_n)-skew(w_n)*inertia));
        
        w_n_1=w_n-inv(jw_n)*fw_n;
        
        w_n=w_n_1;
        
        norma=norm(fw_n,1);
    end
    wk1=w_n_1;
    
  prediction_sigma(1:7,l)=[qk1;wk1];
    
end


end