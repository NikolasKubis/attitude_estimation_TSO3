function predicted_output=prediction_output(output_sigma_points,weight_m,L)

sum=zeros(6,1);

for k=1:2*L
    
    term=weight_m(k)*output_sigma_points(:,k);
    
    sum=sum+term;
    
end


predicted_output=sum;



end
