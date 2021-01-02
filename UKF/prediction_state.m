function predicted_state=prediction_state(predicted_sigma_points,weight_m,L)

sum=zeros(7,1);

for k=1:2*L
    
    term=weight_m(k)*predicted_sigma_points(:,k);
    
    sum=sum+term;
    
end

sum(1:4)=sum(1:4)/norm(sum(1:4));

predicted_state=sum;

end