function  state_covariance=covx(weight_c,predicted_sigma_points,predicted_state,process_noise_covariance,L)

sum=zeros(7,7);

for k=1:2*L
    
   
    term=weight_c(k)*(predicted_sigma_points(:,k)-predicted_state)*(predicted_sigma_points(:,k)-predicted_state)';
    
    sum=sum+term;
    
end

state_covariance=sum;%+process_noise_covariance;


end
