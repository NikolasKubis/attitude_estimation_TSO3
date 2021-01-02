function  output_covariance=covy(weight_c,output_sigma_points,predicted_output,output_noise_covariance,L)

sum=zeros(6,6);

for k=1:2*L
    
    term=weight_c(k)*(output_sigma_points(:,k)-predicted_output)*(output_sigma_points(:,k)-predicted_output)';
    
    sum=sum+term;
    
end

output_covariance=sum+output_noise_covariance;


end
