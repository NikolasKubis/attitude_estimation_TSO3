function cross_covariance_matrix=cross(weight_c,predicted_sigma_points,output_sigma_points,...
    predicted_state,predicted_output,L)


sum=zeros(7,7);

for k=1:2*L
    
    term=(predicted_sigma_points(:,k)-predicted_state)*([output_sigma_points(:,k);0]-[predicted_output;0])';
    
    sum=sum+weight_c(k)*term;
    
    
end
cross_covariance_matrix=sum(:,1:6);

end