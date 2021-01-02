function output_point=out(predicted_sigma_points,tracking_direction1,tracking_direction2)


q=predicted_sigma_points(1:4);

out_1=output_matrix(q)'*tracking_direction1;

out_2=output_matrix(q)'*tracking_direction2;

output_point=[out_1;out_2];


end