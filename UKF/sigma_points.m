function [point,weight_c,weight_m]=sigma_points(x,Px,L,process_noise_covariance)

square_root=chol(Px+process_noise_covariance);
point=zeros(7,2*L);

weight_c=zeros(2*L,1);

weight_m=zeros(2*L,1);

kapa=0.01;

alpha =0.9; 

betta= 0.1;

lambda=(alpha^2)*(L+kapa)-L; 
                             

weight_m(1)=lambda/(L+lambda);

weight_c(1)=(lambda/(L+lambda))+(1-alpha^2+betta);

point(:,1)=x;

ll=1;
for k=1:L
    
    point(:,ll)=x+(sqrt(L+lambda)*square_root(:,k));
    
    point(1:4,ll)=point(1:4,ll)/norm(point(1:4,ll));
    
    point(:,ll+L)=x-(sqrt(L+lambda)*square_root(:,k));
    
    point(1:4,ll+L)=point(1:4,ll+L)/norm(point(1:4,ll+L));
    
    ll=ll+1;
end


for k=2:2*L
    weight_c(k)=1/(2*(L+lambda));
    
    weight_m(k)=1/(2*(L+lambda));
end






end