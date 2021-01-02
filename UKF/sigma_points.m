function [point,weight_c,weight_m]=sigma_points(x,Px,L,process_noise_covariance)

square_root=chol(Px+process_noise_covariance);
point=zeros(7,2*L);

weight_c=zeros(2*L,1);

weight_m=zeros(2*L,1);

kapa=0.01;

alpha =0.9; % ayti einai i simantikiparametros. Apo edoelegxeis to pos tha katanemontai ta sigma points
            % An i avevaiothta poy exeis einai mikri tote de xreiazetai na
            % aploseis poli ta sigma points. An i avevaiothta einai megali
            % prepei na ta aploseis polukai etsi to alpha prepei na to
            % pareis megalo.

betta= 0.1;

lambda=(alpha^2)*(L+kapa)-L; % prosexe oti to lambda pollaplasiazei ti stili tou covariance matrix
                             % diladi tis kateythinsis avevaiotitas-kai
                             % etsi ftiaxneis tin katanomi. An exeis pola
                             % shmeia tote prepei na pareis magalyteres
                             % apostaseis gia na einai katanemimena ta
                             % s-points omoiomorfa. Etsi loipon to lambda
                             % ayxanei me to dimension. Epeita, analoga me
                             % tin avevaiotita poy exeis dialegeis to alpha
                             

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