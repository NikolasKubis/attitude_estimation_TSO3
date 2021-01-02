clear all;
clc;

inertia=diag([6*10^(0);7*10^(0);9*10^(0)]);

dt=0.01;

simulation_time=100;

samples=simulation_time/dt;

measurement_noise =((0.01))*randn(6,samples+1);

mean13=mean(measurement_noise(1:3,:),2);

mean46=mean(measurement_noise(1:3,:),2);

model_uncertainty=(0.)*rand(3,samples+1);

mean=mean(model_uncertainty(1:3,:),2);

ll=1;
 for i=0:dt:simulation_time

tracking_direction1(:,ll)=[0,0,1];

tracking_direction2(:,ll)=[0,1,0];

   model_error(:,ll)=0*[cos(((2*pi)/5)*i);sin(((2*pi)/5)*i);cos(((2*pi)/5)*i)*sin(((2*pi)/5)*i)];

ll=ll+1;
 end
 

 

time=zeros(samples,1);

R=zeros(3,3,samples);

R_est=zeros(3,3,samples);

errorR2=zeros(samples,1);

errorw=zeros(3,samples);

R0=eye(3);

R(:,:,1)=rodrigues(1.8,[1 1 1]);

w0=[0.3 0.2 0.1]';

values=zeros(6,samples);

trace_k=zeros(samples);

R_est0=eye(3);

R_est(:,:,1)=R_est0;

a1 = zeros(3,samples);

a2 = zeros (3,samples);

a1_n= zeros(3,samples);

a2_n= zeros(3,samples);

a1(:,1)=R0'*tracking_direction1(:,1);

a2(:,1)=R0'*tracking_direction2(:,1);

w_est=zeros(3,samples);

w=zeros(3,samples);

w(:,1)=w0;





k=1;
for i=0:dt:simulation_time
    
    
    nominal_input=[sin(((2*pi)/3)*i);cos(((2*pi)/1)*i);sin(((2*pi)/5)*i)];
    
    c_k=alpha_integration(-dt*w(:,k))'*inertia*w(:,k)+dt*(nominal_input+model_error(:,k));

    
    norma=1;
    
    w_n=w(:,k);
    
    iter=0;
    
    while norma>(10^(-12))
        
        w_n_1=w_n-((jacobian_e_integration(inertia,w_n,dt))\e_integration(w_n,c_k,inertia,dt));
        
        norma=norm(e_integration(w_n_1,c_k,inertia,dt));
        
        iter=iter+1;
        
        w_n=w_n_1;
        
    end
    w(:,k+1)=w_n_1;
    
    R(:,:,k+1) = R(:,:,k)*expm(dt*skew(w(:,k)));
    
    a1(:,k)=R(:,:,k)'*tracking_direction1(:,k);
    
    a2(:,k)=R(:,:,k)'*tracking_direction2(:,k);
    
    a1_n(:,k)=a1(:,k)+measurement_noise(1:3,k)-mean13;
    
    a2_n(:,k)=a2(:,k)+measurement_noise(4:6,k)-mean46;
    
    k=k+1;
end

alpha=0.01;

w_est(:,1)=[0;0;0];

k_n=eye(6);

schur=zeros(3,3,samples);

% second connection function

% k=1;
% 
% for i=0:dt:simulation_time
%     
%     k11=k_n(1:3,1:3);
%     
%     k21=k_n(4:6,1:3);
%     
%     k12=k_n(1:3,4:6);
%     
%     k22=k_n(4:6,4:6);
%     
%     
%     nominal_input=[sin(((2*pi)/10)*i);cos(((2*pi)/2)*i);sin(((2*pi)/5)*i)];
%     
%     rR=-(skew(R_est(:,:,k)'*tracking_direction1(:,k))*(a1_n(:,k)))-(skew(R_est(:,:,k)'*tracking_direction2(:,k))*(a2_n(:,k)));
%     
%     c_k2=alpha_integration(-dt*w_est(:,k))'*inertia*w_est(:,k)+dt*(nominal_input+(inertia*k21*rR));
%     
%     norma2=1;
%     
%     w_n2=w_est(:,k);
%     
%     iter2=0;
%     
%     while norma2>(10^(-10))
%         
%         w_n_12=w_n2-((jacobian_e_integration(inertia,w_n2,dt))\e_integration(w_n2,c_k2,inertia,dt));
%         
%         norma2=norm(e_integration(w_n_12,c_k2,inertia,dt),2);
%         
%         iter2=iter2+1;
%         
%         w_n2=w_n_12;
%         
%     end
%     w_est(:,k+1)=w_n2;
%     
%     R_est(:,:,k+1) = R_est(:,:,k)*expm(dt*skew(w_est(:,k)+k11*rR));
%     
%     
%     A=[zeros(3,3) eye(3);zeros(3,3) inv(inertia)*(skew(inertia*w_est(:,k))-skew(w_est(:,k))*inertia)];
%     
%     sum= -(skew(a1_n(:,k))*skew(R_est(:,:,k)'*tracking_direction1(:,k)))-(skew(a2_n(:,k))*skew(R_est(:,:,k)'*tracking_direction2(:,k)));
%     
%     E=[sum zeros(3,3);zeros(3,3) zeros(3,3)];
%     
%     BRB=[zeros(3,3) zeros(3,3);zeros(3,3) eye(3)];
%     
%     W=zeros(6,6);
%     
%     k_n_1=k_n+dt*((-alpha*k_n)+(A*k_n)+(k_n*A')- (k_n*E*k_n)+(BRB)-(W*k_n)-(k_n*W'));
%     
%     k_n=k_n_1;
%     
%     tr=trace(eye(3)-R(:,:,k)'*R_est(:,:,k));
%     
%     eR=acosd(1-(tr/2));
%     
%     errorR(k)=acosd(1-(tr/2));
%     
%     ew=w_est(:,k)-w(:,k);
%     
%     errorw(:,k)=(w_est(:,k)-w(:,k));
%     
%     norm_2_errorw(k)=norm(errorw(:,k),2);
%     
%     schur(:,:,k)=k11-k12*inv(k22)*k21;
%     
%     [v,d]=eig(schur(:,:,k));
%     
%     value(:,k)=diag(d);
%     
%     vectorrr(:,:,k)=v;
%     
%     k=k+1;
%     
%     
% end

alpha2=0.01;

k_n2=eye(6);

k=1;

rrr=1;

w_est2=zeros(3,samples);

R_est2=zeros(3,3,samples);

R_est2(:,:,1)=eye(3);

errorR2=zeros(samples,1);

errorw2=zeros(3,samples);

norm_2_errorw2=zeros(samples);
q1=0.1;
q2=0.1;

for i=0:dt:simulation_time
    
    k112=k_n2(1:3,1:3);
    
    k212=k_n2(4:6,1:3);
    
    nominal_input=[sin(((2*pi)/3)*i);cos(((2*pi)/1)*i);sin(((2*pi)/5)*i)];
    
    rR2=-(skew(R_est2(:,:,k)'*tracking_direction1(:,k))*(a1_n(:,k)))-(skew(R_est2(:,:,k)'*tracking_direction2(:,k))*(a2_n(:,k)));
    
    c_k22=alpha_integration(-dt*w_est2(:,k))'*inertia*w_est2(:,k)+dt*(nominal_input+(inertia*k212*rR2));
    
    norma22=1;
    
    w_n22=w_est2(:,k);
    
    iter22=0;
    
    while norma22>(10^(-10))
        
        w_n_122=w_n22-((jacobian_e_integration(inertia,w_n22,dt))\e_integration(w_n22,c_k22,inertia,dt));
        
        norma22=norm(e_integration(w_n_122,c_k22,inertia,dt),2);
        
        iter22=iter22+1;
        
        w_n22=w_n_122;
        
    end
    w_est2(:,k+1)=w_n22;
    
    R_est2(:,:,k+1) = R_est2(:,:,k)*expm(dt*skew(w_est2(:,k)+k112*rR2));
    
    A2=[-skew(w_est2(:,k)) eye(3);zeros(3,3) inv(inertia)*(skew(inertia*w_est2(:,k))-skew(w_est2(:,k))*inertia)];
    
    sum1=-0.5*q1*(skew(R_est2(:,:,k)'*tracking_direction1(:,k))*skew(a1_n(:,k))+(skew(a1_n(:,k))*skew(R_est2(:,:,k)'*tracking_direction1(:,k))));
    
    sum2=-0.5*q2*(skew(R_est2(:,:,k)'*tracking_direction2(:,k))*skew(a2_n(:,k))+(skew(a2_n(:,k))*skew(R_est2(:,:,k)'*tracking_direction2(:,k))));
    
    sum=sum1+sum2;
    
    E=[sum zeros(3,3);zeros(3,3) zeros(3,3)];
    
    BRB=[zeros(3,3) zeros(3,3);zeros(3,3) 1/rrr*eye(3)];
    
    W=[0.5*skew(k112*rR2) zeros(3,3); zeros(3,3) zeros(3,3)];
    
    k_n_12=k_n2+dt*((-alpha*k_n2)+(A2*k_n2)+(k_n2*A2')- (k_n2*E*k_n2)+(BRB)-(W*k_n2)-(k_n2*W'));
    
    k_n2=k_n_12;
    
    tr=trace(eye(3)-R(:,:,k)'*R_est2(:,:,k));
    
    eR2=acosd(1-(tr/2));
    
    errorR2(k)=acosd(1-(tr/2));
    
    errorw2(:,k)=(w_est2(:,k)-w(:,k));
    
    norm_2_errorw2(k)=norm(errorw2(:,k),2);
    
    
    delta(:,k)=k212*rR2;
    
    n(k)=norm(delta(:,k));
    
    k=k+1;
    
    
end


figure(1)

t=linspace(0,simulation_time,samples);
plot(t,errorR2(1:samples),'-','Linewidth',0.5,'color','blue')
ylabel('e_R (deg)')
set(gca,'TickLabelInterpreter','latex');
set(gca,'fontweight','bold','fontsize',40)
ylabel('$e_R (deg)$','interpreter','latex')
xlabel('$ time (s)$','interpreter','latex')
set(gca,'FontSize',40)
grid on


figure(6)

subplot(3,1,1)
t=linspace(0,simulation_time,samples);
plot(t,errorw2(1,1:samples).*57.2957795,'-','Linewidth',0.9,'color','blue')
set(gca,'TickLabelInterpreter','latex');
set(gca,'fontweight','bold','fontsize',40)
ylabel('$e^x_W \ (deg/s)$','interpreter','latex')
set(gca,'FontSize',40)
grid on

subplot(3,1,2)
t=linspace(0,simulation_time,samples);
plot(t,errorw2(2,1:samples).*57.2957795,'-','Linewidth',0.9,'color','blue')
grid on
set(gca,'TickLabelInterpreter','latex');
set(gca,'fontweight','bold','fontsize',40)
ylabel('$e^y_W \  (deg/s)$','interpreter','latex')

set(gca,'FontSize',40)
grid on
subplot(3,1,3)
t=linspace(0,simulation_time,samples);
plot(t,errorw2(3,1:samples).*57.2957795,'-','Linewidth',0.9,'color','blue')
set(gca,'TickLabelInterpreter','latex');
set(gca,'fontweight','bold','fontsize',40)
ylabel('$e^z_W \ (deg/s)$','interpreter','latex')
set(gca,'FontSize',40)
xlabel('$ time (s)$','interpreter','latex')
grid on
% 
% 
% figure(8)
% t=linspace(0,simulation_time,samples);
% plot(t,delta(1,1:samples),'-','Linewidth',1,'color','blue')
% hold on;
% plot(t,model_error(1,1:samples),'-','Linewidth',1,'color','red')
% set(gca,'TickLabelInterpreter','latex');
% set(gca,'fontweight','bold','fontsize',40) 
% ylabel('$model\ error$','interpreter','latex')
% set(gca,'FontSize',40)
% xlabel('$ time (s)$','interpreter','latex')
% ylim([-1.5,1.5])
% grid on
% 
% figure(9)
% t=linspace(0,simulation_time,samples);
% plot(t,delta(2,1:samples),'-','Linewidth',1,'color','blue')
% hold on;
% plot(t,model_error(2,1:samples),'-','Linewidth',1,'color','red')
% set(gca,'TickLabelInterpreter','latex');
% set(gca,'fontweight','bold','fontsize',40) 
% ylabel('$model\ error$','interpreter','latex')
% set(gca,'FontSize',40)
% xlabel('$ time (s)$','interpreter','latex')
% ylim([-1.5,1.5])
% grid on
% 
% 
% figure(10)
% t=linspace(0,simulation_time,samples);
% plot(t,delta(3,1:samples),'-','Linewidth',1,'color','blue')
% hold on;
% plot(t,model_error(3,1:samples),'-','Linewidth',1,'color','red')
% set(gca,'TickLabelInterpreter','latex');
% set(gca,'fontweight','bold','fontsize',40) 
% ylabel('$model\ error$','interpreter','latex')
% set(gca,'FontSize',40)
% xlabel('$ time (s)$','interpreter','latex')
% ylim([-1.5,1.5])
% grid on
% 
% figure(11)
% t=linspace(0,simulation_time,samples);
% plot(t,n(1,1:samples).*57.2957795,'-','Linewidth',1,'color','green')
% set(gca,'TickLabelInterpreter','latex');
% set(gca,'fontweight','bold','fontsize',40) 
% ylabel('$||\delta(t)||$','interpreter','latex')
% set(gca,'FontSize',40)
% xlabel('$ time (s)$','interpreter','latex')
% grid on


