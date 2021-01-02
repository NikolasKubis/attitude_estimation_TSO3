
clear all;
clc;


% This script implements a geometric predictive filter on the tangent
% bundle of the special orthogonal group TSO(3). The error signal is formed
% by the cross product rR (see line 153). For more details see 
%(http://resolver.tudelft.nl/uuid:23eac88f-a75a-4d8d-aaa5-dcfe32d7c909)
% Necessary functions are rodrigues.m and skew.m. The former one is
% utilized for initializing the orientation as it provides more intuition
% for the starting attitude point. The latter implements the isomorphism
% from the R^3 to the lie algebra so(3). 

inertia=diag([60*10^(0);70*10^(0);90*10^(0)]);

dt=0.01;

simulation_time=100;

samples=simulation_time/dt;

errorw=zeros(3,samples);

tracking_direction1=zeros(3,samples);

tracking_direction2=zeros(3,samples);

ll=1;
for i=0:dt:simulation_time
    
    tracking_direction1(:,ll)=[0;1;0];
    tracking_direction2(:,ll)=[0;0;1];
    
    % model_error(:,ll)=0.0*[cos(((2*pi)/5)*i);sin(((2*pi)/5)*i);cos(((2*pi)/5)*i)*sin(((2*pi)/5)*i)];
    model_error(:,ll)=0.0*[cos(((2*pi)/5)*i);sin(((2*pi)/5)*i);cos(((2*pi)/5)*i)*sin(((2*pi)/5)*i)];
    
    % nm(ll)=norm(model_error(:,ll));
    
    
    ll=ll+1;
end

errorR=zeros(samples,1);

w_est=zeros(3,samples);

R_est=zeros(3,3,samples);

R_est(:,:,1)=eye(3);

delta=zeros(3,samples);

sigma_delta=1;

Q=diag(10*10^(3)*[1 1 1]);

RRR=diag(4.0*10^(-3)*[1 1 1]);

w_est(:,1)=[0 0 0]';

w=zeros(3,samples);

R=zeros(3,3,samples);

y1=zeros(3,samples);

y2=zeros(3,samples);

w(:,1)=[0.3 0.2 0.1]';

R(:,:,1)=rodrigues(1.7,[1 1 1]);

ww_n=zeros(3,samples);

k=1;

sum1=0;

sum2=0;

index=zeros(1,samples);

for i=0:dt:simulation_time
    
    y1(:,k)=R(:,:,k)'*tracking_direction1(:,k)+((0.001))*randn(3,1);
    
    y2(:,k)=R(:,:,k)'*tracking_direction2(:,k)+((0.001))*randn(3,1);
    
    y(:,k)=[y1(:,k);y2(:,k)];
    
    nominal_input=[sin(((2*pi)/3)*i);cos(((2*pi)/1)*i);sin(((2*pi)/5)*i)];
    
    w_imp=w(:,k)+((dt/2)*inv(inertia)*((skew(inertia*w(:,k))*w(:,k)+nominal_input+inertia*model_error(:,k))));
    
    R(:,:,k+1) = R(:,:,k)*expm(dt*skew(w_imp));
    
    norma=100;
    
    w_n=w_imp;
    
    ww_n(:,k)=w_n;
    
    while norma>(10^(-10))
        
        fw_n=-w_n+w_imp+(dt/2)*inv(inertia)*(skew(inertia*w_n)*w_n+nominal_input+inertia*model_error(:,k));
        
        jw_n= -eye(3)+((dt/2)*inv(inertia)*(skew(inertia*w_n)-skew(w_n)*inertia));
        
        w_n_1=w_n-inv(jw_n)*fw_n;
        
        w_n=w_n_1;
        
        norma=norm(fw_n,1);
    end
    w(:,k+1)=w_n_1;
    
    k=k+1;
end

y1_est=zeros(3,samples);

y2_est=zeros(3,samples);

ww_n2=zeros(3,samples);

sum=zeros(6,6);

k=1;

for i=0:dt:simulation_time
    
    y1_est(:,k)=R_est(:,:,k)'*tracking_direction1(:,k);
    
    y2_est(:,k)=R_est(:,:,k)'*tracking_direction2(:,k);
    
    y_est(:,k)=[y1_est(:,k);y2_est(:,k)];
    
    nominal_input=[sin(((2*pi)/3)*i);cos(((2*pi)/1)*i);sin(((2*pi)/5)*i)];
    
    A=[-(skew(w_est(:,k))*R_est(:,:,k)'*tracking_direction1(:,k))*dt; -(skew(w_est(:,k))*R_est(:,:,k)'*tracking_direction2(:,k))*dt];
    
    G=[(skew(R_est(:,:,k)'*tracking_direction1(:,k))*inv(inertia)*((skew(inertia*w_est(:,k))*w_est(:,k))+nominal_input)+...
        ((skew(w_est(:,k))^2)*R_est(:,:,k)'*tracking_direction1(:,k)))*(((dt)^2)/2);...
        (skew(R_est(:,:,k)'*tracking_direction2(:,k))*inv(inertia)*((skew(inertia*w_est(:,k))*w_est(:,k))+nominal_input)+...
        ((skew(w_est(:,k))^2)*R_est(:,:,k)'*tracking_direction2(:,k)))*(((dt)^2)/2)];
    
    zetta(:,k)=A+G;
    
    zeta1(:,k)=zetta(1:3,k);
    
    zeta2(:,k)=zetta(4:6,k);
    
    B=sigma_delta*eye(3);
    
    Lambda=(((dt)^2)/2);

    rR=-(skew(R_est(:,:,k)'*tracking_direction1(:,k))*(y1(:,k)))-(skew(R_est(:,:,k)'*tracking_direction2(:,k))*(y2(:,k)));
    
    W=[skew(R_est(:,:,k)'*tracking_direction1(:,k))*B;skew(R_est(:,:,k)'*tracking_direction2(:,k))*B];
    
    W1=skew(R_est(:,:,k)'*tracking_direction1(:,k))*B;
    
    W2=skew(R_est(:,:,k)'*tracking_direction2(:,k))*B;
    
    betta=skew(y1(:,k))*Lambda*W1+skew(y2(:,k))*Lambda*W2;
    
    alpha=skew(y1(:,k))*(R_est(:,:,k)'*tracking_direction1(:,k))+skew(y2(:,k))*(R_est(:,:,k)'*tracking_direction2(:,k))+...
        skew(y1(:,k))*zeta1(:,k)+skew(y2(:,k))*zeta2(:,k);
    
    delta(:,k)=-0.5*inv(betta'*Q'*betta+RRR')*betta'*(2*Q)*alpha;
    
    w_imp2=w_est(:,k)+((dt/2)*inv(inertia)*((skew(inertia*w_est(:,k))*w_est(:,k))+nominal_input+inertia*delta(:,k)));
    
    R_est(:,:,k+1)=R_est(:,:,k)*expm(dt*skew(w_imp2+0.7*rR));
    
    norma=100;
    
    w_n=w_imp2;
    
    ww_n2(:,k)= w_n;
    
    while norma>(10^(-5))
        
        fw_n2=-w_n+w_imp2+(dt/2)*inv(inertia)*((skew(inertia*w_n)*w_n)+(nominal_input)+delta(:,k));
        
        jw_n= -eye(3)+((dt/2)*inv(inertia)*(skew(inertia*w_n)-skew(w_n)*inertia));
        
        w_n_1=w_n-inv(jw_n)*fw_n2;
        
        w_n=w_n_1;
        
        norma=norm(fw_n2,1);
    end
    
    w_est(:,k+1)=w_n_1;
    
    errorw(:,k)=(w(:,k)-w_est(:,k));
    
    tr=trace(eye(3)-R(:,:,k)'*R_est(:,:,k));
    
    errorR(k)=acosd(1-(tr/2));
    
    nd(k)=norm(delta(:,k));
    
    errory= (y_est(:,k)-y(:,k))*(y_est(:,k)-y(:,k))';
    
    sum=sum+errory;
    
    k=k+1;
end

sum=sum/samples;



figure(1)

t=linspace(0,simulation_time,samples);
plot(t,movmean(errorR(1:samples),10),'-','Linewidth',0.5,'color','blue')
set(gca,'TickLabelInterpreter','latex');
set(gca,'fontweight','bold','fontsize',40)
ylabel('$e_R (deg)$','interpreter','latex')
xlabel('$ time (s)$','interpreter','latex')
set(gca,'FontSize',40)
grid on


figure(2)

subplot(3,1,1)
t=linspace(0,simulation_time,samples);
plot(t,movmean(errorw(1,1:samples),2).*57.2957795,'-','Linewidth',1,'color','black')
set(gca,'TickLabelInterpreter','latex');
set(gca,'fontweight','bold','fontsize',40)
ylabel('$e^x_{\Omega} \ (deg/s)$','interpreter','latex')
set(gca,'FontSize',40)
ylim([-50,50])
grid on


subplot(3,1,2)
t=linspace(0,simulation_time,samples);
plot(t,movmean(errorw(2,1:samples),2).*57.2957795,'-','Linewidth',1,'color','black')
set(gca,'TickLabelInterpreter','latex');
set(gca,'fontweight','bold','fontsize',40)
ylabel('$e^y_{\Omega} \ (deg/s)$','interpreter','latex')
set(gca,'FontSize',40)
ylim([-50,50])
grid on

subplot(3,1,3)
t=linspace(0,simulation_time,samples);
plot(t,movmean(errorw(3,1:samples),2).*57.2957795,'-','Linewidth',1,'color','black')
set(gca,'TickLabelInterpreter','latex');
set(gca,'fontweight','bold','fontsize',40)
ylabel('$e^z_{\Omega} \ (deg/s)$','interpreter','latex')
set(gca,'FontSize',40)
xlabel('$ time (s)$','interpreter','latex')
ylim([-50,50])
grid on
