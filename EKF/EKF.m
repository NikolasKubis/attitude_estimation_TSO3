% The Extended Kalman Filter.

clear all;

% inertia tensor
inertia=diag([6*10^(0);7*10^(0);9*10^(0)]);

dt=0.01;

simulation_time=100;

samples=simulation_time/dt;


% establishing a deterministic model error
ll=1;
for i=0:dt:simulation_time
    
    model_error(:,ll)=0.0*[cos(((2*pi)/5)*i);sin(((2*pi)/5)*i);cos(((2*pi)/5)*i)*sin(((2*pi)/5)*i)];
    
    ll=ll+1;
end

% tracking directions are chosen to be constant. However, the filter works
% with time varying reference directions as well.
tracking_direction1=[0 0 1]';

tracking_direction2=[0 1 0]';


% motion initialization
w=zeros(3,samples);

q=zeros(4,samples);

theta=1.8;

l=[1 1 1];
ll=l/norm(l);
l1=ll(1);
l2=ll(2);
l3=ll(3);

q(:,1)=[cos(theta/2) l1*sin(theta/2) l2*sin(theta/2) l3*sin(theta/2)]; % initial quaternion

y1=zeros(3,samples);

y2=zeros(3,samples);

% initial angular rate in rad/sec
w(:,1)=[0.3 0.2 0.1]';

y=zeros(6,samples);

d=0.01; % measurement noise variance (the two sensors are considered uncorelated)
k=1;

for i=0:dt:simulation_time
    
    nominal_input=[sin(((2*pi)/3)*i);cos(((2*pi)/1)*i);sin(((2*pi)/5)*i)];
    
    % Lie-verlet integration scheme
    w_imp=w(:,k)+((dt/2)*(inv(inertia)*((skew(inertia*w(:,k))*w(:,k))+nominal_input)+inertia*model_error(:,k)));
    
    omega=0.5*[0 -w_imp'; w_imp -skew(w_imp)];
    
    q(:,k+1)=expm(dt*omega)*q(:,k);
    
    
    % newton solver
    norma=100;
    
    w_n=w_imp;
    
    while norma>(10^(-10))
        
        fw_n=-w_n+w_imp+(dt/2)*inv(inertia)*(skew(inertia*w_n)*w_n+nominal_input+inertia*model_error(:,k));
        
        jw_n= -eye(3)+((dt/2)*inv(inertia)*(skew(inertia*w_n)-skew(w_n)*inertia));
        
        w_n_1=w_n-inv(jw_n)*fw_n;
        
        w_n=w_n_1;
        
        norma=norm(fw_n,1);
    end
    w(:,k+1)=w_n_1;
    
    
    y1(:,k)=output_matrix(q(:,k))'*tracking_direction1+((d))*randn(3,1);
    
    
    y2(:,k)=output_matrix(q(:,k))'*tracking_direction2+((d))*randn(3,1);
    
    
    y(:,k)=[y1(:,k);y2(:,k)];
    
    k=k+1;
end

% filter's initialization

% //model error covariance (trial and error tunning)
r=5;

R1=r*eye(3);

R2=r*eye(3);
%//



Q=0.01*eye(7,7);

R=[R1 zeros(3,3);zeros(3,3) R2];

y1_est=zeros(3,samples);

y2_est=zeros(3,samples);

y_est=zeros(6,samples);

PPY=zeros(6,6,samples);

P=zeros(7,7,samples);

PX=zeros(7,6,samples);

K=zeros(7,6,samples);

prediction=zeros(7,samples);

correction=zeros(7,samples);

correction(1:4,1)=[1 0 0 0]';

correction(5:7,1)=[0 0 0]';

y1_est(:,1)=output_matrix(correction(1:4,1))'*tracking_direction1;

y2_est(:,1)=output_matrix(correction(1:4,1))'*tracking_direction2;

y_est(:,1)=[y1_est(:,1);y2_est(:,1)];

w_est=zeros(3,samples);

errorR=zeros(samples);

errorw=zeros(3,samples);

P(:,:,1)=0.01*eye(7);

k=1;
for i=0:dt:simulation_time
    
    q_est(:,k)=correction(1:4,k);
    
    w_est(:,k)=correction(5:7,k);
    
    nominal_input=[sin(((2*pi)/3)*i);cos(((2*pi)/1)*i);sin(((2*pi)/5)*i)];
    
    w_imp2=w_est(:,k)+((dt/2)*inv(inertia)*((skew(inertia*w_est(:,k))*w_est(:,k))+nominal_input));
    
    omega2=0.5*[0 -w_imp2'; w_imp2 -skew(w_imp2)];
    
    q_est(:,k+1)=expm(dt*omega2)*(q_est(:,k));
    
    norma=100;
    
    w_n2=w_imp2;
    
    while norma>(10^(-10))
        
        fw_n2=-w_n2+w_imp2+(dt/2)*inv(inertia)*(skew(inertia*w_n2)*w_n2+nominal_input);
        
        jw_n2= -eye(3)+((dt/2)*inv(inertia)*(skew(inertia*w_n2)-skew(w_n2)*inertia));
        
        w_n_12=w_n2-inv(jw_n2)*fw_n2;
        
        w_n2=w_n_12;
        
        norma=norm(fw_n2,1);
    end
    
    w_est(:,k+1)=w_n_12;
    
    prediction(:,k)=[q_est(:,k+1);w_est(:,k+1)];
    
    P(:,:,k)=(j_sys_EKF_att(prediction(:,k),inertia)*P(:,:,k))+(P(:,:,k)*j_sys_EKF_att(prediction(:,k),inertia)')+Q;
    
    y1_est(:,k)=output_matrix(q_est(:,k+1))'*tracking_direction1;
    
    y2_est(:,k)=output_matrix(q_est(:,k+1))'*tracking_direction2;
    
    y_est(:,k)=[y1_est(:,k);y2_est(:,k)];
    
    PPY(:,:,k)=j_out_EKF_att(q_est(:,k),tracking_direction1,tracking_direction2)*P(:,:,k)*j_out_EKF_att(q_est(:,k),tracking_direction1,tracking_direction2)'+R;
    
    PX(:,:,k)=P(:,:,k)*j_out_EKF_att(q_est(:,k+1),tracking_direction1,tracking_direction2)';
    
    K(:,:,k)=PX(:,:,k)/(PPY(:,:,k));
    
    correction(:,k+1)=prediction(:,k) + K(:,:,k)*(y(:,k)-y_est(:,k));
    
    P(:,:,k+1)=P(:,:,k)-K(:,:,k)*PPY(:,:,k)*K(:,:,k)';
    
    correction(1:4,k+1)=correction(1:4,k+1)/norm(correction(1:4,k+1));
    
    errorw(:,k)=(w_est(:,k)-w(:,k));
    
    DD=output_matrix(q(:,k));
    
    DD_est=output_matrix(q_est(:,k));
    
    tr=trace(eye(3)-(DD'*DD_est));
    
    errorR(k)=acosd(1-(tr/2));
    
    k=k+1;
end


figure(1)

t=linspace(0,simulation_time,samples);
plot(t,errorR(1:samples),'-','Linewidth',0.5,'color','red')
set(gca,'TickLabelInterpreter','latex');
set(gca,'fontweight','bold','fontsize',16)
ylabel('$e_R (deg)$','interpreter','latex')
xlabel('$ time (sec)$','interpreter','latex')
set(gca,'FontSize',16)
grid on





figure(2)
subplot(3,1,1)
t=linspace(0,simulation_time,samples);
plot(t,errorw(1,1:samples).*57.2957795,'-','Linewidth',0.9,'color','red')
hold on
set(gca,'TickLabelInterpreter','latex');
set(gca,'fontweight','bold','fontsize',16)
ylabel('$e_{\omega}^x \ (deg/sec)$','interpreter','latex')
set(gca,'FontSize',16)
grid on


subplot(3,1,2)
t=linspace(0,simulation_time,samples);
plot(t,errorw(2,1:samples).*57.2957795,'-','Linewidth',0.9,'color','red')
hold on
set(gca,'TickLabelInterpreter','latex');
set(gca,'fontweight','bold','fontsize',16)
ylabel('$e_{\omega}^y \ (deg/sec)$','interpreter','latex')
set(gca,'FontSize',16)
grid on

subplot(3,1,3)
t=linspace(0,simulation_time,samples);
plot(t,errorw(3,1:samples).*57.2957795,'-','Linewidth',0.9,'color','red')
hold on
set(gca,'TickLabelInterpreter','latex');
set(gca,'fontweight','bold','fontsize',16)
xlabel('$ time (sec)$','interpreter','latex')
ylabel('$e_{\omega}^z \ (deg/sec)$','interpreter','latex')
set(gca,'FontSize',16)
grid on



